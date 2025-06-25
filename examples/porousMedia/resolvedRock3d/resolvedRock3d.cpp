/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Jan E. Marquardt, Mathias J. Krause
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

// Experimental: Setup hasn't been exntensively validated.

/* How to:
 *
 * 1. Provide or download a medium.
 *    For example: https://www.digitalrocksportal.org/projects/92
 * 2. Convert the TIFF series to VTI using ParaView
 *    (Point data are mandatory, cell data do not work)
 *     - Open TIFF series in ParaView
 *     - Apply "cell data to point data" filter (if it's grayed out, applying the filter isn't necessary)
 *     - Save the result as VTI
 *     - Use ASCII formatting
 *     - Solely export one array and use the name as an input when running the simulation (see <arrayname>)
 * 3. Start the simulation with the following parameters (in that order):
 *    <filename> <arrayname> <scaling-factor> <time-scaling-factor> <resolution> <pressure_drop>
 *    Example: mpirun -np 6 ./resolvedRock3d rock.vti "Tiff Scalars" 2.5e-6 1.0 200 1.0
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T          = FLOATING_POINT_TYPE;
using DESCRIPTOR = D3Q19<>;

// Parameters for the simulation setup
int          N; // resolution of the model
constexpr T  tau = 0.75;
T            maxPhysT; // max. simulation time in s, SI unit
T            pressureDrop;
T            scalingFactor;
Vector<T, 3> extent;

// Fluid parameters
constexpr T fluidDensity       = 1000;
constexpr T kinematicViscosity = 1e-6;

// Rock geometry
std::shared_ptr<BlockVTIreader3D<T, T>> vtiReader; // to ensure that the data persists

// Results
Vector<T, 3> currentAverageVelocity;
T            currentPermeability;

std::shared_ptr<IndicatorBlockData3D<T>>
generateIndicatorFromVTI(const std::string vtiFile, const std::string arrayName)
{
  const T sourceScale = scalingFactor;

  vtiReader = std::make_shared<BlockVTIreader3D<T, T>>(vtiFile, arrayName);

  auto           cuboidSample     = vtiReader->getCuboid();
  T              deltaRsample     = cuboidSample.getDeltaR() * sourceScale;
  Vector<int, 3> extentSample     = cuboidSample.getExtent();
  Vector<T, 3>   originSamplePhys = cuboidSample.getOrigin() * sourceScale;
  Vector<T, 3>   extentSamplePhys = {deltaRsample * T(extentSample[0] + 0.5),
                                     deltaRsample * T(extentSample[1] + 0.5),
                                     deltaRsample * T(extentSample[2] + 0.5)};
  for (unsigned i = 0; i < 3; ++i) {
    extent[i] = extentSamplePhys[i];
  }


  std::shared_ptr<IndicatorBlockData3D<T>> ind(
      new IndicatorBlockData3D<T>(vtiReader->getBlockData(), extentSamplePhys,
                                  originSamplePhys, deltaRsample, false));

  return ind;
}

void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter,
                     IndicatorF3D<T>&                    layer,
                     SuperGeometry<T, 3>&                superGeometry)
{

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  // Set the cells inside of the layer indicator to MN 2
  superGeometry.rename(0, 2, layer);

  // Set material number 1 for cells in the rock cavity with 1 cell offset in each direction
  superGeometry.rename(2, 1, {1,1,1});

  Vector<T, 3> origin = superGeometry.getStatistics().getMinPhysR(2);
  Vector<T, 3> extend = superGeometry.getStatistics().getMaxPhysR(2);

  // Set material number for inflow
  origin[0] = superGeometry.getStatistics().getMinPhysR(2)[0] -
              T {0.5} * converter.getPhysDeltaX();
  extend[0] = T {1} * converter.getPhysDeltaX();
  IndicatorCuboid3D<T> inflow(extend, origin);
  // MN 3 is set only if the neighbor cell has MN 1
  superGeometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = superGeometry.getStatistics().getMaxPhysR(2)[0] -
              T {0.5} * converter.getPhysDeltaX();
  extend[0] = T {1} * converter.getPhysDeltaX();
  IndicatorCuboid3D<T> outflow(extend, origin);
  // MN 4 is set only if the neighbor cell has MN 1
  superGeometry.rename(2, 4, 1, outflow);

  superGeometry.clean();

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(SuperLattice<T, DESCRIPTOR>&        sLattice,
                    UnitConverter<T, DESCRIPTOR> const& converter,
                    SuperGeometry<T, 3>&                superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 3, 4});
  sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  // Setting of the boundary conditions
  // local boundary conditions
  boundary::set<boundary::LocalPressure>(sLattice, superGeometry, 3);
  boundary::set<boundary::LocalPressure>(sLattice, superGeometry, 4);

  // interpolated boundary conditions
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  // Initial conditions
  AnalyticalConst3D<T, T> rhoF(T {1});
  const Vector<T, 3>      velocityV {T {0}, T {0}, T {0}};
  AnalyticalConst3D<T, T> uF(velocityV);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(bulkIndicator, rhoF, uF);
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues(SuperLattice<T, DESCRIPTOR>&        sLattice,
                       UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                       SuperGeometry<T, 3>& superGeometry, int iTmaxStart)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  // No of time steps for smooth start-up
  int iTupdate = std::max(iTmaxStart / 200, 1);

  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T, int> StartScale(iTmaxStart, T(1));

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T   frac[1]  = {};
    StartScale(frac, iTvec);
    const T                 currentPressure = frac[0] * pressureDrop;
    AnalyticalConst3D<T, T> rhoFromPressure(
        converter.getLatticeDensityFromPhysPressure(currentPressure));

    sLattice.defineRho(superGeometry, 3, rhoFromPressure);

    clout << "step=" << iT << "; currentPressureDifference=" << currentPressure
          << std::endl;

    sLattice.setProcessingContext<Array<momenta::FixedDensity::RHO>>(
        ProcessingContext::Simulation);
  }
}

void getResults(SuperLattice<T, DESCRIPTOR>&        sLattice,
                UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T, 3>& superGeometry, util::Timer<T>& timer,
                util::ValueTracer<T>& converge, int iTmaxStart)
{

  OstreamManager clout(std::cout, "getResults");

  SuperVTMwriter3D<T>                       vtmWriter("resolvedRock3d");
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);

  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  const int vtkIter  = converter.getLatticeTime(maxPhysT * T {0.05});
  const int statIter = converter.getLatticeTime(maxPhysT * T {0.01});

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR>   cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR>     rank(sLattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);

    vtmWriter.createMasterFile();
  }

  // Writes output on the console
  if (iT % statIter == 0 || converge.hasConverged()) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));

    SuperAverage3D<T> avgVel(velocity, superGeometry, 1);
    int               input[1] {};
    avgVel(currentAverageVelocity.data(), input);
    clout << "Average x-velocity: " << currentAverageVelocity[0] << " m/s"
          << std::endl;

    currentPermeability = currentAverageVelocity[0] * kinematicViscosity *
                          fluidDensity * extent[0] / pressureDrop;

    clout << "Permeability: " << currentPermeability << " m^2" << std::endl;

    if (iT > iTmaxStart) {
      converge.takeValue(currentAverageVelocity[0], true);
    }
  }

  // Writes the vtk files
  if (iT % vtkIter == 0 || converge.hasConverged()) {
    vtmWriter.write(iT);

    {
      SuperEuklidNorm3D<T>  normVel(velocity);
      BlockReduction3D2D<T> planeReduction(normVel, Vector<T, 3>({0, 0, 1}));
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }
  }
}

int main(int argc, char* argv[])
{

  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  std::string vtiFile;
  std::string arrayName;
  T           timeScalingFactor {1};

  std::vector<std::string> cmdInput;
  if (argc > 1) {
    cmdInput.assign(argv + 1, argv + argc);
    vtiFile = cmdInput[0];
  }
  if (argc > 2) {
    arrayName = cmdInput[1];
  }
  if (argc > 3) {
    scalingFactor = std::stod(cmdInput[2]);
  }
  if (argc > 4) {
    timeScalingFactor = std::stod(cmdInput[3]);
  }
  if (argc > 5) {
    N = std::stoi(cmdInput[4]);
  }
  if (argc > 6) {
    pressureDrop = std::stod(cmdInput[5]);
  }
  else {
    clout
        << "Define <filename> <arrayname> <scaling-factor> "
           "<time-scaling-factor> <resolution> <pressure_drop> (in that order)"
        << std::endl;
    return 1;
  }

  clout << "Using the following input:" << std::endl;
  clout << "VTI file: " << vtiFile << std::endl;
  clout << "Array name: " << arrayName << std::endl;
  clout << "Pressure drop: " << pressureDrop << " Pa" << std::endl;

  std::shared_ptr<IndicatorBlockData3D<T>> rock =
      generateIndicatorFromVTI(vtiFile, arrayName);
  clout << "Rock x-length: " << extent[0] << " m" << std::endl;
  clout << "Rock y-length: " << extent[1] << " m" << std::endl;
  clout << "Rock z-length: " << extent[2] << " m" << std::endl;

  constexpr T permeabilityForVelocityEstimation = T {1e-7};
  // Estimating velocity using Darcy's law
  const T charPhysVelocityEstimation =
      permeabilityForVelocityEstimation * pressureDrop /
      (kinematicViscosity * fluidDensity * extent[0]);
  // Assumption. Likely needs an update after some experience is gained
  maxPhysT = timeScalingFactor * extent[0] / charPhysVelocityEstimation;

  clout << "Timeframe to be simulated: " << maxPhysT << " s" << std::endl;

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
      int {N}, // resolution: number of voxels per charPhysL
      (T)tau, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T)extent[1], // charPhysLength: reference length of simulation geometry
      (T)charPhysVelocityEstimation, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)kinematicViscosity, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)fluidDensity // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("rock");

  // === 2nd Step: Prepare Geometry ===
  // Instantiation of a cuboidDecomposition with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif

  //layer of 1 cell around the rock geometry for correct definition of bondary matreial numbers
  IndicatorLayer3D<T> layer(*rock, converter.getPhysDeltaX());

  CuboidDecomposition3D<T> cuboidDecomposition(
      layer, converter.getPhysDeltaX(), noOfCuboids);
  cuboidDecomposition.setPeriodicity({false, true, true});

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // Instantiation of a superGeometry
  SuperGeometry<T, 3> superGeometry(cuboidDecomposition, loadBalancer);

  prepareGeometry(converter, layer, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);

  //prepareLattice and set boundaryCondition
  prepareLattice(sLattice, converter, superGeometry);

  // === 4th Step: Main Loop with Timer ===
  int                  iTmaxStart = converter.getLatticeTime(maxPhysT * 0.4);
  int                  interval = converter.getLatticeTime(maxPhysT / T {100});
  constexpr T          epsilon  = 1e-6;
  util::ValueTracer<T> converge(interval, epsilon);
  clout << "starting simulation..." << std::endl;
  clout << "MaxIT: " << converter.getLatticeTime(maxPhysT) << std::endl;
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT),
                       superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(sLattice, converter, iT, superGeometry, iTmaxStart);

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer, converge,
               iTmaxStart);

    if (converge.hasConverged()) {
      getResults(sLattice, converter, iT, superGeometry, timer, converge,
                 iTmaxStart);
      break;
    }
  }

  timer.stop();
  timer.printSummary();
}

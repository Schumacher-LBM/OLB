/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2006-2025 Mathias J. Krause,
 *  David Heidenthal, Adrian Kummerländer, Michael Grinschewski,
 *  Fedor Bukreev
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

/* airfoil2d.cpp:
 * This example examines a steady flow past a NACA airfoil placed in a channel.
 * At the inlet, a random Turbulent profile is imposed on the velocity,
 * whereas the outlet implements a Dirichlet pressure condition set by p = 0.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T          = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<POROSITY,VELOCITY,TENSOR>;

// Parameters for the simulation setup
const int N              = 50;      // resolution of the model
const T   Re             = 1400000; // Reynolds number
const T   intensity      = 0.02;
const T   maxPhysT       = 16.;      // max. simulation time in s, SI unit
const T   L              = 0.1 / N; // latticeL
const T   lengthX        = 6.;
const T   lengthY        = 2. + L;
const T   centerAirfoilX = lengthX / 4.;
const T   centerAirfoilY = lengthY / 2.;
const T   chordLength    = 1.; // Length of airfoil from leading edge to trailing edge
const T   angleOfAttack  = 5.*std::numbers::pi/90.;

// The four digits in the index are the digits of the following numbers so NACA 1410 would be
// camber = 0.01, camberPos = 0.4, thicknessPercentage = 0.10
const T camber              = 0.01;
const T camberPos           = 0.4;
const T thicknessPercentage = 0.1; // max thickness to chord length ratio

/// Wallfunction parameters
//  Used method for density reconstruction
//  0: use local density
//  1: extrapolation (Guo)
//  2: constant (rho = 1.)
const int rhoMethod = 0;
//  Used wall profile
//  0: power law profile
//  1: Spalding profile
const int wallFunctionProfile = 1;
// check if descriptor with body force is used
const bool bodyForce = false;
// interpolate sampling velocity along given normal between lattice voxels
const bool interpolateSampleVelocity = true;
// use van Driest damping function for turbulent viscosity in boundary cell
const bool useVanDriest = false;
//  distance from cell to real wall in lattice units if no geometry indicator is given as input
const T latticeWallDistance = 0.5;
//  distance from cell to velocity sampling point in lattice units
const T samplingCellDistance = 3.5;
const bool movingWall = false;
const bool averageVelocity = false;

template <typename T, typename S>
class TurbulentVelocity2D : public AnalyticalF2D<T, S> {

protected:
  T u0, intensity;

public:
  TurbulentVelocity2D(T _u0, T _intensity)
      : AnalyticalF2D<T, S>(2)
  { // ======================================
    u0              = _u0;
    intensity       = _intensity;
    this->getName() = "turbulentVelocity2d";
  };

  bool operator()(T output[2], const S input[2]) override
  {
    T u_calc = u0;

    std::random_device          rd;
    std::mt19937                generator(rd());
    std::normal_distribution<T> dist(-1, 1.);

    output[0] = u_calc + intensity * u0 * dist(generator);
    output[1] = intensity * u0 * dist(generator);

    return true;
  };
};

// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter,
                     SuperGeometry<T, 2>&                superGeometry,
                     IndicatorF2D<T>&                    airfoil)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T, 2> extend(lengthX, lengthY);
  Vector<T, 2> origin;

  superGeometry.rename(0, 2);

  superGeometry.rename(2, 1, {1, 1});

  // Set material number for inflow
  extend[0] = 2. * L;
  origin[0] = -L;
  IndicatorCuboid2D<T> inflow(extend, origin);
  superGeometry.rename(2, 3, 1, inflow);
  // Set material number for outflow
  origin[0] = lengthX - L;
  IndicatorCuboid2D<T> outflow(extend, origin);
  superGeometry.rename(2, 4, 1, outflow);
  // Set material number for airfoil
  superGeometry.rename(1, 5, airfoil);

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>&        sLattice,
                    UnitConverter<T, DESCRIPTOR> const& converter,
                    SuperGeometry<T, 2>&                superGeometry,
                    IndicatorF2D<T>&                    airfoil,
                    WallModelParameters<T>&             wallModelParameters)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  // auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  // sLattice.defineDynamics<SmagorinskyBGKdynamics>(bulkIndicator);
  setTurbulentWallModelDynamics(sLattice, superGeometry, 1,
                                wallModelParameters);

  // Material=2 -->Full Slip
  boundary::set<boundary::FullSlip>(sLattice, superGeometry, 2); // <-

  // Setting of the boundary conditions

  //if boundary conditions are chosen to be local

  //if boundary conditions are chosen to be interpolated
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  setBouzidiBoundary(sLattice, superGeometry, 5, airfoil);
  setTurbulentWallModel(sLattice, superGeometry, 5, wallModelParameters);
  // Initial conditions
  AnalyticalConst2D<T, T> rhoF(1);
  AnalyticalConst2D<T, T> rho0(0);
  std::vector<T>          velocity(2, T(0));
  AnalyticalConst2D<T, T> uF(velocity);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry.getMaterialIndicator({1,2,3,4,5}), rhoF, uF);
  sLattice.iniEquilibrium(superGeometry.getMaterialIndicator({1,2,3,4,5}), rhoF, uF);
  sLattice.defineField<VELOCITY>(superGeometry.getMaterialIndicator({1,2,3,4,5}), uF);
  sLattice.defineField<POROSITY>(superGeometry.getMaterialIndicator({1,2,3,4}), rhoF);
  sLattice.defineField<POROSITY>(superGeometry, 5, rho0);

  sLattice.setParameter<descriptors::OMEGA>(omega);
  sLattice.setParameter<collision::LES::SMAGORINSKY>(0.3);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues(SuperLattice<T, DESCRIPTOR>&        sLattice,
                       UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                       SuperGeometry<T, 2>& superGeometry)
{

  OstreamManager clout(std::cout, "setBoundaryValues");

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime(maxPhysT * 0.1);
  int iTupdate   = 5;

  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T, T> StartScale(iTmaxStart, T(1));

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T(iT)};
    T frac[1]  = {};
    StartScale(frac, iTvec);
    T maxVelocity = converter.getCharLatticeVelocity() * frac[0];
    TurbulentVelocity2D<T, T> poiseuilleU(maxVelocity, intensity);

    sLattice.defineU(superGeometry, 3, poiseuilleU);

    sLattice.setProcessingContext<
        Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults(SuperLattice<T, DESCRIPTOR>&        sLattice,
                UnitConverter<T, DESCRIPTOR> const& converter, std::size_t iT,
                SuperGeometry<T, 2>& superGeometry, util::Timer<T>& timer)
{
  OstreamManager clout(std::cout, "getResults");

  // Gnuplot constructor (must be static!)
  // for real-time plotting: gplot("name", true) // experimental!
  static Gnuplot<T> gplot("drag");

  SuperVTMwriter2D<T>                       vtmWriter("airfoil2d");
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, converter);
  SuperGeometryF<T, 2>                      materials(superGeometry);
  SuperLatticeField2D<T, DESCRIPTOR, descriptors::Y1> y1(sLattice);
  SuperLatticePhysField2D<T, DESCRIPTOR, descriptors::WMVELOCITY> wmvelocity(sLattice, converter.getConversionFactorVelocity());
  wmvelocity.getName() = "wmvelocity";

  vtmWriter.addFunctor(materials);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);
  vtmWriter.addFunctor(y1);
  vtmWriter.addFunctor(wmvelocity);

  const int vtkIter  = converter.getLatticeTime(.3);
  const int statIter = converter.getLatticeTime(.1);

  T point[2] = {};
  point[0]   = centerAirfoilX + chordLength;
  point[1]   = centerAirfoilY;
  AnalyticalFfromSuperF2D<T> intpolateP(pressure, true);
  T                          p;
  intpolateP(&p, point);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  if (iT % statIter == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF2D<T>            intpolatePressure(pressure, true);
    SuperLatticePhysDrag2D<T, DESCRIPTOR> drag(sLattice, superGeometry, 5,
                                               converter);

    T point1[2] = {};
    T point2[2] = {};

    point1[0] = centerAirfoilX - chordLength / 2;
    point1[1] = centerAirfoilY;

    point2[0] = centerAirfoilX + chordLength / 2;
    point2[1] = centerAirfoilY;

    T p1, p2;
    intpolatePressure(&p1, point1);
    intpolatePressure(&p2, point2);

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1 - p2;
    clout << "; pressureDrop=" << pressureDrop;

    int input[3] = {};
    T   _drag[drag.getTargetDim()];
    drag(_drag, input);
    clout << "; drag=" << _drag[0] << "; lift=" << _drag[1] << std::endl;

    // set data for gnuplot: input={xValue, yValue(s), names (optional), position of key (optional)}
    gplot.setData(converter.getPhysTime(iT), {_drag[0], 5.58},
                  {"drag(openLB)", "drag(schaeferTurek)"}, "bottom right",
                  {'l', 'l'});

    // every (iT%vtkIter) write an png of the plot
    if (iT % (vtkIter) == 0) {
      // writes pngs: input={name of the files (optional), x range for the plot (optional)}
      gplot.writePNG(iT, maxPhysT);
    }
  }

  // Writes the vtk files
  if (iT % vtkIter == 0) {
    vtmWriter.write(iT);
  }

  // write pdf at last time step
  if (iT == converter.getLatticeTime(maxPhysT) - 1) {
    // writes pdf
    gplot.writePNG();
  }
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> const converter(
      int {N},        // resolution: number of voxels per charPhysL
      (T)0.03,        // Max cell speed?
      (T)chordLength, // charPhysLength: reference length of simulation geometry
      (T)Re * 1.e-5 /
          chordLength, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)1.e-5, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)1.0    // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("airfoil2d");

  // === 2rd Step: Prepare Geometry ===
  Vector<T, 2>         extend(lengthX, lengthY);
  Vector<T, 2>         origin;
  IndicatorCuboid2D<T> cuboid(extend, origin);

  // Instantiation of a cuboidGeometry with weights
  CuboidDecomposition2D<T> cuboidGeometry(cuboid, L, singleton::mpi().getSize());

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  // Instantiation of a superGeometry
  SuperGeometry<T, 2> superGeometry(cuboidGeometry, loadBalancer);

  Vector<T, 2>          center(centerAirfoilX, centerAirfoilY);
  IndicatorAirfoil2D<T> airfoilI = IndicatorAirfoil2D<T>(
      center, chordLength, camber, camberPos, thicknessPercentage);
  IndicatorRotate<T,2> airfoil(center, angleOfAttack, airfoilI);
  prepareGeometry(converter, superGeometry, airfoil);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);
  WallModelParameters<T>      wallModelParameters;
  wallModelParameters.bodyForce                 = bodyForce;
  wallModelParameters.rhoMethod                 = rhoMethod;
  wallModelParameters.samplingCellDistance      = samplingCellDistance;
  wallModelParameters.interpolateSampleVelocity = interpolateSampleVelocity;
  wallModelParameters.useVanDriest              = useVanDriest;
  wallModelParameters.wallFunctionProfile       = wallFunctionProfile;
  wallModelParameters.latticeWallDistance       = latticeWallDistance;
  wallModelParameters.movingWall                = movingWall;
  wallModelParameters.averageVelocity           = averageVelocity;

  //prepareLattice and set boundaryConditions
  prepareLattice(sLattice, converter, superGeometry, airfoil,
                 wallModelParameters);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT),
                       superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(sLattice, converter, iT, superGeometry);

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer);
  }

  timer.stop();
  timer.printSummary();
}

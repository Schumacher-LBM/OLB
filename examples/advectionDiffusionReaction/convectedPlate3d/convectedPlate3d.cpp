/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2023 Shota Ito
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

// This app is used for the validation of the ADE-LBM and of the Dirchlet and Zero-Gradient BC regarding the transported concentration.
// The simulation setup and analytic solution is taken from Krüger et al. 2017 page 324 (doi: https://doi.org/10.1007/978-3-319-44649-3).

#include <olb.h>
#include <cmath>
using namespace olb;
using namespace olb::descriptors;

using S = FLOATING_POINT_TYPE;
using U = util::ADf<S,1>;
using RADDESCRIPTOR = D3Q7<VELOCITY,OMEGA>; // the scalar field contains material values if a cell should be considered for the zero gradient BC or not

// parameters for the simulation setup
const int N = 40;                                 // resolution
const S height = 160;                             // height of the channel (y-dir)
const S dx = height / N;                          // lattice length
const S length = N * dx * 10;                     // length of the channel (x-dir)
const S width = dx;                               // width of the channel (z-dir)
const S u0 = 0.05;                                // Bulk velocity
const S Pe = 108.005;                             // Peclet number
const S maxPhysT = 100000.0;                      // max. physical simulation time
const S physInterval = maxPhysT * 0.001;          // interval for convergence check
const S residuum = 1e-10;                         // residuum for convergence check

// Analytic concentration field for the convected plate
template<typename T>
class ConvPlate3D : public AnalyticalF3D<T,T> {

protected:
  T diff;
  T bulkVel;
  T dx;
public:
  ConvPlate3D(AdeUnitConverter<T,RADDESCRIPTOR> converter) :
    AnalyticalF3D<T,T>(1),
    diff(converter.getPhysDiffusivity()),
    bulkVel(converter.getCharPhysVelocity()) {};

  bool operator()(T output[], const T input[]) override {
    T x = input[0];
    T y = input[1];

    if (x >= 0.0 && y >= 0.0) {
      output[0] = erfc(y / sqrt(4 * diff * x / bulkVel));
    }
    else {
      output [0] = 0.0;
    }
    return true;
  };
};

// write geometry file to check geometry
template <typename T>
void writeGeometry3D(SuperVTMwriter3D<T> vtmWriter,
                     SuperGeometry<T,3>& superGeometry,
                     SuperLattice<T,RADDESCRIPTOR>& adLattice)
{

  vtmWriter.createMasterFile();
}

// prepare numerical domain
template<typename T>
void prepareGeometry(AdeUnitConverter<T,RADDESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 1);

  // === MATERIAL NUMBERS === //
  // solid: MN = 0            //
  // bulk: MN = 1             //
  // inlet: MN = 2            //
  // outlet: MN = 3           //
  // topWall: MN = 4          //
  // bottomWall: MN = 5       //
  // ======================== //

  T eps = dx * T(0.5);
  {
  Vector<T,3>origin(-eps, -eps, -eps - dx);
  Vector<T,3>extend(eps * T(2.), height + T(2.) * eps, width + T(2.) * eps + T(2.) * dx);
  IndicatorCuboid3D<T> inflow(extend, origin);
  superGeometry.rename(1, 2, inflow);

  origin[0] += length;
  IndicatorCuboid3D<T> outflow(extend, origin);
  superGeometry.rename(1, 3, outflow);
  }

  {
  Vector<T,3>origin(eps, height - eps, -eps - dx);
  Vector<T,3>extend(length - eps * T(2.), eps * T(2.), width + T(2.) * eps + T(2.) * dx);
  IndicatorCuboid3D<T> topWall(extend, origin);
  superGeometry.rename(1, 4, topWall);

  origin[1] -= height;
  IndicatorCuboid3D<T> bottomWall(extend, origin);
  superGeometry.rename(1, 5, bottomWall);
  }

  // Clean gemetry
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

template<typename T>
void prepareLattice(AdeUnitConverter<T,RADDESCRIPTOR> const& converter,
                    SuperLattice<T,RADDESCRIPTOR>& adLattice,
                    SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // get the relaxation time from unitConverter
  const T omega = 1.0 / converter.getLatticeRelaxationTime();

  // define dynamics and boundary conditions
  const auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3, 4, 5});

  adLattice.template defineDynamics<NoDynamics>(superGeometry, 0);
  adLattice.template defineDynamics<AdvectionDiffusionBGKdynamics>(bulkIndicator);

  boundary::set<boundary::AdvectionDiffusionDirichlet>(adLattice, superGeometry.getMaterialIndicator({2}));
  setZeroGradientBoundary<T, RADDESCRIPTOR>(adLattice, superGeometry.getMaterialIndicator({3}));
  setZeroGradientBoundary<T, RADDESCRIPTOR>(adLattice, superGeometry.getMaterialIndicator({4}));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(adLattice, superGeometry.getMaterialIndicator({5}));

  // Initial conditions
  AnalyticalConst3D<T,T> rho0(1.e-8);
  AnalyticalConst3D<T,T> rho_plate(1.0);
  Vector<T,3> velocity(converter.getCharLatticeVelocity(), 0.0, 0.0);
  AnalyticalConst3D<T,T> u(velocity);

  // set the convective velocity field
  adLattice.template defineField<descriptors::VELOCITY>(bulkIndicator, u);

  // set initail population to equilibrium
  adLattice.defineRho(superGeometry.getMaterialIndicator({1, 2, 3, 4}), rho0);
  adLattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 2, 3, 4}), rho0, u);
  adLattice.defineRho(superGeometry, 5, rho_plate);
  adLattice.iniEquilibrium(superGeometry, 5, rho_plate, u);

  adLattice.template setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  adLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// compute error to analytical solution
template<typename T>
void error(SuperGeometry<T,3>& superGeometry,
           AdeUnitConverter<T,RADDESCRIPTOR> const& converter,
           SuperLattice<T,RADDESCRIPTOR>& adLattice,
           T& relError)
{
  OstreamManager clout(std::cout, "Error2Analytic");

  T result[3] = {T(), T(), T()};
  int tmp[] = {int()};

  // Compute L2-norm error to the analytic concentration field
  SuperLatticeDensity3D<T,RADDESCRIPTOR> concentration(adLattice);
  ConvPlate3D<T> analyticSol(converter);
  auto indicatorF = superGeometry.getMaterialIndicator({1});

  SuperRelativeErrorL2Norm3D<T> relCErrorL2Norm(concentration, analyticSol, indicatorF);
  SuperAbsoluteErrorL2Norm3D<T> absCErrorL2Norm(concentration, analyticSol, indicatorF);

  relCErrorL2Norm(result, tmp);
  relError = result[0];
}

// compute results and generate output
template<typename T>
void getResults(SuperLattice<T,RADDESCRIPTOR>& adLattice,
                AdeUnitConverter<T,RADDESCRIPTOR> const& converter,
                std::size_t iT,
                SuperGeometry<T,3>& superGeometry,
                util::Timer<T>& timer,
                util::ValueTracer<T> converge,
                T& relError)
{
  OstreamManager clout(std::cout, "Results");

  // generate vtm files for visualisation
  SuperVTMwriter3D<T> vtmWriter("convectedPlate");
  SuperLatticeDensity3D<T,RADDESCRIPTOR> concentration(adLattice);
  SuperLatticePhysField3D<T,RADDESCRIPTOR,VELOCITY> velocity(adLattice, converter.getConversionFactorVelocity());
  ConvPlate3D<T> analyticSol(converter);
  SuperLatticeFfromAnalyticalF3D<T,RADDESCRIPTOR> analyticalC(analyticSol, adLattice);

  concentration.getName() = "C";
  velocity.getName() = "u0";
  analyticalC.getName() = "analyticC";
  vtmWriter.addFunctor(concentration);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(analyticalC);

  const std::size_t vtkIter = converter.getLatticeTime(maxPhysT / 100);
  const std::size_t statIter = converter.getLatticeTime(maxPhysT / 100);

  if (iT == 0) {
    // write the geometry
    writeGeometry3D(vtmWriter, superGeometry, adLattice);

    // write inital field
    vtmWriter.write(iT);
  }

  if (iT%statIter == 0) {
    adLattice.setProcessingContext(ProcessingContext::Evaluation);

    // print L2 error
    error<T>(superGeometry, converter, adLattice, relError);

    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    adLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    clout << "concentration-L2-relative-error = " << relError << std::endl;

    converge.takeValue( relError, false );
  }

  // write the vtk files
  if (iT % vtkIter == 0 && iT > 0) {
    vtmWriter.write(iT);
  }
}

// main function executed to run simulation
template<typename T>
int simulateFlow(T D)
{
  // === 1. Initialization ===
  OstreamManager clout(std::cout, "main");

  AdeUnitConverter<T, RADDESCRIPTOR> const converter(
    (T)   height / N,                        // physDeltaX
    (T)   util::pow(height / N, 2),          // physDeltaY
    (T)   height,                            // charPhysLength: reference length of simulation geometry
    (T)   Pe * D / height,                   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   D,                                 // physDiffusiviy
    (T)   1.0                                // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  // === 2. Prepare geometry ===
  const Vector<T,3> origin;
  const Vector<T,3> extend(length, height, width);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  // Instantiation of a cuboidDecomposition with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidDecomposition3D<T> cuboidDecomposition(cuboid, dx, noOfCuboids);

  // periodic in the z dir
  cuboidDecomposition.setPeriodicity({false, false, true});

  // load balancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // superGeometry
  const int overlap = 3;
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, overlap);

  prepareGeometry<T>(converter, superGeometry);

  // === 3. Prepare lattice ====
  SuperLattice<T,RADDESCRIPTOR> adLattice(superGeometry);

  prepareLattice<T>(converter, adLattice, superGeometry);

  // === 4. Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());

  // Convergence check
  util::ValueTracer<T> converge(converter.getLatticeTime(physInterval), residuum);
  timer.start();
  T relError{};

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    // === 5. Set boundary values ===

    // === 6. Collide and stream ===
    adLattice.collideAndStream();

    // === 7. Computation and output of the results ===
    getResults<T>(adLattice, converter, iT, superGeometry, timer,  converge, relError);
    if (converge.hasConverged()) {
      clout << "Simulation stops." << std::endl;
      break;
    }
  }

  timer.stop();
  timer.printSummary();

  return 0;
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  if constexpr(true) {
    const S diffusivity = 0.07407;
    simulateFlow<S>(diffusivity);
  }

  if constexpr(false) {
    U diffusivity = 0.07407;
    diffusivity.setDiffVariable(0);
    simulateFlow<U>(diffusivity);
  }
}

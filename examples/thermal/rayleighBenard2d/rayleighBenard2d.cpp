/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
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

/* rayleighBenard2d.cpp:
 * Rayleigh-Benard convection rolls in 2D, simulated with
 * the thermal LB model by Z. Guo e.a., between a hot plate at
 * the bottom and a cold plate at the top.
 */


#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;

using NSDESCRIPTOR = D2Q9<FORCE>;
using TDESCRIPTOR = D2Q5<VELOCITY>;

// Parameters for the simulation setup
const T lx  = 2.;          // length of the channel
const T ly  = 1.;          // height of the channel
const int N = 10;          // resolution of the model
const T Ra = 1e4;          // Rayleigh number
const T Pr = 0.71;         // Prandtl number
const T maxPhysT = 1000.; // max. simulation time in s, SI unit
const T epsilon = 1.e-5;   // precision of the convergence (residuum)

const T Thot = 274.15;     // temperature of the lower wall in Kelvin
const T Tcold = 273.15;    // temperature of the fluid in Kelvin
const T Tperturb = 1./5. * Tcold + 4./5. * Thot; // temperature of the perturbation

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry<T,2>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2);
  superGeometry.rename(2,1,{0,1});

  std::vector<T> extend( 2, T(0) );
  extend[0] = lx;
  extend[1] = converter.getPhysLength(1);
  std::vector<T> origin( 2, T(0) );
  IndicatorCuboid2D<T> bottom(extend, origin);

  origin[1] = ly-converter.getPhysLength(1);
  IndicatorCuboid2D<T> top(extend, origin);

  origin[0] = lx/2.;
  origin[1] = converter.getPhysLength(1);
  extend[0] = converter.getPhysLength(1);
  extend[1] = converter.getPhysLength(1);
  IndicatorCuboid2D<T> perturbation(extend, origin);

  /// Set material numbers for bottom, top and pertubation
  superGeometry.rename(2,2,1,bottom);
  superGeometry.rename(2,3,1,top);
  superGeometry.rename(1,4,perturbation);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                     SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice<T, TDESCRIPTOR>& ADlattice,
                     SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");

  T Tomega  = converter.getLatticeThermalRelaxationFrequency();
  T NSomega = converter.getLatticeRelaxationFrequency();

  NSlattice.defineDynamics<ForcedBGKdynamics>(superGeometry, 1);
  boundary::set<boundary::BounceBack>(NSlattice, superGeometry, 2);
  boundary::set<boundary::BounceBack>(NSlattice, superGeometry, 3);
  NSlattice.defineDynamics<ForcedBGKdynamics>(superGeometry, 4);

  ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(superGeometry, 1);
  ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(superGeometry, 2);
  ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(superGeometry, 3);
  ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(superGeometry, 4);

  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, superGeometry, 2);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, superGeometry, 3);

  /// define initial conditions
  AnalyticalConst2D<T,T> rho(1.);
  AnalyticalConst2D<T,T> u0(0.0, 0.0);
  AnalyticalConst2D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst2D<T,T> T_hot(converter.getLatticeTemperature(Thot));
  AnalyticalConst2D<T,T> T_perturb(converter.getLatticeTemperature(Tperturb));

  /// for each material set Rho, U and the Equilibrium
  auto indicator = superGeometry.getMaterialIndicator({1,2,3,4});
  NSlattice.defineRhoU(indicator, rho, u0);
  NSlattice.iniEquilibrium(indicator, rho, u0);

  ADlattice.defineRho(superGeometry, 1, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 1, T_cold, u0);
  ADlattice.defineRho(superGeometry, 2, T_hot);
  ADlattice.iniEquilibrium(superGeometry, 2, T_hot, u0);
  ADlattice.defineRho(superGeometry, 3, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 3, T_cold, u0);
  ADlattice.defineRho(superGeometry, 4, T_perturb);
  ADlattice.iniEquilibrium(superGeometry, 4, T_perturb, u0);

  NSlattice.setParameter<descriptors::OMEGA>(NSomega);
  ADlattice.setParameter<descriptors::OMEGA>(Tomega);

  /// Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                SuperLattice<T, TDESCRIPTOR>& ADlattice, int iT,
                SuperGeometry<T,2>& superGeometry,
                util::Timer<T>& timer,
                bool converged)
{
  SuperVTMwriter2D<T> vtkWriter("rayleighBenard2d");

  const int statIter = converter.getLatticeTime(10.0);
  const int saveIter = converter.getLatticeTime(10.0);

  if (iT == 0) {
    /// Writes the converter log file
    // writeLogFile(converter,"rayleighBenard2d");

    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  if (iT%statIter == 0 || converged) {
    /// Timer console output
    timer.update(iT);
    timer.printStep();

    /// Lattice statistics console output
    NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));
  }

  /// Writes the VTK files and prints statistics
  if (iT%saveIter == 0 || converged) {
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity(NSlattice, converter);
    SuperLatticePhysPressure2D<T, NSDESCRIPTOR> presure(NSlattice, converter);
    SuperLatticePhysTemperature2D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);
    vtkWriter.addFunctor( presure );
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( temperature );

    vtkWriter.write(iT);

    BlockReduction2D2D<T> planeReduction(temperature, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tcold-0.1, Thot+0.1, iT, "temperature");
  }

}

int main(int argc, char *argv[])
{
  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> converter(
    (T) 0.1/N, // physDeltaX
    (T) 0.1 / (1e-5 / 0.1 * util::sqrt( Ra / Pr)) * 0.1 / N, // physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX
    (T) 0.1,  // charPhysLength
    (T) 1e-5 / 0.1 * util::sqrt( Ra / Pr ), // charPhysVelocity
    (T) 1e-5,  // physViscosity
    (T) 1.0, // physDensity
    (T) 0.03, // physThermalConductivity
    (T) Pr * 0.03 / 1e-5 / 1.0,    // physSpecificHeatCapacity
    (T) Ra * 1e-5 * 1e-5 / Pr / 9.81 / (Thot - Tcold) / util::pow(0.1, 3), // physThermalExpansionCoefficient
    (T) Tcold, // charPhysLowTemperature
    (T) Thot // charPhysHighTemperature
  );
  converter.print();

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  extend[0] = lx;
  extend[1] = ly;
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of a cuboidDecomposition with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), noOfCuboids);

  cuboidDecomposition.setPeriodicity({true, false});

  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer);

  prepareGeometry(superGeometry, converter);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice<T, NSDESCRIPTOR> NSlattice(superGeometry);
  SuperLattice<T, TDESCRIPTOR> ADlattice(superGeometry);

  prepareLattice(converter, NSlattice, ADlattice, superGeometry);

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
  // This coupling must be necessarily be put on the Navier-Stokes lattice!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                               converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();
  SuperLatticeCoupling coupling(
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSlattice,
    names::Temperature{}, ADlattice);
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
    converter.getLatticeTemperature(Tcold));
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});

  /// === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  util::ValueTracer<T> converge(converter.getLatticeTime(50.),epsilon);
  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    if (converge.hasConverged()) {
      clout << "Simulation converged." << std::endl;
      getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());

      clout << "Time " << iT << "." << std::endl;

      break;
    }

    /// === 6th Step: Collide and Stream Execution ===
    NSlattice.collideAndStream();
    ADlattice.collideAndStream();

    coupling.execute();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
    converge.takeValue(ADlattice.getStatistics().getAverageEnergy(),true);
  }

  timer.stop();
  timer.printSummary();
}

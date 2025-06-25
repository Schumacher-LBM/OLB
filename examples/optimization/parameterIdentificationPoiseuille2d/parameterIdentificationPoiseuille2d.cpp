/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012, 2022 Jonas Latt, Mathias J. Krause
 *  Vojtech Cvrcek, Peter Weisbrod, Julius Jeßberger
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

/** \file A simple two-dimensional fluid flow optimization problem is solved.
 * The setup is a planar channel flow, similar to the example poiseuille2d.
 * For a given pressure drop, the steady flow is simulated and the mass flow
 * rate is computed. In the optimization problem, the inlet pressure is
 * determined s.t. a pre-defined mass flow rate is achieved.
 */

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

using namespace olb;
using namespace olb::opti;

using S = FLOATING_POINT_TYPE;
using U = util::ADf<S,1>;
using DESCRIPTOR = descriptors::D2Q9<>;

template <typename T>
using VectorHelp = Vector<T,1>;

const int N = 50;             // resolution
const S lx  = 2.;             // length of the channel
const S ly  = 1.;             // height of the channel
const S Re = 10.;             // Reynolds number
const S maxPhysT = 30.;       // max. simulation time in s, SI unit
const S physInterval = 0.25;  // interval for the convergence check in s
const S residuum = 1e-9;      // residuum for the convergence check
const S wantedMassFlow = 0.00026159;


template<typename T>
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry )
{
  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1} );

  const T physSpacing = converter.getPhysDeltaX();
  const Vector<T,2> extend    {physSpacing / T(2), ly};
  const Vector<T,2> originIn  {-physSpacing / T(4), 0};
  const Vector<T,2> originOut {lx-physSpacing / T(4), 0};

  IndicatorCuboid2D<T> inflow( extend, originIn );
  superGeometry.rename( 2,3,1,inflow );

  IndicatorCuboid2D<T> outflow( extend, originOut );
  superGeometry.rename( 2,4,1,outflow );

  superGeometry.clean(false);
  superGeometry.innerClean(false);
  superGeometry.checkForErrors(false);
}

template<typename T>
void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry,
                     T inletPressure)
{
  const T omega = converter.getLatticeRelaxationFrequency();

  sLattice.template defineDynamics<BGKdynamics<T,DESCRIPTOR>>(superGeometry, 1);
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  // Initial conditions
  AnalyticalLinear2D<T,T> rho(
    - inletPressure / lx * descriptors::invCs2<T,DESCRIPTOR>(),
    0,
    inletPressure * descriptors::invCs2<T,DESCRIPTOR>() + 1 );

  const T Lx = converter.getLatticeLength( lx ) - 1;
  const T Ly = converter.getLatticeLength( ly ) - 1;
  const T maxVelocity = inletPressure * Ly * Ly
    / (8.0 * converter.getLatticeViscosity() * Lx);
  const T radius = T(0.5) * (ly - converter.getPhysDeltaX());
  const std::vector<T> axisPoint { lx/T(2), ly/T(2) };
  const std::vector<T> axisDirection { 1, 0 };
  Poiseuille2D<T> u( axisPoint, axisDirection, maxVelocity, radius );

  const std::vector<T> zero(2, T());
  AnalyticalConst2D<T, T> u0(zero);

  sLattice.defineRhoU(superGeometry, 0, rho, u0);
  sLattice.iniEquilibrium(superGeometry, 0, rho, u0);

  const auto domain = superGeometry.getMaterialIndicator({1,2,3,4});
  sLattice.defineRhoU( domain, rho, u );
  sLattice.iniEquilibrium( domain, rho, u );

  sLattice.template setParameter<descriptors::OMEGA>(omega);

  sLattice.initialize();
}

template<typename T>
T getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer)
{
  OstreamManager clout( std::cout,"getResults" );
  sLattice.communicate();
  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  // Output on the console
  timer.update( iT );
  timer.printStep();
  sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

  SuperLatticeVelocity2D velocity(sLattice);
  SuperLatticeDensity2D density(sLattice);
  SuperPlaneIntegralFluxMass2D<T> massFlowRate(
    velocity, density, superGeometry, converter.getConversionFactorMass(),
    converter.getConversionFactorTime(), Vector<T,2>({T(0.5)*lx, T(0.5)*ly}),
    Vector<T,2>({0, 1}), BlockDataReductionMode::Analytical
  );
  const int input[3] = {0};
  T mFlow[4] = {0.};
  massFlowRate(mFlow, input);
  clout << "Mass flow rate = " << mFlow[0] << std::endl;
  return mFlow[0];
}

template <typename T>
void writeVTK( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,2>& superGeometry, bool hasConverged)
{
  SuperVTMwriter2D<T> vtmWriter( "poiseuille2d" );
  const bool lastTimeStep
   = ( hasConverged || (iT + 1 == converter.getLatticeTime( maxPhysT )) );
  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );

  SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  if (iT == 0) {
      sLattice.communicate();
      SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
      SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
      vtmWriter.write( cuboid );
      vtmWriter.write( rank );
      vtmWriter.createMasterFile();
  }
  if ( iT%vtmIter==0 || lastTimeStep )
  {
    sLattice.communicate();
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write( iT );
  }
}

/** Perform a single flow simulation.
 * \param inletPressure: The pressure at the inlet. Pressure at outlet is fixed
 * to 0.
 * \param optiMode: this is true when the method is called from an Optimization
 * loop. In that case, less vtk output is produced.
 */
template<typename T, bool optiMode>
T simulatePoiseuille( T inletPressure )
{
  OstreamManager clout( std::cout,"simulatePoiseuille" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   0.8,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,     // charPhysLength: reference length of simulation geometry
    (T)   1,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );

  // === 2nd Step: Prepare Geometry ===
  const Vector<T,2> extend( lx, ly );
  const Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidDecomposition2D<T> cuboidDecomposition(
    cuboid, converter.getPhysDeltaX(), noOfCuboids );

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

  SuperGeometry<T,2> superGeometry( cuboidDecomposition, loadBalancer, 3 );

  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  // Prepare lattice and set boundary conditions
  prepareLattice<T>(converter, sLattice, superGeometry, inletPressure);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation with pressure = " << inletPressure << "..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ),
    superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();

  std::size_t iT {0};
  for ( ; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      if constexpr (! optiMode) {
        writeVTK<T>(sLattice, converter, iT, superGeometry, true);
      }
      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    if constexpr (! optiMode) {
      writeVTK<T>(sLattice, converter, iT, superGeometry, converge.hasConverged());
    }
    converge.takeValue( sLattice.getStatistics().getMaxU(), false );
  }
  const T massFlow = getResults( sLattice, converter, iT, superGeometry, timer );

  timer.stop();

  return massFlow;
}

/// Compute squared error between simulated and wanted mass flow rate
template<typename T, bool optiMode>
T poiseuilleMassFlowError(Vector<T,1> inletPressure)
{
  const T res = simulatePoiseuille<T,optiMode>(inletPressure[0]);
  const T wantedRes {wantedMassFlow};
  return 0.5 * (res - wantedRes) * (res - wantedRes);
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );

  if constexpr (false) {  // direct (standard) simulation
    S pressure {0.000659};
    simulatePoiseuille<S,false>(pressure);
  }

  if constexpr (false) {  // direct simulation with ADf -> compute derivatives
    U pressure {0.000659};
    pressure.setDiffVariable(0);
    simulatePoiseuille<U,false>(pressure);
  }

  if constexpr (true) {  // Optimization (fit mass flow rate) with ADf
    OptiCaseAD<S,1,VectorHelp> optiCase(
      poiseuilleMassFlowError<S,true>,
      poiseuilleMassFlowError<U,true>);
    /*
    // Perform optimization with steepest descent method
    OptimizerSteepestDescent<S,Vector<S,1>> optimizer(
      1, 1.e-7, 20, 1., 10, "Armijo", true, "", "log",
      true, 0.01, true, 0., false, 0.,
      {OptimizerLogType::value, OptimizerLogType::control, OptimizerLogType::derivative});
    */
    // Perform optimization with LBFGS method.
    OptimizerLBFGS<S,Vector<S,1>> optimizer(
      1, 1.e-7, 20, 1., 10, "StrongWolfe", 20, 1.e-4, true, "", "log",
      true, 0.01, true, 0., false, 0., true,
      {OptimizerLogType::value, OptimizerLogType::control, OptimizerLogType::derivative});

    Vector<S,1> startValue {0.0001};
    optimizer.setControl(startValue);
    optimizer.optimize(optiCase);
  }
}

/*Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Luiz Eduardo Czelusniak
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

/* dropletSplashing2d.cpp
 * In this example a droplet falls in a liquid film.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<FORCE,          // Interaction force
                        EXTERNAL_FORCE, // External force such as gravity
                        STATISTIC,      // Store the density field
                        SCALAR>;        // Store the EOS pressure field
using COUPLING = PseudopotentialForcedPostProcessor<interaction::Polinomial>;

// Parameters in physical units
const T radius = 0.5e-4; // radius of liquid phase [m]
const T Lx = 8*radius; // length of the domain [m]
const T Ly = 4*radius; // height of the domain [m]
const T U_droplet = 3.0; // droplet velocity [m.s-1]

// Properties of R134 vapor and liquid
// at saturation condition 30 degree Celsius

const T nu_vapor   = 3.3e-6; // vapor kinematic viscosity [m2.s-1]
const T nu_liquid = 1.55e-6;   // liquid kinematic viscosity [m2.s-1]

const T SurfTension = 7.58e-3; // surface tension for water-air [N.m-1]

const T rho_vapor   = 37.5; // vapor density [kg.m-3]
const T rho_liquid = 1187;   // liquid density [kg.m-3]

// Resolution
const int Nx = 400;
const int Ny = 200;

// Interface parameters
const T thickness = 1.0; // interface thickness in grid nodes; physThickness = delta_x * thickness

// Time and plot parameters
const int maxIter  = 1500; // amount of time steps to complete simulation
const int vtkIter  = 10;   // amount of time steps to save vtk files
const int statIter = 10;   // amount of time steps to display simulation parameters

/**
 * Correction for initial population
 * Take into account the force effect in initial population
 */
struct InitialPopulationCorrectionO {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    const auto force = cell.template getField<descriptors::FORCE>();
    const T rho = cell.computeRho();

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      T c_F{};
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_F += descriptors::c<DESCRIPTOR>(iPop,iD)*force[iD];
      }
      c_F *= descriptors::invCs2<T,DESCRIPTOR>();
      cell[iPop] -= descriptors::t<T,DESCRIPTOR>(iPop) * 0.5 * rho * c_F;
    }

  }
};

template<typename COUPLING_FORCE>
void initialCondition( SuperLattice<T, DESCRIPTOR>& sLattice,
                       COUPLING_FORCE& couplingInteractionForce )
{
  OstreamManager clout( std::cout,"initialCondition" );
  clout << "Initial Condition ..." << std::endl;

  // Compute initial force
  sLattice.executePostProcessors(stage::PreCoupling());
  sLattice.getCommunicator(stage::PreCoupling()).communicate();
  couplingInteractionForce.execute();

  // Correct velocity
  sLattice.addPostProcessor<stage::PostCoupling>(meta::id<InitialPopulationCorrectionO>());
  sLattice.executePostProcessors(stage::PostCoupling());

  clout << "Initial Condition ... OK" << std::endl;

}

void prepareGeometry( SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{0,1} );
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry,
                     MultiPhaseUnitConverterFromRelaxationTime<T, DESCRIPTOR> const& converter )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // define lattice Dynamics
  using BulkDynamics = MultiphaseForcedBGKdynamics<T, DESCRIPTOR>;
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  std::vector<T> dropletVelocity_x = {0.};
  AnalyticalConst2D<T,T> _dropletVelocity_x( dropletVelocity_x );

  SmoothIndicatorFactoredCircle2D<T,T> dropletVelocity_y( {Lx/2., Ly/2.}, radius,
                                                           sqrt(2.)*thickness*converter.getConversionFactorLength(),
                                                           0, {0,0}, 0,
                                                           -U_droplet/converter.getConversionFactorVelocity() );
  AnalyticalIdentity2D<T,T> _dropletVelocity_y( dropletVelocity_y );

  AnalyticalComposed2D<T,T> fluidVelocity(_dropletVelocity_x,_dropletVelocity_y);


  AnalyticalConst2D<T,T> vapor ( rho_vapor/converter.getConversionFactorDensity() );
  SmoothIndicatorFactoredCircle2D<T,T> liquid ( {Lx/2., Ly/2.}, radius,
                                                 sqrt(2.)*thickness*converter.getConversionFactorLength(),
                                                 0, {0,0}, 0,
                                                 ( rho_liquid - rho_vapor )/converter.getConversionFactorDensity() );
  SmoothIndicatorFactoredCuboid2D<T,T> film( {Lx/2., 0.}, 2.*Lx, radius,
                                                 sqrt(2.)*thickness*converter.getConversionFactorLength(),
                                                 0, {0,0}, 0,
                                                 ( rho_liquid - rho_vapor )/converter.getConversionFactorDensity() );

  AnalyticalIdentity2D<T,T> fluidDensity( vapor + liquid + film );

  auto allGeometry = superGeometry.getMaterialIndicator({1,2});
  sLattice.defineRhoU( allGeometry, fluidDensity, fluidVelocity );
  sLattice.iniEquilibrium( allGeometry, fluidDensity, fluidVelocity );

  std::cout << "External force field = 0" << std::endl;
  AnalyticalConst2D<T,T> f( 0 );
  sLattice.defineField<descriptors::EXTERNAL_FORCE>( superGeometry, 1, f );

  // global relaxation frequency (it can be initialized as one)
  AnalyticalConst2D<T,T> one( 1. );
  sLattice.defineField<descriptors::OMEGA>( superGeometry, 1, one );
  sLattice.setParameter<OMEGA>( 1. );
  sLattice.setParameter<multiphase::RHO_VAPOR>( rho_vapor/converter.getConversionFactorDensity() );
  sLattice.setParameter<multiphase::RHO_LIQUID>( rho_liquid/converter.getConversionFactorDensity() );
  sLattice.setParameter<multiphase::OMEGA_VAPOR>( 1./converter.computeRelaxationTimefromPhysViscosity( nu_vapor ) );
  sLattice.setParameter<multiphase::OMEGA_LIQUID>( 1./converter.computeRelaxationTimefromPhysViscosity( nu_liquid ) );


  sLattice.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

  {
    auto& communicator = sLattice.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.exchangeRequests();
  }

  sLattice.initialize();
  sLattice.getCommunicator(stage::PreCoupling()).communicate();
  sLattice.executePostProcessors(stage::PreCoupling());

  clout << "Prepare Lattice ... OK" << std::endl;
}

//std::vector<T>
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice1,
                 int iT, SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer,
                 MultiPhaseUnitConverterFromRelaxationTime<T, DESCRIPTOR> const& converter )
{
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "dropletSplashing2d" );
  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice1 );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice1 );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }
  // Get statistics
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLattice1.getStatistics().print( iT, iT );
  }

  // Save vtk files
  if ( iT%vtkIter==0 ) {
    sLattice1.setProcessingContext(ProcessingContext::Evaluation);

    // rho_lat -> density in lattice units
    SuperLatticeDensity2D<T, DESCRIPTOR> rhoL_lat( sLattice1 );
    SuperIdentity2D<T,T> rho_lat( rhoL_lat );
    rho_lat.getName() = "rho_lat";

    // rho_phs -> density in physical units
    AnalyticalConst2D<T,T> _C_rho( converter.getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> __C_rho( _C_rho, sLattice1 );
    SuperIdentity2D<T,T> rho_phs( rho_lat * __C_rho );
    rho_phs.getName() = "rho_phs";

    // velocity_lat -> velocity in lattice units
    SuperLatticeVelocity2D<T, DESCRIPTOR> velocityL_lat( sLattice1 );
    SuperIdentity2D<T,T> velocity_lat( velocityL_lat );
    velocity_lat.getName() = "velocity_lat";

    // velocity_phs -> velocity in physical units
    AnalyticalConst2D<T,T> _C_U( converter.getConversionFactorVelocity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> __C_U( _C_U, sLattice1 );
    SuperIdentity2D<T,T> velocity_phs( velocity_lat * __C_U );
    velocity_phs.getName() = "velocity_phs";

    // force_lat -> force in lattice units
    SuperLatticeField2D<T, DESCRIPTOR, FORCE> forceL_lat( sLattice1 );
    SuperIdentity2D<T,T> force_lat( forceL_lat );
    force_lat.getName() = "force_lat";

    // force_phs -> force in physical units
    AnalyticalConst2D<T,T> _C_F( converter.getConversionFactorForce() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> __C_F( _C_F, sLattice1 );
    SuperIdentity2D<T,T> force_phs( force_lat * __C_F );
    force_phs.getName() = "force_phs";

    // p_lat -> pressure in lattice units
    SuperLatticeField2D<T, DESCRIPTOR, SCALAR> pL_lat( sLattice1 );
    SuperIdentity2D<T,T> p_lat( pL_lat );
    p_lat.getName() = "p_lat";

    // p_phs -> physical in physical units
    AnalyticalConst2D<T,T> _C_P( converter.getConversionFactorPressure() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> __C_P( _C_P, sLattice1 );
    SuperIdentity2D<T,T> p_phs( p_lat * __C_P );
    p_phs.getName() = "p_phs";

    // omega
    SuperLatticeField2D<T, DESCRIPTOR, OMEGA> omega( sLattice1 );
    SuperIdentity2D<T,T> _omega( omega );
    _omega.getName() = "Omega";

    vtmWriter.addFunctor( rho_lat );
    vtmWriter.addFunctor( rho_phs );
    vtmWriter.addFunctor( velocity_lat );
    vtmWriter.addFunctor( velocity_phs );
    vtmWriter.addFunctor( force_lat );
    vtmWriter.addFunctor( force_phs );
    vtmWriter.addFunctor( p_lat );
    vtmWriter.addFunctor( p_phs );
    vtmWriter.addFunctor( _omega );
    vtmWriter.write( iT );

  }

}


int main( int argc, char *argv[] )
{
  // Initialization for the parallel processing
  initialize( &argc, &argv );
  OstreamManager clout(std::cout, "main");

  // === 1st Step: Unit Converter ===
  MultiPhaseUnitConverterFromRelaxationTime<T,DESCRIPTOR> converter(
    (T)   Nx,                        // resolution
    (T)   0.7,                       // lattice relaxation time
    (T)   rho_liquid/1000,           // lattice density
    (T)   Lx,                        // charPhysLength: reference length of simulation geometry
    (T)   nu_liquid,                 // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   rho_liquid                 // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  // Instantiation of a cuboidDecomposition with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidDecomposition2D<T> cuboidDecomposition(0, converter.getPhysDeltaX(), {Nx, Ny}, noOfCuboids);
  // set periodic boundaries to the domain
  cuboidDecomposition.setPeriodicity({true, false});
  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  loadBalancer.print();
  // Instantiation of superGeometry
  SuperGeometry<T,2> superGeometry( cuboidDecomposition,loadBalancer, 3 );
  prepareGeometry( superGeometry );


  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  prepareLattice( sLattice, superGeometry, converter );

  SuperLatticeCoupling couplingInteractionForce(
    COUPLING{},
    names::Component1{}, sLattice
  );

  couplingInteractionForce.template setParameter<interaction::Polinomial::RHOV>(
                                                  rho_vapor/converter.getConversionFactorDensity());
  couplingInteractionForce.template setParameter<interaction::Polinomial::RHOL>(
                                                  rho_liquid/converter.getConversionFactorDensity());
  couplingInteractionForce.template setParameter<interaction::Polinomial::THICKNESS>(thickness);
  couplingInteractionForce.template setParameter<interaction::Polinomial::SURFTENSION>(
                                                  SurfTension/converter.getConversionFactorSurfaceTension());

  // Compute the interaction parameters
  interaction::Polinomial::computeParameters<T>(couplingInteractionForce);

  // Display the value of surface tension parameter
  // recommended kappaP not much larger than 1
  auto kappaP = couplingInteractionForce.template getParameter<interaction::Polinomial::KAPPAP>();
  clout << "Surface tension parameter: " << kappaP[0] << std::endl;

  initialCondition( sLattice, couplingInteractionForce );

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  std::cout << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();
  std::vector<T> output = {1./radius,0};
  for ( iT=0; iT<=maxIter; ++iT ) {
    getResults( sLattice, iT, superGeometry, timer, converter );

    // Collide and stream (and coupling) execution
    sLattice.collideAndStream();

    sLattice.executePostProcessors(stage::PreCoupling());
    sLattice.getCommunicator(stage::PreCoupling()).communicate();
    couplingInteractionForce.execute();

  }
  timer.stop();
  timer.printSummary();

  std::cout << "Finish" << std::endl;

  return 0;
}

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

/* isothermalEvaporation2d.cpp
 * Flat interface evaporates due to a lower pressure assigned
 * to the boundary
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

typedef double T;
using DESCRIPTOR = D2Q9<STATISTIC, // Store density field
                        CHEM_POTENTIAL, // Chemical potential for multiphase system
                        FORCE, // Thermodynamic forced based on chemical potential
                        EXTERNAL_FORCE, // Other forces such as gravity (optional in this case)
                        OMEGA, // Store density dependent relaxation frequency
                        SCALAR>; // Store laplacian of the density
using BulkDynamics = ForcedWagnerBGKdynamics<T, DESCRIPTOR>;
using COUPLING1 = ChemicalPotentialPostProcessor<EOS::Landau>;
using COUPLING2 = FreeEnergyForcedPostProcessor;

// Parameters in physical units
const T Lx = 400;     // domain length [m]
const T Ly = 10;      // domain height [m]
const T length = 300; // liquid phase length [m]

const T nu = 0.1; // interface viscosity [m2.s-1]

const T surfTension = 1.e-3; // surface tension [N.m-1]

const T rho_vapor = 1.; // vapor density [kg.m-3]
const T rho_liquid = 2.; // liquid density [kg.m-3]

// Boundary conditions in physical units
/**
 * Boundary conditions for: p_sat/p_c = 0.99,
 * we use the analytical solution to initialize
 * the system close to the final solution
 */
const T rho_V = 0.9809107952180959; // analytical solution vapor density
const T rho_L = 1.9926932521825274; // analytical solution liquid density
const T mu_b = -3.4274274827508194e-05; // boundary chemical potential
const T Ui_sol = 0.0021753492521563193; // interface velocity - analytical solution

// Resolution
const int Nx = 400;
const int Ny = 10;

// Parameters in lattice units
const T tau = 0.8; // relaxation time in lattice units
const T thickness = 5.; // interface thickness in lattice units; physThickness = delta_x * thickness

// Time parameters in lattice units
const int maxIter  = 30000;
const int vtkIter  = 100;
const int statIter = 100;


void prepareGeometry( SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  superGeometry.rename( 0,3 );
  superGeometry.rename( 3,2,{1,0} );
  superGeometry.rename( 2,1,{1,0} );
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
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);
  sLattice.defineDynamics<BulkDynamics>(superGeometry, 2);
  sLattice.defineDynamics<NoDynamicsWithFixedDensity>(superGeometry, 3);

  // Set boundary pressure condition
  boundary::set<boundary::LocalPressure<T, DESCRIPTOR, BulkDynamics>>(sLattice, superGeometry, 2);

  std::vector<T> zeroV = {0., 0.};
  AnalyticalConst2D<T,T> zeroVelocity( zeroV );

  // Computing the density
  AnalyticalConst2D<T,T> vapor ( rho_V/converter.getConversionFactorDensity() );
  SmoothIndicatorFactoredCuboid2D<T,T> liquid ( {Lx/2., Ly/2.}, length, Lx,
                                                 sqrt(2.)*thickness*converter.getConversionFactorLength(),
                                                 0, {0,0}, 0,
                                                 ( rho_L - rho_V )/converter.getConversionFactorDensity() );
  AnalyticalIdentity2D<T,T> rho ( vapor + liquid );

  // sign function
  AnalyticalConst2D<T,T> left ( -1. );
  SmoothIndicatorFactoredCuboid2D<T,T> right ( {Lx, Ly/2.}, Lx, Lx,
                                                sqrt(2.)*thickness*converter.getConversionFactorLength(),
                                                0, {0,0}, 0, 2. );
  AnalyticalIdentity2D<T,T> sign ( left + right );

  // Computing the velocity
  AnalyticalConst2D<T,T> _rho_L ( rho_L/converter.getConversionFactorDensity() );
  AnalyticalConst2D<T,T> _Ui ( Ui_sol/converter.getConversionFactorVelocity() );
  AnalyticalIdentity2D<T,T> solVelocity ( sign * _Ui * _rho_L / rho - sign * _Ui );

  sLattice.defineRhoU( superGeometry, 1, rho, solVelocity );
  sLattice.iniEquilibrium( superGeometry, 1, rho, solVelocity );
  sLattice.defineRhoU( superGeometry, 2, rho, solVelocity );
  sLattice.iniEquilibrium( superGeometry, 2, rho, solVelocity );
  sLattice.defineRhoU( superGeometry, 3, rho, solVelocity );
  sLattice.iniEquilibrium( superGeometry, 3, rho, solVelocity );

  // Set Chemical Potential
  AnalyticalConst2D<T,T> chemical( mu_b/converter.getConversionFactorChemicalPotential() );
  sLattice.defineField<descriptors::CHEM_POTENTIAL>( superGeometry, 2, chemical );
  sLattice.defineField<descriptors::CHEM_POTENTIAL>( superGeometry, 3, chemical );

  // Set OMEGA
  AnalyticalConst2D<T,T> omega( T( 1./converter.computeRelaxationTimefromPhysViscosity(nu) ) );
  sLattice.defineField<descriptors::OMEGA>( superGeometry, 1, omega );

  sLattice.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

  {
    auto& communicator = sLattice.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.exchangeRequests();
  }

  sLattice.initialize();
  sLattice.getCommunicator(stage::PreCoupling()).communicate();
  sLattice.executePostProcessors(stage::PreCoupling());

  clout << "Prepare Lattice ... OK" << std::endl;
}

/**
 * Correction for lattice velocity
 * Lattice velocity Ulat in this model is not the real velocity
 * Ureal = Ulat + 0.5*F
 * We compute Ulat = Ureal - 0.5*F as initial condition
 */
struct initialVelocityCorrection  {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    T U[descriptors::d<DESCRIPTOR>()];
    const auto force = cell.template getField<descriptors::FORCE>();

    cell.computeU(U);

    for(int iD = 0; iD < DESCRIPTOR::d; ++iD){
      U[iD] = U[iD] - 0.5*force[iD];
    }

    cell.defineU(U);

  }
};

template<typename COUPLING1, typename COUPLING2>
void initialCondition( SuperLattice<T, DESCRIPTOR>& sLattice,
                       SuperGeometry<T,2>& superGeometry,
                       COUPLING1& chemicalCoupling,
                       COUPLING2& forceCoupling )
{
  OstreamManager clout( std::cout,"initialCondition" );
  clout << "Initial Condition ..." << std::endl;


  // Compute force
  sLattice.executePostProcessors(stage::PreCoupling());
  sLattice.getCommunicator(stage::PreCoupling()).communicate();
  chemicalCoupling.execute();
  sLattice.executePostProcessors(stage::PreCoupling());
  sLattice.getCommunicator(stage::PreCoupling()).communicate();
  forceCoupling.execute();

  // Correct velocity
  sLattice.addPostProcessor<stage::PostCoupling>(meta::id<initialVelocityCorrection>());
  sLattice.executePostProcessors(stage::PostCoupling());

  clout << "Initial Condition ... OK" << std::endl;

}

void boundaryCondition( SuperLattice<T, DESCRIPTOR>& sLattice,
                        SuperGeometry<T,2>& superGeometry,
                        MultiPhaseUnitConverterFromRelaxationTime<T, DESCRIPTOR> const& converter )
{
  // Compute density at MN1
    SuperLatticeDensity2D<T, DESCRIPTOR> _density( sLattice );
    AnalyticalFfromSuperF2D<T,T> interpolRHO( _density, true, 1 );
    T dx = converter.getConversionFactorLength();
    T pos_MN1[2] = {dx, 0.};
    T rho_MN1;
    interpolRHO( &rho_MN1, pos_MN1 );

    // Set boundary density at MN2
    AnalyticalConst2D<T,T> densityV ( rho_V );
    AnalyticalConst2D<T,T> densityMN1 ( rho_MN1 );
    AnalyticalConst2D<T,T> half ( T(0.5) );
    AnalyticalIdentity2D<T,T> densityMN2 ( half * densityV + half * densityMN1 );
    AnalyticalConst2D<T,T> chemicalMN2 ( mu_b );

    sLattice.defineRho( superGeometry, 2, densityMN2 );
    sLattice.defineField<descriptors::CHEM_POTENTIAL>( superGeometry, 3, chemicalMN2 );

}

/**
 * Variable relaxation frequency increase stability and
 * decrease oscillations due to initialization
 * Bulk: omega = 1
 * Interface: omega = 1/tau
*/
void computeRelaxationFrequency( SuperLattice<T, DESCRIPTOR>& sLattice,
                                 SuperGeometry<T,2>& superGeometry,
                                 MultiPhaseUnitConverterFromRelaxationTime<T, DESCRIPTOR> const& converter )
{
    // Get density
    SuperLatticeDensity2D<T, DESCRIPTOR> _density( sLattice );

    // Find interface position
    T rhom = 0.5*(rho_V+rho_L)/converter.getConversionFactorDensity();
    T dx = converter.getConversionFactorLength();
    AnalyticalFfromSuperF2D<T,T> interpolRHO( _density, true, 1 );
    T pos1[2] = {0., 0.};
    T pos2[2] = {0., 0.};
    T h1 = 0, h2 = 0;
    for (int ix=1; ix<Nx-1; ix++) {
      pos1[0] = ix * dx;
      pos2[0] = ( ix + 1 ) * dx;
      T rho1, rho2;
      interpolRHO( &rho1, pos1 );
      interpolRHO( &rho2, pos2 );
      if ( rho1 < rhom && rho2 > rhom )
        h1 = ix * dx;
      if ( rho1 > rhom && rho2 < rhom )
        h2 = ix * dx;
    }

    // Computing indicator
    AnalyticalConst2D<T,T> bulk ( 1./tau );
    SmoothIndicatorFactoredCuboid2D<T,T> left ( {0., Ly/2.}, 2.*h1-10.*thickness*dx,
                      Lx, 1., 0, {0,0}, 0,
                      1. - 1./tau );
    SmoothIndicatorFactoredCuboid2D<T,T> center ( {Nx/2., Ly/2.}, h2-h1-10.*thickness*dx,
                      Lx, 1., 0, {0,0}, 0,
                      1. - 1./tau );
    SmoothIndicatorFactoredCuboid2D<T,T> right ( {Nx, Ly/2.}, 2.*h1-10.*thickness*dx,
                      Lx, 1., 0, {0,0}, 0,
                      1. - 1./tau );
    AnalyticalIdentity2D<T,T> _omega ( bulk + left + center + right );

    // Setting Omega field
    //SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> __omega( _omega, sLattice );
    sLattice.defineField<descriptors::OMEGA>( superGeometry, 1, _omega );
    sLattice.defineField<descriptors::OMEGA>( superGeometry, 2, _omega );

}

//std::vector<T>
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 int iT, SuperGeometry<T,2>& superGeometry,
                 MultiPhaseUnitConverterFromRelaxationTime<T, DESCRIPTOR> const& converter,
                 util::Timer<T>& timer )
{
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "isothermalEvaporation2D" );
  if ( iT==0 ) {
    // Writes the cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }
  // Get statistics
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT, iT );
  }

  //T delta_p = 0.;
  //T final_radius = 1.;
  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticeDensity2D<T, DESCRIPTOR> density_lat( sLattice );
    SuperIdentity2D<T,T> _density_lat( density_lat );
    _density_lat.getName() = "lattice density";

    AnalyticalConst2D<T,T> C_rho( converter.getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> _C_rho( C_rho, sLattice );
    SuperIdentity2D<T,T> _density_phs( density_lat * _C_rho );
    _density_phs.getName() = "physical density";

    SuperLatticeVelocity2D<T, DESCRIPTOR> velocityEq_lat( sLattice );
    SuperIdentity2D<T,T> _velocityEq_lat( velocityEq_lat );
    _velocityEq_lat.getName() = "equibibrium distribution lattice velocity";

    SuperLatticeField2D<T, DESCRIPTOR, FORCE> force_lat( sLattice );
    SuperIdentity2D<T,T> _force_lat( force_lat );
    _force_lat.getName() = "lattice force";

    AnalyticalConst2D<T,T> C_F( converter.getConversionFactorForce() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> _C_F( C_F, sLattice );
    SuperIdentity2D<T,T> _force_phs( force_lat * _C_F );
    _force_phs.getName() = "physical force";

    SuperLatticeExternalScalarField2D<T, DESCRIPTOR, CHEM_POTENTIAL> chemical_lat( sLattice );
    SuperIdentity2D<T,T> _chemical_lat( chemical_lat );
    _chemical_lat.getName() = "lattice chemical potential";

    AnalyticalConst2D<T,T> C_MU( converter.getConversionFactorChemicalPotential() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> _C_MU( C_MU, sLattice );
    SuperIdentity2D<T,T> _chemical_phs( chemical_lat * _C_MU );
    _chemical_phs.getName() = "physical chemical potential";

    SuperLatticeExternalScalarField2D<T, DESCRIPTOR, OMEGA> omega( sLattice );
    SuperIdentity2D<T,T> _omega( omega );
    _omega.getName() = "omega";

    AnalyticalConst2D<T,T> half( T(0.5) );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> _half( half, sLattice );
    SuperIdentity2D<T,T> _velocity_lat( velocityEq_lat + _half * force_lat );
    _velocity_lat.getName() = "lattice velocity";

    AnalyticalConst2D<T,T> C_U( converter.getConversionFactorVelocity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> _C_U( C_U, sLattice );
    SuperIdentity2D<T,T> _velocity_phs( _velocity_lat * _C_U );
    _velocity_phs.getName() = "physical velocity";

    vtmWriter.addFunctor( _density_lat );
    vtmWriter.addFunctor( _density_phs );
    vtmWriter.addFunctor( _force_lat );
    vtmWriter.addFunctor( _force_phs );
    vtmWriter.addFunctor( _chemical_lat );
    vtmWriter.addFunctor( _chemical_phs );
    vtmWriter.addFunctor( _omega );
    vtmWriter.addFunctor( _velocity_lat );
    vtmWriter.addFunctor( _velocity_phs );
    vtmWriter.write( iT );

    // Compute error in Interface velocity
    // first: compute vapor phase velocity
    AnalyticalFfromSuperF2D<T,T> interpolVelocity( _velocity_phs, true, 1 );
    T dx = converter.getConversionFactorLength();
    T pos[2] = {dx, 0.};
    T Ui_num[2] = {0., 0.};
    interpolVelocity( &Ui_num[0], pos );
    // Compute interface velocity from vapor phase velocity
    Ui_num[0] = rho_V * Ui_num[0] / ( rho_V - rho_L );
    clout << "Numerical interface velocity [m/s]: " << Ui_num[0] << std::endl;
    clout << "Analytical interface velocity [m/s]: " << Ui_sol << std::endl;
    clout << "Relative error [%]: " << 100. * abs( Ui_num[0] - Ui_sol ) / abs( Ui_sol ) << std::endl;

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
    (T)   tau,                       // lattice relaxation time
    (T)   rho_liquid,                // lattice density
    (T)   Lx,                        // charPhysLength: reference length of simulation geometry
    (T)   nu,                        // physViscosity: physical kinematic viscosity in __m^2 / s__
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
  CuboidDecomposition2D<T> cuboidDecomposition({0, 0}, converter.getPhysDeltaX(), {Nx, Ny}, noOfCuboids );
  // set periodic boundaries to the vertical direction
  cuboidDecomposition.setPeriodicity({ false, true });
  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  loadBalancer.print();
  // Instantiation of superGeometry
  SuperGeometry<T,2> superGeometry( cuboidDecomposition,loadBalancer, 3 ); // checar esse 2
  prepareGeometry( superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  prepareLattice( sLattice, superGeometry, converter );

  // === 4th Step: Prepare Couplings ===
  // Chemical potential coupling
  SuperLatticeCoupling chemicalCoupling(
    COUPLING1{},
    names::Component1{}, sLattice
  );
  chemicalCoupling.restrictTo(superGeometry.getMaterialIndicator({1,2}));

  chemicalCoupling.template setParameter<EOS::Landau::RHOV>(
                  rho_vapor/converter.getConversionFactorDensity());
  chemicalCoupling.template setParameter<EOS::Landau::RHOL>(
                  rho_liquid/converter.getConversionFactorDensity());
  chemicalCoupling.template setParameter<EOS::Landau::THICKNESS>(thickness);
  chemicalCoupling.template setParameter<EOS::Landau::SURFTENSION>(
                  surfTension/converter.getConversionFactorSurfaceTension());

  // Compute the EOS parameters
  EOS::Landau::computeParameters<T>(chemicalCoupling);

  // Force coupling
  SuperLatticeCoupling forceCoupling(
    COUPLING2{},
    names::Component1{}, sLattice
  );
  forceCoupling.restrictTo(superGeometry.getMaterialIndicator({1,2}));

  // Initial Condition
  initialCondition( sLattice, superGeometry, chemicalCoupling, forceCoupling );

  // === 5th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();
  std::vector<T> output = {1./length,0};
  for ( iT=0; iT<=maxIter; ++iT ) {

    // Apply boundary conditions
    boundaryCondition( sLattice, superGeometry, converter );

    // Compute relaxation frequency
    computeRelaxationFrequency( sLattice, superGeometry, converter );

    // Computation and output of the results
    getResults( sLattice, iT, superGeometry, converter, timer );
    if ( std::isnan( sLattice.getStatistics().getAverageEnergy() ) ) {
      clout << "unstable!" << std::endl;
      break;
    }

    // Collide and stream (and coupling) execution
    sLattice.collideAndStream();

    // Compute chemical potential
    sLattice.executePostProcessors(stage::PreCoupling());
    sLattice.getCommunicator(stage::PreCoupling()).communicate();
    chemicalCoupling.execute();

    // Compute thermodynamic force
    sLattice.executePostProcessors(stage::PreCoupling());
    sLattice.getCommunicator(stage::PreCoupling()).communicate();
    forceCoupling.execute();

  }
  timer.stop();
  timer.printSummary();

  clout << "Finish" << std::endl;

  return 0;
}

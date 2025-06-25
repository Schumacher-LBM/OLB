/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2023 Tim Bingert
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

/* waterAirFlatInterface2d.cpp
 * Benchmark for thermodynamic consistency of the multi-component-
 * multi-phase Shan-Chen model at standard conditions.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;

using DESCRIPTOR = D2Q9<VELOCITY,FORCE,EXTERNAL_FORCE,STATISTIC,SCALAR,PSI>;
using BulkDynamics = MultiComponentForcedBGKdynamics<T, DESCRIPTOR>;

constexpr unsigned N_COMPONENTS = 3;
using COUPLING = MCMPForcedPostProcessor<N_COMPONENTS>;
using STATISTICS = RhoPsiStatistics<interaction::MCPRpseudoPotential<N_COMPONENTS>,N_COMPONENTS>;

// Parameters for the simulation setup
const int N = 4;                                                // domain resolution
const int heightFactor = 25;
int phaseLength = 50;                                           // [lattice units]
const T L_char = 1e-6;                                          // charPhysLength [physical units]
const T Re = 0.;                                                // definition: Reynolds number of continuous phase
const T tau_nuH2O = 0.6;                                        // relaxation time H2O lattice [lattice units]
const T tau_nuAir = (tau_nuH2O-0.5)*15.32+0.5;//4.415;          // relaxation time N2+O2 lattice [lattice units]
const T viscosityH2O = 1.0e-6;                                  // physViscosity H2O lattice [physical units]
const T viscosityAir = viscosityH2O*15.32;                      // physViscosity N2+O2 lattice [physical units]
const T g = 0;                                                  // gravitational acceleration [physical units]
T a_0L = 3./245.;                                               // tune for stability/accuracy
T initEpsilon = 2.2355;                                         // tune for chemical potential equilibrium in phases
const T pressure = 1.013e5;                                     // For Vapor-Liquid-Equilibrium [physical units] Pascal
const T temperature = 298.15;                                   // For Vapor-Liquid-Equilibrium [physical units] Kelvin
const T surfaceTension = 0.07;                                  // For unit conversion [physical units], not crucial for phase composition

/*==============================================================
 * go to src/dynamics/shanChenForcedPostProcessor.h and enable *
 * the THIRD_COMPONENT by uncommenting #define THIRD_COMPONENT *
 ==============================================================*/
//H2O, N2, O2
const std::vector<T> z = {0.99, 0.0079, 0.0021};                // feed composition for Vapor-Liquid-Equilibrium as molar fraction
const std::vector<T> a = {0.5995808, 0.1480650, 0.1506765};
const std::vector<T> b = {1.8955853e-5, 2.4010114e-5, 1.9893672e-5};
const std::vector<T> M = {0.01802, 0.02801, 0.03200};
const std::vector<T> T_c = {647.3, 126.2, 155.0};
const std::vector<T> p_c{22089000, 3400000, 5040000};
const std::vector<T> omega{0.34, 0.0377, 0.025};
const std::vector<T> devi = {0.867805648, 0.432399567, 0.41302780};

//H2OH2O, H2ON2, H2OO2, N2H2O, N2N2, N2O2, O2H2O, O2N2, O2O2
const std::vector<T> alpha = {0., 0.199222317, 0.193233601, 0.199222317, 0., 0., 0.193233601, 0., 0.};
const std::vector<T> g_I = {0., -4.29088111e3, -1.95640777e2, 3.02126911e4, 0., 5.65244934e2, 6.13396078e4, -5.01392189e2, 0.};
const std::vector<T> g_II = {0., 3.47847412e1, 2.10776021e1, -3.70834075e1, 0., 0., -1.22744109e2, 0., 0.};
std::vector<T> rho0L(3), rho0V(3);
const int maxIter  = 250000;
const int vtkIter  = 1000;
const int statIter = 1000;

void prepareGeometry( SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  superGeometry.rename( 0,1 );
  // bulk, MN=1

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

template <typename SuperLatticeCoupling, typename StatisticsCoupling>
void prepareLattice( SuperLattice<T, DESCRIPTOR>& sLattice1,
                     SuperLattice<T, DESCRIPTOR>& sLattice2,
                     SuperLattice<T, DESCRIPTOR>& sLattice3,
                     SuperLatticeCoupling& coupling,
                     StatisticsCoupling& statistics,
                     MultiPhaseUnitConverter<T, DESCRIPTOR> const& converter,
                     SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;
  clout << "epsilon: " << initEpsilon << std::endl;
  // define lattice Dynamics
  sLattice1.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice2.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice3.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice1.defineDynamics<BulkDynamics>(superGeometry, 1);
  sLattice2.defineDynamics<BulkDynamics>(superGeometry, 1);
  sLattice3.defineDynamics<BulkDynamics>(superGeometry, 1);

  //thermodynamic initial conditions in lattice units
  T p_L = pressure/converter.getConversionFactorPressure();
  T T_L = temperature/converter.getConversionFactorTemperature();
  int componentNumber = z.size(); //maximal number of 10 components!
  std::vector<T> a_L(10), b_L(10), M_L(10), Tc_L(10), pc_L(10), omega_L(10), devi_L(10), alpha_L(100), gI_L(100), gII_L(100);
  for(int i = 0; i < componentNumber; i++){
    a_L[i]     = a[i]/converter.getConversionFactorEoSa();
    b_L[i]     = b[i]/converter.getConversionFactorEoSb();
    M_L[i]     = M[i]/converter.getConversionFactorMolarMass();
    Tc_L[i]    = T_c[i]/converter.getConversionFactorTemperature();
    pc_L[i]    = p_c[i]/converter.getConversionFactorPressure();
    omega_L[i] = omega[i];
    devi_L[i]  = devi[i];
  }
  T C_R = 8.314462618;
  T C_temp = converter.getConversionFactorTemperature();
  for(int i = 0; i < componentNumber*componentNumber; i++){
    alpha_L[i] = alpha[i];
    gI_L[i]    = g_I[i]/(C_R*C_temp);
    gII_L[i]   = g_II[i]/C_R;
  }

  MultiComponentPengRobinson VLE(p_L, T_L, z, a_L, b_L, M_L, Tc_L, pc_L, omega, devi, alpha, gI_L, gII_L);
  T beta0 = (z[1]+z[2]);
  std::vector<T> vxVLE = VLE.iterate_VLE(1e-11, beta0);   // Molar volumes and fractions of equilibrium phases [lattice units]
  clout <<"VLE: "<<vxVLE[0]<<", "<<vxVLE[1]<<", "<<vxVLE[2]<<", "<<vxVLE[3]<<", "<<vxVLE[4]<<", "<<vxVLE[5]<<", "<<vxVLE[6]<<", "<<vxVLE[7]<<std::endl;
  std::vector<T> chi = VLE.getChis(3);                             // Force split factor from VLE
  clout <<"Chis: "<<chi[0]<<", "<<chi[1]<<", "<<chi[2]<<std::endl;

  // bulk initial conditions
  std::vector<T> v = {0., 0.};
  AnalyticalConst2D<T,T> zeroVelocity( v );

  AnalyticalConst2D<T,T> liquidH2O ( 1./vxVLE[0]*vxVLE[2]*M_L[0] );
  AnalyticalConst2D<T,T> liquidN2  ( 1./vxVLE[0]*vxVLE[3]*M_L[1] );
  AnalyticalConst2D<T,T> liquidO2  ( 1./vxVLE[0]*vxVLE[4]*M_L[2] );

  SmoothIndicatorFactoredCuboid2D<T,T> vaporH2O( {N/2., 100./2.}, 0, phaseLength, 4., 0, {0,0}, 0, (1./vxVLE[1]*vxVLE[5]*M_L[0] - 1./vxVLE[0]*vxVLE[2]*M_L[0]) );
  SmoothIndicatorFactoredCuboid2D<T,T> vaporN2 ( {N/2., 100./2.}, 0, phaseLength, 4., 0, {0,0}, 0, (1./vxVLE[1]*vxVLE[6]*M_L[1] - 1./vxVLE[0]*vxVLE[3]*M_L[1]) );
  SmoothIndicatorFactoredCuboid2D<T,T> vaporO2 ( {N/2., 100./2.}, 0, phaseLength, 4., 0, {0,0}, 0, (1./vxVLE[1]*vxVLE[7]*M_L[2] - 1./vxVLE[0]*vxVLE[4]*M_L[2]) );

  AnalyticalIdentity2D<T,T> rhoH2O( liquidH2O + vaporH2O );
  AnalyticalIdentity2D<T,T> rhoN2( liquidN2 + vaporN2 );
  AnalyticalIdentity2D<T,T> rhoO2( liquidO2 + vaporO2 );

  for(int i = 0; i < (int)z.size(); i++){
    rho0L[i] = 1./vxVLE[0]*vxVLE[i+2]*M_L[i]*converter.getConversionFactorDensity();
    rho0V[i] = 1./vxVLE[1]*vxVLE[i+5]*M_L[i]*converter.getConversionFactorDensity();
  }

  sLattice1.defineRhoU( superGeometry, 1, rhoH2O, zeroVelocity );
  sLattice2.defineRhoU( superGeometry, 1, rhoN2, zeroVelocity );
  sLattice3.defineRhoU( superGeometry, 1, rhoO2, zeroVelocity );

  sLattice1.iniEquilibrium( superGeometry, 1, rhoH2O, zeroVelocity );
  sLattice2.iniEquilibrium( superGeometry, 1, rhoN2, zeroVelocity );
  sLattice3.iniEquilibrium( superGeometry, 1, rhoO2, zeroVelocity );

  // Only works properly after implementing walls at the top and bottom
  std::vector<T> F( 2,T() );
  F[1] = -g/(converter.getPhysDeltaX()/converter.getPhysDeltaT()/converter.getPhysDeltaT());
  clout << "Gravitational acceleration in lattice units: " << g/(converter.getPhysDeltaX()/converter.getPhysDeltaT()/converter.getPhysDeltaT()) << std::endl;
  AnalyticalConst2D<T,T> f( F );
  sLattice1.defineField<descriptors::EXTERNAL_FORCE>( superGeometry, 1, f );
  sLattice2.defineField<descriptors::EXTERNAL_FORCE>( superGeometry, 1, f );
  sLattice3.defineField<descriptors::EXTERNAL_FORCE>( superGeometry, 1, f );

  sLattice1.setParameter<descriptors::OMEGA>( 1./tau_nuH2O );
  sLattice2.setParameter<descriptors::OMEGA>( 1./tau_nuAir );
  sLattice3.setParameter<descriptors::OMEGA>( 1./tau_nuAir );

  T sigma = converter.getLatticeSurfaceTension()*(1.);
  clout << "Sigma for correct unit conversion: " << sigma << std::endl;
  coupling.template setParameter<COUPLING::CHI>(chi);
  coupling.template setParameter<COUPLING::G>(-1.);
  coupling.template setParameter<COUPLING::SIGMA>(sigma);
  coupling.template setParameter<COUPLING::EPSILON>(initEpsilon);

  statistics.template setParameter<STATISTICS::TEMPERATURE>(T_L);
  statistics.template setParameter<STATISTICS::G>(-1.);
  statistics.template setParameter<STATISTICS::K>(1.);
  statistics.template setParameter<STATISTICS::A>(a_L);
  statistics.template setParameter<STATISTICS::B>(b_L);
  statistics.template setParameter<STATISTICS::MM>(M_L);
  statistics.template setParameter<STATISTICS::TCRIT>(Tc_L);
  statistics.template setParameter<STATISTICS::DEVI>(devi);
  statistics.template setParameter<STATISTICS::ALPHA>(alpha);
  statistics.template setParameter<STATISTICS::GI>(gI_L);
  statistics.template setParameter<STATISTICS::GII>(gII_L);

  {
    auto& communicator = sLattice1.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.requestField<PSI>();
    communicator.exchangeRequests();
  }

  {
    auto& communicator = sLattice2.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.exchangeRequests();
  }

  {
    auto& communicator = sLattice3.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.exchangeRequests();
  }

  sLattice1.initialize();
  sLattice2.initialize();
  sLattice3.initialize();
  statistics.execute();

  clout << "Prepare Lattice ... OK" << std::endl;
}

T getResults( SuperLattice<T, DESCRIPTOR>& sLattice1,
              SuperLattice<T, DESCRIPTOR>& sLattice2,
              SuperLattice<T, DESCRIPTOR>& sLattice3,
              int iT, SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer,
              MultiPhaseUnitConverter<T, DESCRIPTOR> const& converter)
{
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "waterAirFlatInterface2d" );
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
    sLattice1.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLattice2.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLattice3.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  T error = 0;
  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    sLattice1.setProcessingContext(ProcessingContext::Evaluation);
    sLattice2.setProcessingContext(ProcessingContext::Evaluation);
    sLattice3.setProcessingContext(ProcessingContext::Evaluation);

    //Factors for conversion back to physical units
    AnalyticalConst2D<T,T> _C_rho( converter.getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> C_rho(_C_rho, sLattice1);
    AnalyticalConst2D<T,T> _C_u( converter.getConversionFactorVelocity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> C_u(_C_u, sLattice2);
    AnalyticalConst2D<T,T> _C_p( converter.getConversionFactorPressure() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> C_p(_C_p, sLattice1);

    SuperLatticeDensity2D<T, DESCRIPTOR> density1L( sLattice1 );
    SuperIdentity2D<T,T> density1( C_rho*density1L );
    density1.getName() = "rhoH2O";
    SuperLatticeDensity2D<T, DESCRIPTOR> density2L( sLattice2 );
    SuperIdentity2D<T,T> density2( C_rho*density2L );
    density2.getName() = "rhoN2";
    SuperLatticeDensity2D<T, DESCRIPTOR> density3L( sLattice3 );
    SuperIdentity2D<T,T> density3( C_rho*density3L );
    density3.getName() = "rhoO2";
    SuperIdentity2D<T,T> density( density1+density2+density3 );
    density.getName() = "rho";

    SuperLatticeVelocity2D<T, DESCRIPTOR> velocityL( sLattice1 );
    SuperIdentity2D<T,T> velocity( C_u*velocityL );
    velocity.getName() = "velocity";

    SuperLatticeExternalScalarField2D<T, DESCRIPTOR, SCALAR> bulkPressureL( sLattice1 );
    SuperIdentity2D<T,T> bulkPressure( C_p*bulkPressureL );
    bulkPressure.getName() = "bulkPressure";

    vtmWriter.addFunctor( density1 );
    vtmWriter.addFunctor( density2 );
    vtmWriter.addFunctor( density3 );
    vtmWriter.addFunctor( density );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( bulkPressure );
    vtmWriter.write( iT );

    // calculate total error of phase compositions
    AnalyticalFfromSuperF2D<T, T> interpolRho1( density1, true, 1);
    AnalyticalFfromSuperF2D<T, T> interpolRho2( density2, true, 1);
    AnalyticalFfromSuperF2D<T, T> interpolRho3( density3, true, 1);
    double position[2] = { 0.5*N, 0.5*N*heightFactor+1 };
    double liquidDensities[3] = {0.,0.,0.};
    double vaporDensities[3] = {0.,0.,0.};
    interpolRho1(&vaporDensities[0], position);
    interpolRho2(&vaporDensities[1], position);
    interpolRho3(&vaporDensities[2], position);
    position[1] = 1.;
    interpolRho1(&liquidDensities[0], position);
    interpolRho2(&liquidDensities[1], position);
    interpolRho3(&liquidDensities[2], position);
    double liquidError[3] = {abs(liquidDensities[0]-rho0L[0])/rho0L[0],
                             abs(liquidDensities[1]-rho0L[1])/rho0L[1],
                             abs(liquidDensities[2]-rho0L[2])/rho0L[2]};
    double vaporError[3] = {abs(vaporDensities[0]-rho0V[0])/rho0V[0],
                            abs(vaporDensities[1]-rho0V[1])/rho0V[1],
                            abs(vaporDensities[2]-rho0V[2])/rho0V[2]};
    for(int i = 0; i < (int)z.size(); i++){
      error += vaporError[i]*vaporError[i] + liquidError[i]*liquidError[i];
    }
    error = util::pow(error/6.,0.5);
    clout << "Average error of all component densities: " << std::setprecision(6) << error << std::endl;
  }
  return error;
}

std::vector<T> simulate()
{
  // === 1st Step: Initialization ===
  OstreamManager clout( std::cout,"main" );
  MultiPhaseUnitConverter<T,DESCRIPTOR> const converter(
    int   {phaseLength},            // resolution
    (T)   L_char,                   // charPhysLength: reference length of simulation geometry
    (T)   Re/L_char*viscosityH2O,   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   viscosityH2O,             // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   a[0],                     // physEoSa: H2O energy parameter in __kg m^5 / mol^2 s^2__
    (T)   a_0L,                     // latticeEoSa: first component's energy parameter in lattice units
    (T)   b[0],                     // physEoSb: H2O co-volume parameter in __m^3 / mol__
    (T)   M[0],                     // physMolarMass: H2O molar mass for EoS in __kg / mol__
    (T)   surfaceTension,           // physSurfaceTension: physical surface tension of mixture in __kg / s^2__
    (T)   temperature,              // charPhysTemperature: temperature of VLE in __K__
    (T)   pressure                  // charPhysPressure: pressure of VLE in __kg / m s^2__
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
  CuboidDecomposition2D<T> cuboidDecomposition({0, 0}, 1, {N, heightFactor*N}, noOfCuboids);
  // set periodic boundaries to the domain
  cuboidDecomposition.setPeriodicity({ true, true });
  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  loadBalancer.print();
  // Instantiation of superGeometry
  SuperGeometry<T,2> superGeometry( cuboidDecomposition,loadBalancer,2 );
  prepareGeometry( superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice1( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLattice2( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLattice3( superGeometry );
  SuperLatticeCoupling coupling(
    COUPLING{},
    names::Component1{}, sLattice1,
    names::Component2{}, sLattice2,
    names::Component3{}, sLattice3);

  SuperLatticeCoupling statistics(
    STATISTICS{},
    names::Component1{}, sLattice1,
    names::Component2{}, sLattice2,
    names::Component3{}, sLattice3);

  prepareLattice( sLattice1, sLattice2, sLattice3, coupling, statistics, converter, superGeometry );

  sLattice1.writeSummary("lattice1");
  sLattice2.writeSummary("lattice2");
  sLattice3.writeSummary("lattice3");

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();
  std::vector<T> output = {initEpsilon,0};
  for ( iT=0; iT<=maxIter; ++iT ) {
    // Collide and stream (and coupling) execution
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();
    sLattice3.collideAndStream();

    statistics.execute();
    sLattice1.getCommunicator(stage::PreCoupling()).communicate();
    sLattice2.getCommunicator(stage::PreCoupling()).communicate();
    sLattice3.getCommunicator(stage::PreCoupling()).communicate();
    coupling.execute();
    // Computation and output of the results
    output[1] = getResults( sLattice1, sLattice2, sLattice3, iT, superGeometry, timer, converter );
    if ( std::isnan( sLattice1.getStatistics().getAverageEnergy() ) ) {
      break;
    }
  }
  timer.stop();
  timer.printSummary();
  return output;
}

int main( int argc, char *argv[] )
{
  initialize( &argc, &argv );
  if (argc > 1) {
    initEpsilon = atof(argv[1]);
  }
  if (argc > 2) {
    a_0L = atof(argv[2]);
  }

  singleton::directories().setOutputDir( "./tmp/" );
  std::vector<T> out = simulate();
  std::cout << "Final average error of all component densities: " << out[1] << std::endl;
}

/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Tim Bingert, Michael Rennick
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

/** flatInterfaceEOC2d.cpp
 * In this example a flat interface test is performed for several multi-
 * phase models such as incompressible Allen-Cahn or well-balanced Cahn-
 * Hilliard. A rectangular domain of fluid 2 is immersed in fluid 1.
 * A diffusive interface forms with the profile of a hyperbolic tangent
 * whose accuracy is measured for multiple resolutions in order to test
 * the models experimental order of convergence. The equilibrium pressure
 * is also investigated in a similar way. The grid refinement can be
 * performed with either a constant or a decreasing Cahn number. Multiple
 * different error norms can then be written into an output EOC.dat file.
 *
 * This example shows the simplest application of the hybrid (and local)
 * phase field Allen-Cahn model with periodic boundaries, based on:
 *
 * Liu, Xi, Zhenhua Chai, and Baochang Shi. "Improved hybrid Allen-Cahn
 * phase-field-based lattice Boltzmann method for incompressible two-phase
 * flows." Physical Review E 107.3 (2023): 035308.
 *
 * It also shows the same application of the well-balanced Cahn-Hilliard
 * model that has no spurious currents (machine accuracy), from:
 *
 * Ju, Long, et al. "A well-balanced lattice Boltzmann model for binary
 * fluids based on the incompressible phase-field theory." arXiv preprint
 * arXiv:2311.10827 (2023).
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

/// Available models
#define ALLEN_CAHN
#ifdef ALLEN_CAHN
#define LOCAL
//#define HYBRID
#endif
//#define CAHN_HILLIARD

using T = FLOATING_POINT_TYPE;
using NSDESCRIPTOR = D2Q9<RHO,NABLARHO,FORCE,EXTERNAL_FORCE,TAU_EFF,STATISTIC>;
using NSBulkDynamics = MPIncBGKdynamics<T,NSDESCRIPTOR>;
#ifdef CAHN_HILLIARD
using PFDESCRIPTOR = D2Q9<FORCE,SOURCE,SOURCE_OLD,VELOCITY,STATISTIC,CHEM_POTENTIAL,PHIWETTING>;
using PFBulkDynamics = WellBalancedCahnHilliardBGKdynamics<T,PFDESCRIPTOR>;
using Coupling = WellBalancedCahnHilliardPostProcessor;
#endif
#ifdef ALLEN_CAHN
  #ifdef LOCAL
using PFDESCRIPTOR = D2Q9<FORCE,SOURCE,VELOCITY,OLD_PHIU,STATISTIC>;
using Coupling = LiangPostProcessor;
  #endif
  #ifdef HYBRID
using PFDESCRIPTOR = D2Q9<FORCE,SOURCE,SOURCE_OLD,VELOCITY,OLD_PHIU,STATISTIC,PSI,NORMGRADPSI,SCALAR,PSI0,TOP,BOTTOM>;
using Coupling = AllenCahnPostProcessor;
using Helper = AllenCahnNonLocalHelper;
  #endif
using PFBulkDynamics = AllenCahnBGKdynamics<T,PFDESCRIPTOR>;
#endif

// Parameters for the simulation setup
const int Nx = 2;                                        // domain resolution x
int Ny = 100;                                            // domain resolution y
int phaseLength = Ny/2;                                  // [lattice units]
const T charPhysLength = 100e-6;                         // charPhysLength [physical units]
const T Re = 0.;                                         // definition: Reynolds number of continuous phase
const T surfaceTension = 0.072;                          // surface tension [physical units]
const T phys_pressure = 1e5;                             // system pressure [physical units]
const T viscosityH2O = 9e-7;                             // physViscosity H2O liquid [physical units]
const T densityH2O = 1e3;                                // physDensity H2O liquid [physical units]
const T tau_l = 0.8;                                     // relaxation time Water lattice [lattice units]
const T tau_v = 10.*(tau_l-0.5)+0.5;                     // relaxation time Air lattice [lattice units]
const T tau_mobil = 0.8;                                 // relaxation time for interface mobility [lattice units]
T sigma = 0.01;                                          // surface tension [lattice units]
T w = 5.;                                                // interface thickness [lattice units]
std::vector<T> rhos = {0.001, 1.};
T DeltaRho = densityH2O/rhos[1];
const int maxIter  = 100000000;
const int vtkIter  = 20000;
const int statIter = 20000;

// variables for eoc analysis
const bool Cahn_const = true;
T orderParameterL1RelError = 0;
T orderParameterL2RelError = 0;
T orderParameterLinfRelError = 0;
T pressureL1RelError = 0;
T pressureL2RelError = 0;
T pressureLinfRelError = 0;
T velocityLinfAbsError = 0;

void prepareGeometry( SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  superGeometry.rename( 0,1 );
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

#ifdef HYBRID
template <typename STAGE>
void signedDistanceFunction( SuperLattice<T,PFDESCRIPTOR>& sLatticeAC,
                             SuperGeometry<T,2>& superGeometry, T w )
{
  T max = 1.;
  while ( max <= 3*w ) {
    int in[2];
    sLatticeAC.getCommunicator(STAGE{}).communicate();
    sLatticeAC.executePostProcessors(STAGE{});
    SuperLatticeExternalScalarField2D<T, PFDESCRIPTOR, PSI> psi( sLatticeAC );
    SuperMax2D<T,T> Max_psi_(psi, superGeometry, 1);
    Max_psi_(&max, in);
  }
}
#endif

template <typename SuperLatticeCoupling>
void prepareLattice( SuperLattice<T,NSDESCRIPTOR>& sLatticeNS,
                     SuperLattice<T,PFDESCRIPTOR>& sLatticePF,
                     SuperLatticeCoupling& coupling,
                     UnitConverter<T,NSDESCRIPTOR> const& converter,
                     SuperGeometry<T,2>& superGeometry, int Ny, int phaseLength, T w, T sigma )
{
  OstreamManager clout( std::cout,"prepareLattice" );
#ifdef HYBRID
  clout << "Prepare Lattice with hybrid Allen-Cahn model ..." << std::endl;
#endif
#ifdef LOCAL
  clout << "Prepare Lattice with local Allen-Cahn model ..." << std::endl;
#endif
#ifdef CAHN_HILLIARD
  clout << "Prepare Lattice with well-balanced Cahn-Hilliard model ..." << std::endl;
#endif

  // define lattice Dynamics
  sLatticeNS.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticeNS.defineDynamics<NSBulkDynamics>(superGeometry, 1);
  sLatticePF.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticePF.defineDynamics<PFBulkDynamics>(superGeometry, 1);

  // bulk initial conditions
  T deltaX = converter.getPhysDeltaX();
  Vector<T,2> u(0., 0.);
  AnalyticalConst2D<T,T> zeroVelocity( u );

  AnalyticalConst2D<T,T> one ( 1. );
  AnalyticalConst2D<T,T> zero ( 0. );
  T lattice_pressure = phys_pressure/converter.getConversionFactorPressure();
  clout << "lattice pressure: " << lattice_pressure << std::endl;
  AnalyticalConst2D<T,T> pres ( lattice_pressure );
  AnalyticalConst2D<T,T> rhov ( rhos[0] );
  AnalyticalConst2D<T,T> rhol ( rhos[1] );
  AnalyticalConst2D<T,T> tauv ( tau_v );
  AnalyticalConst2D<T,T> taul ( tau_l );

  SmoothIndicatorFactoredCuboid2D<T,T> interfaceAC( {Nx/2.*deltaX, (Ny/2.+0.5)*deltaX}, 0, phaseLength*deltaX, w*deltaX/2, 0, {0,0}, 0, -1. );
  AnalyticalIdentity2D<T,T> phi( one + interfaceAC );
  AnalyticalIdentity2D<T,T> rho( rhov + (rhol-rhov)*phi );
  AnalyticalIdentity2D<T,T> tau( tauv + (taul-tauv)*phi );
  AnalyticalIdentity2D<T,T> pressure( pres );

  sLatticeNS.defineField<descriptors::RHO>( superGeometry, 1, rho );
  sLatticeNS.defineField<descriptors::TAU_EFF>( superGeometry, 1, tau );
#ifdef ALLEN_CAHN
  sLatticePF.defineField<descriptors::OLD_PHIU>( superGeometry, 1, zeroVelocity );
#endif
#ifdef CAHN_HILLIARD
  sLatticePF.defineField<descriptors::CHEM_POTENTIAL>( superGeometry, 1, zero );
#endif

  sLatticeNS.defineRhoU( superGeometry, 1, pressure, zeroVelocity );
  sLatticeNS.iniEquilibrium( superGeometry, 1, pressure, zeroVelocity );
  sLatticePF.defineRhoU( superGeometry, 1, phi, zeroVelocity );
  sLatticePF.iniEquilibrium( superGeometry, 1, phi, zeroVelocity );

#ifdef ALLEN_CAHN
  sLatticePF.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());
  #ifdef HYBRID
  sLatticePF.addPostProcessor<stage::PostPostProcess>(meta::id<Helper>{});
  sLatticePF.addPostProcessor<stage::PreCoupling>(meta::id<initialPsi>{});
  sLatticePF.addPostProcessor<stage::IterativePostProcess>(meta::id<normGradPsi>{});
  sLatticePF.addPostProcessor<stage::IterativePostProcess>(meta::id<psiEvolve>{});
  sLatticePF.setParameter<psiEvolve::DELTAT>(0.5);
  sLatticePF.setParameter<descriptors::EPSILON>( 3.0*w );
  {
    auto& communicator = sLatticePF.getCommunicator(stage::IterativePostProcess());
    communicator.requestOverlap(1);
    communicator.requestField<PSI>();
    communicator.exchangeRequests();
  }
  coupling.template setParameter<Coupling::TAUS>({tau_v,tau_l,tau_mobil});
  #else
  coupling.template setParameter<Coupling::TAUS>({tau_v,tau_l});
  #endif
  coupling.template setParameter<Coupling::SIGMA>(sigma);
  coupling.template setParameter<Coupling::W>(w);
  coupling.template setParameter<Coupling::RHOS>(rhos);
#endif
#ifdef CAHN_HILLIARD
  coupling.template setParameter<Coupling::TAUS>({tau_v,tau_l});
  coupling.template setParameter<Coupling::RHOS>(rhos);

  sLatticePF.addPostProcessor<stage::PreCoupling>(meta::id<RhoWettingStatistics>());
  sLatticePF.addPostProcessor<stage::ChemPotCalc>(meta::id<ChemPotentialPhaseFieldProcessor>());
  sLatticePF.setParameter<descriptors::SCALAR>( sigma );
#endif

  sLatticeNS.setParameter<descriptors::OMEGA>( 1./tau_l );
  sLatticePF.template setParameter<descriptors::OMEGA>( 1./tau_mobil );
  sLatticePF.template setParameter<descriptors::INTERFACE_WIDTH>( w );

  {
    auto& communicator = sLatticePF.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.template requestField<STATISTIC>();
#ifdef HYBRID
    communicator.template requestField<PSI>();
#endif
#ifdef CAHN_HILLIARD
    communicator.template requestField<CHEM_POTENTIAL>();
#endif
    communicator.exchangeRequests();
  }

  sLatticePF.executePostProcessors(stage::PreCoupling());
#ifdef HYBRID
  signedDistanceFunction<stage::IterativePostProcess>(sLatticePF,superGeometry,w);
#endif
#ifdef CAHN_HILLIARD
  sLatticePF.getCommunicator(stage::PreCoupling()).communicate();
  sLatticePF.executePostProcessors(stage::ChemPotCalc());
#endif
  sLatticeNS.initialize();
  sLatticePF.initialize();
  sLatticePF.getCommunicator(stage::PreCoupling()).communicate();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void error( SuperGeometry<T,2>& superGeometry,
            SuperLattice<T, PFDESCRIPTOR>& sLatticePF,
            SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
            UnitConverter<T,NSDESCRIPTOR> const& converter,
            int Nx, int Ny, int phaseLength, T w )
{
  OstreamManager clout( std::cout,"error" );
  int tmp[]= { };
  T result[2]= { };

  T deltaX = converter.getPhysDeltaX();
  AnalyticalConst2D<T,T> p_phys ( phys_pressure );
  AnalyticalConst2D<T,T> zero ( 0. );
  AnalyticalConst2D<T,T> one ( 1. );
  Vector<T,2> u(0., 0.);
  AnalyticalConst2D<T,T> zeroVelocity( u );
  SmoothIndicatorFactoredCuboid2D<T,T> interface_diffuse( {Nx/2.*deltaX, (Ny/2.+0.5)*deltaX}, 0, phaseLength*deltaX, w*deltaX/2, 0, {0,0}, 0, -1. );
  AnalyticalIdentity2D<T,T> phi0_diffuse( one + interface_diffuse );
  Vector<T,2> extend( 2./100.*charPhysLength, 0.25*charPhysLength );
  Vector<T,2> origin1( 0., 0.);
  Vector<T,2> origin2( 0., 3./4.*charPhysLength);
  IndicatorCuboid2D<T> ind1( extend, origin1 );
  IndicatorCuboid2D<T> ind2( extend, origin2 );
  SmoothIndicatorCuboid2D<T,T> interface_sharp1( ind1, T(0) );
  SmoothIndicatorCuboid2D<T,T> interface_sharp2( ind2, T(0) );
  AnalyticalIdentity2D<T,T> phi0_sharp( zero + interface_sharp1 + interface_sharp2 );

  SuperLatticeDensity2D<T, PFDESCRIPTOR> phi( sLatticePF );
  SuperLatticePhysIncPressure2D<T, NSDESCRIPTOR> p_hydro( sLatticeNS, converter );
  SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity( sLatticeNS, converter );

  auto indicatorF = superGeometry.getMaterialIndicator(1);

  //for order parameter
  if (Cahn_const) {
    SuperRelativeErrorL1Norm2D<T> relPhiErrorNormL1(phi, phi0_diffuse, indicatorF);
    relPhiErrorNormL1(result, tmp);
    clout << "phi-L1-error(rel)=" << result[0];
    orderParameterL1RelError = result[0];

    SuperRelativeErrorL2Norm2D<T> relPhiErrorNormL2(phi, phi0_diffuse, indicatorF);
    relPhiErrorNormL2(result, tmp);
    clout << "; phi-L2-error(rel)=" << result[0];
    orderParameterL2RelError = result[0];

    SuperRelativeErrorLinfNorm2D<T> relPhiErrorNormLinf(phi, phi0_diffuse, indicatorF);
    relPhiErrorNormLinf(result, tmp);
    clout << "; phi-Linf-error(rel)=" << result[0] << std::endl;
    orderParameterLinfRelError = result[0];
  }
  else {
    SuperRelativeErrorL1Norm2D<T> relPhiErrorNormL1(phi, phi0_sharp, indicatorF);
    relPhiErrorNormL1(result, tmp);
    clout << "phi-L1-error(rel)=" << result[0];
    orderParameterL1RelError = result[0];

    SuperAbsoluteErrorL2Norm2D<T> relPhiErrorNormL2(phi, phi0_sharp, indicatorF);
    relPhiErrorNormL2(result, tmp);
    clout << "; phi-L2-error(rel)=" << result[0];
    orderParameterL2RelError = result[0];

    SuperAbsoluteErrorLinfNorm2D<T> relPhiErrorNormLinf(phi, phi0_sharp, indicatorF);
    relPhiErrorNormLinf(result, tmp);
    clout << "; phi-Linf-error(rel)=" << result[0] << std::endl;
    orderParameterLinfRelError = result[0];
  }
  //for pressure
  SuperRelativeErrorL1Norm2D<T> relPErrorNormL1(p_hydro, p_phys, indicatorF);
  relPErrorNormL1(result, tmp);
  clout << "p-L1-error(rel)=" << result[0];
  pressureL1RelError = result[0];

  SuperRelativeErrorL2Norm2D<T> relPErrorNormL2(p_hydro, p_phys, indicatorF);
  relPErrorNormL2(result, tmp);
  clout << "; p-L2-error(rel)=" << result[0];
  pressureL2RelError = result[0];

  SuperRelativeErrorLinfNorm2D<T> relPErrorNormLinf(p_hydro, p_phys, indicatorF);
  relPErrorNormLinf(result, tmp);
  clout << "; p-Linf-error(rel)=" << result[0] << std::endl;
  pressureLinfRelError = result[0];

  //for velocity
  SuperAbsoluteErrorLinfNorm2D<T> absUErrorNormLinf(velocity, zeroVelocity, indicatorF);
  absUErrorNormLinf(result, tmp);
  clout << "u-Linf-error(abs)=" << result[0] << std::endl;
  velocityLinfAbsError = result[0];
}

void getResults( SuperLattice<T,NSDESCRIPTOR>& sLatticeNS,
              SuperLattice<T,PFDESCRIPTOR>& sLatticePF,
              int iT, SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer,
              UnitConverter<T,NSDESCRIPTOR> converter,
              int Ny, int phaseLength, T w )
{
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "flatInterfaceEOC2d" );
  T deltaX = converter.getPhysDeltaX();
  SmoothIndicatorFactoredCuboid2D<T,T> interface_diffuse( {Nx/2.*deltaX, (Ny/2.+0.5)*deltaX}, 0, phaseLength*deltaX, w*deltaX/2, 0, {0,0}, 0, -1. );
  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid( sLatticeNS );
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank( sLatticeNS );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLatticeNS.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLatticePF.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    SuperLatticePhysIncPressure2D<T, NSDESCRIPTOR> p_hydro( sLatticeNS, converter );
    p_hydro.getName() = "p_hydro";
    SuperLatticeDensity2D<T, PFDESCRIPTOR> phi( sLatticePF );
    phi.getName() = "phi";
#ifdef HYBRID
    // result of signed distance function
    SuperLatticeExternalScalarField2D<T, PFDESCRIPTOR, PSI> psi( sLatticePF );
    psi.getName() = "psi";
    vtmWriter.addFunctor( psi );
#endif
    SuperLatticeExternalScalarField2D<T, NSDESCRIPTOR, RHO> rho_L( sLatticeNS );
    AnalyticalConst2D<T,T> ConversionDensity_( converter.getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR> ConversionDensity(ConversionDensity_, sLatticeNS);
    SuperIdentity2D<T,T> rho( rho_L * ConversionDensity );
    rho.getName() = "rho";
    SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity( sLatticeNS, converter );
    velocity.getName() = "u";

    AnalyticalConst2D<T,T> one ( 1. );
    AnalyticalIdentity2D<T,T> phi_ana( one + interface_diffuse );
    SuperLatticeFfromAnalyticalF2D<T, PFDESCRIPTOR> phi_analytical_diff(phi_ana, sLatticePF);
    phi_analytical_diff.getName() = "phi_analytical_diff";
    vtmWriter.addFunctor( phi_analytical_diff );

    AnalyticalConst2D<T,T> zero ( 0. );
    Vector<T,2> extend( 2./100.*charPhysLength, 0.25*charPhysLength );
    Vector<T,2> origin1( 0., 0.);
    Vector<T,2> origin2( 0., 3./4.*charPhysLength);
    IndicatorCuboid2D<T> ind1( extend, origin1 );
    IndicatorCuboid2D<T> ind2( extend, origin2 );
    SmoothIndicatorCuboid2D<T,T> interface_sharp1( ind1, T(0) );
    SmoothIndicatorCuboid2D<T,T> interface_sharp2( ind2, T(0) );
    AnalyticalIdentity2D<T,T> phi_ana_sharp( zero + interface_sharp1 + interface_sharp2 );
    SuperLatticeFfromAnalyticalF2D<T, PFDESCRIPTOR> phi_analytical_sharp(phi_ana_sharp, sLatticePF);
    phi_analytical_sharp.getName() = "phi_analytical_sharp";
    vtmWriter.addFunctor( phi_analytical_sharp );

    vtmWriter.addFunctor( p_hydro );
    vtmWriter.addFunctor( phi );
    vtmWriter.addFunctor( rho );
    vtmWriter.addFunctor( velocity );
    vtmWriter.write( iT );
  }
}

#ifdef HYBRID
T helper( SuperLattice<T,PFDESCRIPTOR>& sLatticeAC,
          SuperGeometry<T,2>& superGeometry )
{
  sLatticeAC.executePostProcessors(stage::PostPostProcess());
  int in[2];
  T out[1];
  SuperLatticeExternalScalarField2D<T, PFDESCRIPTOR, TOP> top( sLatticeAC );
  SuperLatticeExternalScalarField2D<T, PFDESCRIPTOR, BOTTOM> bottom( sLatticeAC );
  SuperAverage2D<T> top_total_(top, superGeometry, 1);
  SuperAverage2D<T> bottom_total_(bottom, superGeometry, 1);
  top_total_(out, in);
  T top_total = out[0];
  bottom_total_(out, in);
  T bottom_total = out[0];
  return top_total/bottom_total;
}
#endif

void simulate( int Ny, int phaseLength, T w, T sigma )
{
  OstreamManager clout( std::cout,"main" );
  // === 1st Step: Initialization ===
  UnitConverterFromResolutionAndRelaxationTime<T,NSDESCRIPTOR> converter(
    int   {Ny},                     // resolution
    (T)   tau_l,                    // lattice relaxation time
    (T)   charPhysLength,           // charPhysLength: reference length of simulation geometry
    (T)   Re/charPhysLength*viscosityH2O,   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   viscosityH2O,             // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   DeltaRho                  // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  clout << "Physical Surface Tension: " << converter.getConversionFactorMass()/converter.getPhysDeltaT()/converter.getPhysDeltaT()*sigma << std::endl;

  // === 2nd Step: Prepare Geometry ===
  T deltaX = converter.getPhysDeltaX();
  Vector<T,2> extend( Nx*deltaX*Ny/100., charPhysLength );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, deltaX, singleton::mpi().getSize() );
#else
  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, deltaX );
#endif
  // set periodic boundaries to the domain
  cuboidDecomposition.setPeriodicity({ true, true });
  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  loadBalancer.print();
  // Instantiation of superGeometry
  SuperGeometry<T,2> superGeometry( cuboidDecomposition,loadBalancer );
  prepareGeometry( superGeometry );

  // === 3rd Step: Prepare Lattice and Coupling ===
  SuperLattice<T,NSDESCRIPTOR> sLatticeNS( superGeometry );
  SuperLattice<T,PFDESCRIPTOR> sLatticePF( superGeometry );
  SuperLatticeCoupling coupling(
    Coupling{},
    names::NavierStokes{}, sLatticeNS,
    names::Component1{}, sLatticePF);

  prepareLattice( sLatticeNS, sLatticePF, coupling, converter, superGeometry, Ny, phaseLength, w, sigma );

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();
  std::vector<T> output = {static_cast<T>(Ny),0,0,0};
  T orderParameterL1RelError_old = 1.;

#ifdef HYBRID
  T tobo = helper(sLatticePF,superGeometry);
  coupling.template setParameter<Coupling::NONLOCALITY>(tobo);
#endif
  coupling.execute();
  for ( iT=0; iT<=maxIter; ++iT ) {
    // Collide and stream (and coupling) execution
    sLatticeNS.collideAndStream();
    sLatticePF.collideAndStream();

    sLatticePF.executePostProcessors(stage::PreCoupling());
#ifdef HYBRID
    signedDistanceFunction<stage::IterativePostProcess>(sLatticePF,superGeometry,w);
    tobo = helper(sLatticePF,superGeometry);
    coupling.template setParameter<Coupling::NONLOCALITY>(tobo);
#endif
    sLatticePF.getCommunicator(stage::PreCoupling()).communicate();
#ifdef CAHN_HILLIARD
    sLatticePF.executePostProcessors(stage::ChemPotCalc());
    sLatticePF.getCommunicator(stage::PreCoupling()).communicate();
#endif
    coupling.execute();

    // Computation and output of the results
    getResults( sLatticeNS, sLatticePF, iT, superGeometry, timer, converter, Ny, phaseLength, w );
    if ( std::isnan( sLatticeNS.getStatistics().getAverageEnergy() ) ) {
      break;
    }
    if(iT % (vtkIter) == 0) {
      error(superGeometry, sLatticePF, sLatticeNS, converter, Nx, Ny, phaseLength, w);
      clout << "Error: " << abs(orderParameterL1RelError - orderParameterL1RelError_old)/abs(orderParameterL1RelError) << std::endl;
      if(abs(orderParameterL1RelError - orderParameterL1RelError_old)/abs(orderParameterL1RelError) < 1e-5){
        clout << "Converged at: " << abs(orderParameterL1RelError - orderParameterL1RelError_old)/abs(orderParameterL1RelError) << std::endl;
        break;
      }
      orderParameterL1RelError_old = orderParameterL1RelError;
    }
  }
  timer.stop();
  timer.printSummary();
}

int main( int argc, char *argv[] )
{
  initialize( &argc, &argv );
  if (argc > 1) {
    Ny = atof(argv[1]);
  }
  std::ofstream outfile;
  outfile.open ("EOC.dat");

  //Update the geometry, lattice, and other dependent variables for each Ny
  for (size_t i=0; i<3; ++i) {
    singleton::directories().setOutputDir("./tmp/sub" + std::to_string(i) + "/");
    simulate(Ny*util::pow(2,i),
             phaseLength*util::pow(2,i),
             w*util::pow(2,i*(Cahn_const)),
             sigma/util::pow(2,i));
    if (Cahn_const) outfile << Ny*util::pow(2,i) << "," << orderParameterL2RelError <<
                        "," << pressureL2RelError << "," << velocityLinfAbsError << "\n";
    else outfile << Ny*util::pow(2,i) << "," << orderParameterL1RelError <<
             "," << pressureL2RelError << "," << velocityLinfAbsError << "\n";
  }
  outfile.close();
}

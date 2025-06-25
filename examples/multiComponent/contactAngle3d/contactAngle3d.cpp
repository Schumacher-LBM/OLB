/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Michael Rennick, Tim Bingert
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

/* contactAngle3d.cpp
 * This example shows a liquid droplet in contact with a solid wall
 * under a certain contact angle using the L. Ju et al. (2025) Well-Balanced
 * model with an energy contribution to the free energy for the order parameter
 * boundary condition. This model benefits from low spurious velocity and can
 * handle large density and viscosity ratios.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using NSDESCRIPTOR = D3Q19<RHO,NABLARHO,FORCE,EXTERNAL_FORCE,TAU_EFF,STATISTIC>;
using CHDESCRIPTOR = D3Q19<FORCE,SOURCE,SOURCE_OLD,PHIWETTING,VELOCITY,OLD_PHIU,STATISTIC,CHEM_POTENTIAL,BOUNDARY>;
using NSBulkDynamics = MPIncBGKdynamics<T,NSDESCRIPTOR>;
using CHBulkDynamics = WellBalancedCahnHilliardBGKdynamics<T,CHDESCRIPTOR>;
using Coupling = WellBalancedCahnHilliardPostProcessor;

// Parameters for the simulation domain
const int Nx = 120;                   // domain resolution x [lattice units]
const int Ny = 75;                    // domain resolution y [lattice units]
const int Nz = 120;
int diameter = 70;                    // droplet diameter [lattice units]
const int maxIter  = 500000;          // number of iterations to perform
const int vtkIter  = 1000;            // interval to save vtk output
const int statIter = 1000;            // interval to print statistics

// Characteristic physical parameters
const T L_char = 70e-6;               // charPhysLength [physical units]
const T DeltaRho = 1000.;               // physViscosity H2O liquid [physical units]
const T viscosityH2O = 9e-7;          // physViscosity H2O liquid [physical units]

// Lattice parameters for fluid properties
const T tau_mobil = 1.00;             // relaxation time for interface mobility in Cahn-Hilliard equation (mobility=(tau_mobil-0.5)/3) [lattice units]
const T tau_l = 1.0;                  // relaxation time H2O liquid lattice (kinematic viscosity=(tau_l-0.5)/3) [lattice units]
const T tau_g = 0.52;                  // relaxation time gas lattice (kinematic viscosity=(tau_g-0.5)/3) [lattice units]
const T sigma = 0.01;                 // liquid-gas surface tension [lattice units]
const T w = 4.;                       // diffuse interface width [lattice units]
const T g = 0*9.81;                   // gravitational force magnitude [lattice units]
const std::vector<T> rhos = {1., 1.}; // densities of the {gas, liquid} [lattice units]
const T theta = 100.;         // equilibrium contact angle [degrees]

// Labels for boundary and fluid locations
void prepareGeometry( SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  // Fluid nodes labelled 2
  superGeometry.rename( 0,2 );
  // Label edges as 1
  superGeometry.rename( 2, 1, {0, 1, 0} );

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

template <typename SuperLatticeCoupling>
void prepareLattice( SuperLattice<T,NSDESCRIPTOR>& sLatticeNS,
                     SuperLattice<T,CHDESCRIPTOR>& sLatticeCH,
                     SuperLatticeCoupling& coupling,
                     UnitConverter<T,NSDESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry, int diameter )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // lattice dynamics for Navier-Stokes equation, only needed in fluid domain with id 1
  sLatticeNS.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticeNS.defineDynamics<NSBulkDynamics>(superGeometry, 1);
  sLatticeNS.defineDynamics<NoDynamics>(superGeometry, 2);

  // lattice dynamics for Cahn-Hilliard equation, only needed in fluid domain with id 1
  sLatticeCH.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticeCH.defineDynamics<CHBulkDynamics>(superGeometry, 1);
  sLatticeCH.defineDynamics<NoDynamics>(superGeometry, 2);

  // set velocity to zero everywhere initially
  Vector<T,3> u(0., 0., 0.);
  AnalyticalConst3D<T,T> zeroVelocity( u );

  // analytical constants to use in initialisation
  AnalyticalConst3D<T,T> one ( 1. );
  AnalyticalConst3D<T,T> two ( 2. );
  AnalyticalConst3D<T,T> zero ( 0. );
  AnalyticalConst3D<T,T> rhov ( rhos[0] );
  AnalyticalConst3D<T,T> rhol ( rhos[1] );
  AnalyticalConst3D<T,T> tauv ( tau_g );
  AnalyticalConst3D<T,T> taul ( tau_l );

  // regions with different material ids
  auto bulk = superGeometry.getMaterialIndicator( 1 );
  auto walls = superGeometry.getMaterialIndicator( 2 );
  auto all = superGeometry.getMaterialIndicator( {0,1,2} );

  // initialise phi with a circle with centre at (Nx/2,1) and radius (diameter/2), smoothed by w/2
  IndicatorSphere3D<T> sphere( {Nx/2., 0., Nz/2.}, diameter/2. );
  SmoothIndicatorSphere3D<T,T> smoothSphere( sphere, w/2. );
  AnalyticalIdentity3D<T,T> phi( one - smoothSphere );

  // initial values for interpolated rho and tau across interfaces
  AnalyticalIdentity3D<T,T> rho( rhov + (rhol-rhov)*phi );
  AnalyticalIdentity3D<T,T> tau( tauv + (taul-tauv)*phi );

  // initial (hydrodynamic) pressure
  AnalyticalIdentity3D<T,T> pressure( zero );

  // set the initial values for the fields
  sLatticeNS.defineField<descriptors::RHO>( all, rho );                          // density
  sLatticeNS.defineField<descriptors::TAU_EFF>( superGeometry, 1, tau );         // relaxation time for navier-stokes lattice
  sLatticeCH.defineField<descriptors::SOURCE>( superGeometry, 1, zero );         // source term for Cahn-Hilliard lattice
  sLatticeCH.defineField<descriptors::PHIWETTING>( superGeometry, 1, phi );      // fluid concentration used by wetting boundaries
  sLatticeCH.defineField<descriptors::CHEM_POTENTIAL>( superGeometry, 1, zero ); // chemical potential
  sLatticeCH.defineField<descriptors::BOUNDARY>( walls, two );                   // boundary field used by wetting boundaries
  sLatticeCH.defineField<descriptors::BOUNDARY>( bulk, zero );                   // boundary field used by wetting boundaries

  // use interpolated bounce back conditions for the distribution function
  std::vector<T> origin = { -1.5, 0.5 ,-1.5};
  std::vector<T> extend = { Nx+2, Ny-2, Nz+2 };
  IndicatorCuboid3D<T> cuboid( extend, origin ); // rectangle spans from (0,0.5) to (Nx-1,Ny-1.5), in the x direction we add padding for the periodic boundaries
  setBouzidiBoundary( sLatticeNS, superGeometry, 2, cuboid );
  setBouzidiWellBalanced( sLatticeCH, superGeometry, 2, cuboid );

  // apply gravitational force
  std::vector<T> F( 3,T() );
  F[1] = -g/(converter.getPhysDeltaX()/converter.getPhysDeltaT()/converter.getPhysDeltaT());
  AnalyticalConst3D<T,T> f( F );
  sLatticeNS.defineField<descriptors::EXTERNAL_FORCE>( superGeometry, 1, f );

  // initialise distribution functions in equilibrium with initial values of phi, velocity and pressure
  sLatticeCH.defineRhoU( all, phi, zeroVelocity );
  sLatticeCH.iniEquilibrium( all, phi, zeroVelocity );
  sLatticeNS.defineRhoU( all, pressure, zeroVelocity );
  sLatticeNS.iniEquilibrium( all, pressure, zeroVelocity );

  // postprocessor to update rhowetting
  sLatticeCH.addPostProcessor<stage::PostStream>(bulk,meta::id<RhoWettingStatistics>());

  // give coupling the values of tau_l, tau_g, rho_l and rho_g so that it can update the viscosity and density
  coupling.template setParameter<Coupling::TAUS>({tau_l,tau_g});
  coupling.template setParameter<Coupling::RHOS>(rhos);

  // postprocessor to calculate the chemical potential
  sLatticeCH.addPostProcessor<stage::ChemPotCalc>(meta::id<ChemPotentialPhaseFieldProcessor>());

  // parameters needed by models
  sLatticeNS.setParameter<descriptors::OMEGA>( 1./tau_l );     // collision rate for Navier-Stokes lattice
  sLatticeCH.setParameter<descriptors::OMEGA>( 1./tau_mobil ); // collision rate for Cahn-Hilliard lattice
  sLatticeCH.setParameter<descriptors::THETA>( (M_PI-theta*M_PI/180.) );        // contact angle
  sLatticeCH.setParameter<descriptors::INTERFACE_WIDTH>( w );  // diffuse interface width
  sLatticeCH.setParameter<descriptors::SCALAR>( sigma );       // surface tension

  // define fields that must be communicated by mpi (needed for finite difference gradients)
  {
    auto& communicator = sLatticeCH.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(2);
    communicator.requestField<STATISTIC>();
    communicator.requestField<PHIWETTING>();
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.exchangeRequests();
  }

  // some initial steps
  sLatticeCH.executePostProcessors(stage::PreCoupling());
  sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();
  sLatticeCH.executePostProcessors(stage::ChemPotCalc());
  sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();
  sLatticeNS.initialize();
  sLatticeCH.initialize();
  sLatticeCH.iniEquilibrium( all, phi, zeroVelocity );
  sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// for statistic output and vtk saving
T getResults( SuperLattice<T,NSDESCRIPTOR>& sLatticeNS,
              SuperLattice<T,CHDESCRIPTOR>& sLatticeCH,
              int iT, SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
              UnitConverter<T,NSDESCRIPTOR> converter )
{
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "contactAngle3d" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid( sLatticeNS );
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank( sLatticeNS );
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
    sLatticeCH.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  T contact_angle = 0;
  // Writes the VTK files
  if ( iT%vtkIter==0 ) {

    SuperLatticeDensity3D<T, NSDESCRIPTOR> p_hydro( sLatticeNS );
    p_hydro.getName() = "p_hydro";

    SuperLatticeField3D<T, CHDESCRIPTOR, STATISTIC> phi( sLatticeCH );
    phi.getName() = "phi";

    SuperLatticeExternalScalarField3D<T, NSDESCRIPTOR, RHO> rho_L( sLatticeNS );
    AnalyticalConst3D<T,T> ConversionDensity_( converter.getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF3D<T, NSDESCRIPTOR> ConversionDensity(ConversionDensity_, sLatticeNS);
    SuperIdentity3D<T,T> rho ( rho_L * ConversionDensity );
    rho.getName() = "rho";

    SuperLatticeVelocity3D<T, NSDESCRIPTOR> velocity( sLatticeNS );
    velocity.getName() = "u";

    vtmWriter.addFunctor( p_hydro );
    vtmWriter.addFunctor( phi );
    vtmWriter.addFunctor( rho );
    vtmWriter.addFunctor( velocity );
    vtmWriter.write( iT );

    // function to interpolate phi at a given position
    AnalyticalFfromSuperF3D<T,T> interpolPhi( phi, true, 1 );
    T point1 = 0, point2 = 0, point3 = 0;

    // contact angle fitting, will be done by finding three points on the circlular droplet

    // point 1, x = ? y = 2
    T pos[3] = {0., 2., Nz/2.};
    for (int ix=0; ix<0.5*Nx; ix++) {
      T phi1, phi2;
      pos[0] = ix;
      interpolPhi( &phi1, pos );
      // if phi at ix is greater than 0.5, we have passed the interface of the liquid phase
      if (phi1 < 0.5) {
        pos[0] = (ix-1); // check the value at the previous point
        interpolPhi( &phi2, pos );
        point1 = ix - 1 + (0.5-phi2)/(phi1-phi2); // interpolate between these values
        break;
      }
    }

    // point 2, x = ? y = 2
    for (int ix=0.5*Nx; ix<Nx; ix++) {
      T phi1, phi2;
      pos[0] = ix;
      interpolPhi( &phi1, pos );
      if (phi1 > 0.5) {
        pos[0] = (ix-1);
        interpolPhi( &phi2, pos );
        point2 = ix - 1 + (0.5-phi2)/(phi1-phi2);
        break;
      }
    }

    // point 3, x = Nx/2 y = ?
    pos[0] = 0.5*(Nx);
    for (int iy=3; iy<Ny; iy++) {
      T phi1, phi2;
      pos[1] = iy;
      interpolPhi( &phi1, pos );
      if (phi1 > 0.5) {
        pos[0] = (iy-1);
        interpolPhi( &phi2, pos );
        point3 = iy - 1 + (0.5-phi2)/(phi1-phi2);
        break;
      }
    }

    T x1=point1;
    T y1=2;
    T x2=point2;
    T y2=2;
    T x3=0.5*(Nx);
    T y3=point3;

    // estimate centre and radius from three points, we must solve three simulatneous equations
    T s1 = x1*x1 + y1*y1;
    T s2 = x2*x2 + y1*y1;
    T s3 = x3*x3 + y3*y3;
    T M11 = x1*y2 + x2*y3 + x3*y1 - (x2*y1 + x3*y2 + x1*y3);
    T M12 = s1*y2 + s2*y3 + s3*y1 - (s2*y1 + s3*y2 + s1*y3);
    T M13 = s1*x2 + s2*x3 + s3*x1 - (s2*x1 + s3*x2 + s1*x3);
    T xc =  0.5*M12/M11;
    T yc = -0.5*M13/M11;
    T r = sqrt(pow(x2 - xc,2) + pow(y2 - yc,2));

    contact_angle = util::acos(-(yc-0.5) / r) * 180 / M_PI;

    clout << "Numerical Contact angle: " << contact_angle << std::endl;
    clout << "Analytical contact angle: " << theta <<  std::endl;
  }
  return contact_angle;
}

// run the simulation
void simulate( int diameter )
{
  OstreamManager clout( std::cout,"main" );
  // === 1st Step: Initialization ===
  UnitConverterFromResolutionAndRelaxationTime<T,NSDESCRIPTOR> converter(
    int   {diameter},               // resolution
    (T)   tau_l,                    // lattice relaxation time
    (T)   L_char,                   // charPhysLength: reference length of simulation geometry
    (T)   0,                        // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   viscosityH2O,             // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   DeltaRho                  // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  // === 2nd Step: Prepare Geometry ===
  // Instantiation of a cuboidGeometry with weights
  std::vector<T> extend = { Nx, Ny, Nz };
  std::vector<T> origin = { 0., 0., 0. };
  IndicatorCuboid3D<T> cuboid(extend,origin);
#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, 1., singleton::mpi().getSize() );
#else
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, 1. );
#endif

  // Set periodic boundaries to the domain
  cuboidDecomposition.setPeriodicity({ true, false, true });

  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  loadBalancer.print();
  // Instantiation of superGeometry
  SuperGeometry<T,3> superGeometry( cuboidDecomposition,loadBalancer );
  prepareGeometry( superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,NSDESCRIPTOR> sLatticeNS( superGeometry );
  SuperLattice<T,CHDESCRIPTOR> sLatticeCH( superGeometry );
  SuperLatticeCoupling coupling(
      WellBalancedCahnHilliardPostProcessor{},
      names::NavierStokes{}, sLatticeNS,
      names::Component1{}, sLatticeCH);
  coupling.restrictTo(superGeometry.getMaterialIndicator({1}));
  prepareLattice( sLatticeNS, sLatticeCH, coupling, converter, superGeometry, diameter );

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();
  std::ofstream outfile;
  outfile.open ("contactAngleVsTime.dat");
  T old_angle = 1.;

  for ( iT=0; iT<=maxIter; ++iT ) {
    // Collide and stream (and coupling) execution
    sLatticeNS.collideAndStream();
    sLatticeCH.collideAndStream();

    sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();
    sLatticeCH.executePostProcessors(stage::PreCoupling());
    sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();
    sLatticeCH.executePostProcessors(stage::ChemPotCalc());
    sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();

    coupling.execute();

    // Computation and output of the results
    T angle = getResults( sLatticeNS, sLatticeCH, iT, superGeometry, timer, converter );
    if ( iT%vtkIter == 0 ) {
      outfile << iT*converter.getPhysDeltaT() << "," ;
      outfile << angle << "\n";
      if ( fabs(angle-old_angle)/fabs(old_angle) < 5e-4 ) {
        clout << "contact angle converged..." << std::endl;
        break;
      }
      old_angle = angle;
    }
    if ( std::isnan( sLatticeNS.getStatistics().getAverageEnergy() ) ) {
      break;
    }

  }
  outfile.close();
  timer.stop();
  timer.printSummary();
}

int main( int argc, char *argv[] )
{
  initialize( &argc, &argv );
  if (argc > 1) {
    diameter = atof(argv[1]);
  }

  singleton::directories().setOutputDir( "./tmp/" );
  simulate( diameter );
}
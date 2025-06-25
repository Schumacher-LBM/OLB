/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Johannes Nguyen, Stephan Simonis, Sam Avis
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

 /* binaryShearFlow2d.cpp
  * A circular domain of fluid phase II is immersed in a square
  * filled with fluid phase I. The top and bottom walls are
  * moving in opposite directions, such that the droplet consisting
  * of phase II is exposed to shear flow and deforms accordingly.
  * The default parameter setting is given in
  *  A.E. Komrakova et al. International Journal of Multiphase
  *  Flow 59 (2014) 24–43.
  * The present example uses the ternary free energy model by
  *  C Semprebon et al. Physical Review E 93.3 (2016) p. 033305.
  * Reference results are published in
  *  S. Simonis et al. Discrete and Continuous Dynamical
  *  Systems - Series S (2023), doi:10.3934/dcdss.2023069
  *
  * The droplet breakup case should be run on parallel mode
  * because of increased computation time.
  */

#include <olb.h>

#include <chrono>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
#define DESCRIPTOR D2Q9<CHEM_POTENTIAL,FORCE>

// // DEFAULT [lu] Setup Komrakova et al. 2013 "Lattice Boltzmann ... " Fig. 9
// // Simulated Domain LxHxW = 12*radius x 8*radius x 2*radius;
// // Resolution: radius = 30 [lu]
// // Re = 0.0625
// // Ca = 0.42 (breakup for Stokes flow)
// // Ch = 0.0379
// // Pe = 0.43
// // Breakup happens at normTime ~30

// Guideline parameters:
// Ca = 0.25 for steady state
// Ca = 1.0  for drop breakup
// Recommended xi = 1.14, min 1.0

// Set physical input parameters
const T radius = 1.;                          // [pu] physical droplet radius
const T nx = 20 * radius;                     // [pu] domain length H = 12*radius
const T ny = 4 * radius;                      // [pu] domain height H = 4*radius (this makes no difference to komra 8, since Re equal!)
const T contVisc = 1.;                        // [pu] physical(!) kinematic viscosity of continuous phase
const T physDensity = 1.;                     // [pu] physical density

// Set numerical parameters
const T N = 10;                               // resolution N cells per charL = radius/2;
const T latticeRelaxationTime = 1.;           // [lu] tau for first lattice (f)
const T latticeRelaxationTime_g = 1.;         // [lu] caution: not used..., i.e. tau_f = tau_g

// Set (model) dimensionless numbers
const T dropRe = 1.;                          // definition: shearRate * radius^2 / continuousViscosity
const T capillaryNr = 1.;                     // capillary number for droplet breakup
//const T capillaryNr = 0.25;                 // definition: radius * shearRate * viscosity(non kinematic) / surfTen


// Set (numerical) dimensionless numbers
const T cahnNr = 0.057;                       // Cahn Number = xi / radius;
const T pecletNr = .43;                       // Peclet Number = shearRate * radius * xi / (M A)

// Set time interval to simulate and normalized output
const T maxNormTime = 30.;                    // normalized time normT = T * shearRate [either in lu or pu]
const T vtkNormIter = 0.2;                    // normalized vtk output (i.e. per normalized time step, one output!)
const T statNormIter = 0.2;                   // normalized stat output (i.e. per normalized time step, one output!)
const T physInterval = 0.1;                   // interval for the convergence check in s
const T residuum = 1e-5;                      // residuum for the convergence check
const T defResiduum = 0.0015;                 // deformation residuum criterion

// Set contact angles
// no boundary touching here...
// h_i= Parameter related to resulting contact angle of the boundary. [lattice units]
const T h1 = 0.;                              // Contact angle 90 degrees   [lattice units]
const T h2 = 0.;                              // Contact angle 90 degrees   [lattice units]

// Initialize remaining parameters (conversion, computation, free energy)
T shearRate;        // [pu] physical shear rate [1/s]
T shearVelo;        // [pu] physical shear speed corr. to dropRe [m/s]

T surfTen;          // [pu] physical surface tension
T surfTenLatt;      // [lu] lattice surface tension

T kappa;            // [lu] this is kappa1+kappa2, to simplify, and meet komrakova model. this is NOT kappa in Komrakova2013!!
T kappa1;           // [lu] variable for interface tension tuning
T kappa2;           // [lu] variable for interface tension tuning
T alpha;            // [lu] free energy alpha in surface tension
T gama;             // [lu] gama "Diffusivity of the interface" same as in Komrakova2013

T xiThickness;      // [lu] characteristical thickness of interface (xi from Komrakova2013)
T xiThicknessPU;    // [pu] interface thickness in physical units
T kappaKomra;       // [lu] this is kappa from Komrakova2013
T aKomra;           // [lu] this is A from Komrakova2013
T transportCoeff;   // [lu] this is M from Komrakova2013

T maxPhysTime;      // [pu] maximum physical time
int maxIter;        // [lu] maximum lattice time

T physVtkIter;      // [pu] physical time step, when output vtk happens
T physStatIter;     // [pu] physical time step, when output happens
int vtkIter;        // [lu] lattice time step, when output happens
int statIter;       // [lu] lattice timestep when statistics apper in terminal

const T PI_2 = 1.57079632679;   // PI/2

// Interpolation for large axis - smallest circle outside droplet
// This will only work for single ellipsoid droplets with stationary center
// Particularly in cases of droplet breakup, this will not work
std::vector<T> circleApproximationBig(SuperLattice<T, DESCRIPTOR>& sLattice2,
      T minDistance,    // stop criterion: minimum distance for radii
      int points,       // number of points on the circle
      T factor ){       // expansion and retraction of the circle with factor 1+factor or 1-factor
  SuperLattice<T, DESCRIPTOR> *sLatticeTest = &sLattice2;;
  SuperLatticeDensity2D<T, DESCRIPTOR> densityCircle( *sLatticeTest );

  AnalyticalFfromSuperF2D<T> aDensity(densityCircle, true, 1);

  T oldRadius = 0.;
  T newRadius = radius;

  T safety = 0.5;
  T coordinates[2] = {T(0),T(0)};
  T theta = 0.;

  while(std::abs(oldRadius-newRadius) > minDistance){
    T maxAngle = PI_2;      // third quadrant, as shear direction is known beforehand
    if( ny/2. + newRadius >= ny - safety) maxAngle = util::asin( ( ny/2. - safety) / newRadius );       // eliminate points in proximity to no-slip wall (phi=0 there)
    for(int i=0; i<=points-1; i++){
      coordinates[0] = nx/2. + util::cos(maxAngle * (T) i/ ((T) (points-1.)) ) * newRadius;
      coordinates[1] = ny/2. + util::sin(maxAngle * (T) i/ ((T) (points-1.)) ) * newRadius;
      T output[1];

      aDensity( output, coordinates );

      if(output[0] <= -0.0) {
        oldRadius = newRadius;
        newRadius *= 1. + factor;
        //theta is only assigned if the circle intersects the droplet
        theta = 90./PI_2 * util::atan( (coordinates[1] - ny/2.) / (coordinates[0] - nx/2.) );
        break;
      } else if ( i == points - 1){
        oldRadius = newRadius;
        newRadius *= 1. - factor;
      }
    }
    factor *= 0.8;
  }
  return {oldRadius, theta};
}

// Interpolation for small axis - largest circle inside droplet
T circleApproximationSmall(SuperLattice<T, DESCRIPTOR>& sLattice2,
      T minDistance,
      int points,
      T factor){
  SuperLattice<T, DESCRIPTOR> *sLatticeTest = &sLattice2;;
  SuperLatticeDensity2D<T, DESCRIPTOR> densityCircle( *sLatticeTest );

  AnalyticalFfromSuperF2D<T> aDensity(densityCircle, true, 1);

  T oldRadius = 0.;
  T newRadius = radius;

  T coordinates[2] = {T(0),T(0)};

    while(std::abs(oldRadius-newRadius) > minDistance){
    T maxAngle = 2 * PI_2;
    for(int i=0; i<=points-1; i++){
      coordinates[0] = nx/2. + util::cos(maxAngle * (T) i/ ((T) (points-1.)) ) * newRadius;
      coordinates[1] = ny/2. + util::sin(maxAngle * (T) i/ ((T) (points-1.)) ) * newRadius;
      T output[1];

      aDensity( output, coordinates );

      if(output[0] >= -0.0) {
        oldRadius = newRadius;
        newRadius *= 1. - factor;
        break;
      } else if ( i == points - 1){
        oldRadius = newRadius;
        newRadius *= 1. + factor;
      }
    }
    factor *= 0.8;
  }
  return oldRadius;
}

// The accuracy of this is dependent on the resolution
// Additionally droplet deformation is only valid for ellipsoidal droplets
T measureDroplet(SuperLattice<T, DESCRIPTOR>& sLattice2,
                    SuperGeometry<T,2>& superGeometry,
                    UnitConverter<T,DESCRIPTOR>& converter){
  OstreamManager clout( std::cout,"Measurement" );

  //axis approximation
  std::vector<T> circleMeasurementsBig = circleApproximationBig( sLattice2, 0.001, 5000, 0.8 );
  T circleMeasurementsSmall = circleApproximationSmall( sLattice2, 0.001, 5000, 0.8 );

  T deformation = (circleMeasurementsBig[0] - circleMeasurementsSmall) / (circleMeasurementsBig[0] + circleMeasurementsSmall);

  clout << "Circle Measurement:  " << " L: " << circleMeasurementsBig[0] << ", H: " << circleMeasurementsSmall << ", theta: " << circleMeasurementsBig[1] << ", deformation: " << deformation<< std::endl;
  return deformation;
}

void prepareGeometry( SuperGeometry<T,2>& superGeometry,
                      UnitConverter<T,DESCRIPTOR> const& converter )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0, 2 );

  // bulk, MN=1
  superGeometry.rename( 2, 1, {0, 1} );

  T eps = converter.getPhysLength(1);
  T edge = 0.5*eps;
  T safety = 2*eps;

  // top wall, MN=3
  std::vector<T> origin = {T(0) - safety, T(ny)-edge};
  std::vector<T> extend = {nx + 2*safety, edge+safety};
  IndicatorCuboid2D<T> top( extend, origin );
  superGeometry.rename(2, 3, top );

  // bottom wall, MN=4
  origin[1] = 0.0 - safety;
  IndicatorCuboid2D<T> bottom( extend, origin);
  superGeometry.rename(2, 4, bottom);

  // clean up
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( SuperLattice<T, DESCRIPTOR>& sLattice1,
                     SuperLattice<T, DESCRIPTOR>& sLattice2,
                     UnitConverter<T, DESCRIPTOR>& converter,
                     SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // define lattice Dynamics
  sLattice1.defineDynamics<ForcedBGKdynamics>( superGeometry, 1 );
  sLattice2.defineDynamics<FreeEnergyBGKdynamics>( superGeometry, 1 );

  // moving walls (inlet boundary with tangetial velocity condition)
  T omega = converter.getLatticeRelaxationFrequency();

  auto topIndicator = superGeometry.getMaterialIndicator(3);
  boundary::set<boundary::FreeEnergyVelocity>(sLattice1, topIndicator);
  boundary::set<boundary::FreeEnergyOrderParameter>(sLattice2, topIndicator);

  auto bottomIndicator = superGeometry.getMaterialIndicator(4);
  boundary::set<boundary::FreeEnergyVelocity>(sLattice1, bottomIndicator);
  boundary::set<boundary::FreeEnergyOrderParameter>(sLattice2, bottomIndicator);

  sLattice1.setParameter<OMEGA>(omega);
  sLattice2.setParameter<OMEGA>(omega);
  sLattice2.setParameter<collision::FreeEnergy::GAMMA>(gama);

  // bulk initial conditions
  // define circular domain for fluid 2
  std::vector<T> v( 2,T() );
  AnalyticalConst2D<T,T> zeroVelocity( v );

  AnalyticalConst2D<T,T> one ( 1. );

  SmoothIndicatorCircle2D<T,T> circle( {T(nx)/T(2), T(ny)/T(2)}, radius, converter.getPhysLength(alpha) );

  AnalyticalIdentity2D<T,T> rho( one );
  AnalyticalIdentity2D<T,T> phi( one - circle - circle );

  // shear velocity
  AnalyticalConst2D<T,T> uTop   ( converter.getCharLatticeVelocity(), T( 0 ) );
  AnalyticalConst2D<T,T> uBottom( -1.*converter.getCharLatticeVelocity(), T( 0 ) );

  // top
  sLattice1.defineRhoU( topIndicator, rho, uTop );
  sLattice2.defineRho( topIndicator, phi );

  // bottom
  sLattice1.defineRhoU( bottomIndicator, rho, uBottom );
  sLattice2.defineRho( bottomIndicator, phi );

  // equilibrium population initialization
  sLattice1.iniEquilibrium( superGeometry, 1, rho, zeroVelocity );
  sLattice2.iniEquilibrium( superGeometry, 1, phi, zeroVelocity );
  sLattice1.iniEquilibrium( superGeometry, 3, rho, uTop );
  sLattice2.iniEquilibrium( superGeometry, 3, phi, uTop );
  sLattice1.iniEquilibrium( superGeometry, 4, rho, uBottom );
  sLattice2.iniEquilibrium( superGeometry, 4, phi, uBottom );

  sLattice1.initialize();
  sLattice2.initialize();

  sLattice1.communicate();
  sLattice2.communicate();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults( SuperLattice<T, DESCRIPTOR>& sLattice2,
                 SuperLattice<T, DESCRIPTOR>& sLattice1, int iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer,
                 UnitConverter<T, DESCRIPTOR> converter) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "binaryShearFlow2d" );

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
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    sLattice1.setProcessingContext(ProcessingContext::Evaluation);
    AnalyticalConst2D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> half(half_, sLattice1);

    SuperLatticeDensity2D<T, DESCRIPTOR> density1( sLattice1 );
    density1.getName() = "rho";
    SuperLatticeDensity2D<T, DESCRIPTOR> density2( sLattice2 );
    density2.getName() = "phi";

    SuperIdentity2D<T,T> c1 (half*(density1+density2));
    c1.getName() = "density-fluid-1";
    SuperIdentity2D<T,T> c2 (half*(density1-density2));
    c2.getName() = "density-fluid-2";

    // velocity (seperately and combined)
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity1( sLattice1, converter );
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity2( sLattice2, converter );
    SuperIdentity2D<T,T> velocityMittel (half*(velocity1+velocity2));
    velocityMittel.getName() = "velocity-Mittel";
    velocity1.getName() = "velocity1";
    velocity2.getName() = "velocity2";

    vtmWriter.addFunctor( density1 );
    vtmWriter.addFunctor( density2 );
    vtmWriter.addFunctor( c1 );
    vtmWriter.addFunctor( c2 );
    vtmWriter.addFunctor( velocityMittel );
    vtmWriter.write( iT );
  }
}

int main( int argc, char *argv[] )
{
  // === 1st Step: Initialization ===

  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  // use remaining equations from Komrakova2013 to calculate missing parameters:
  // 0. compute velocity from Re (first shearrate, than velocity)
  // (do that before converter, since it is required as argument!)
  shearRate = contVisc * dropRe / util::pow(radius, 2);
  shearVelo = ny * shearRate / 2.;

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    (T)   N,                            // resolution
    (T)   latticeRelaxationTime,        // lattice relaxation time
    (T)   2*util::pow(radius, 2.)/ny,   // charPhysLength: reference length of simulation geometry (dropRe!!)
    (T)   shearVelo,                    // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   contVisc,                     // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physDensity                   // physDensity: physical density in __kg / m^3__
  );

  // 1. compute surface tension from Ca
  surfTen = radius * shearRate * physDensity * contVisc / capillaryNr;
  surfTenLatt = surfTen / (converter.getConversionFactorPressure() * converter.getPhysDeltaX());

  // 2. compute interface thickness from Ch
  xiThickness = cahnNr * converter.getLatticeLength(radius);
  xiThicknessPU = converter.getPhysLength(xiThickness);

  // 3. compute kappa from xi (Komrakova2013 equation (12))
  kappaKomra = (3./4.) * surfTenLatt * xiThickness;

  // 4. compute A from kappa (Komrakova2013 equation (11), neglect minus sign, since solely multiplications here)
  aKomra = 2. * kappaKomra / util::pow(xiThickness, 2);

  // 5. compute M from Pe (Komrakova2013 equation (17), there's a typo in the paper, should be a minus on the right hand side... we dropped that anyway.)
  transportCoeff = shearRate * converter.getConversionFactorTime() * converter.getLatticeLength(radius) * xiThickness / (pecletNr * aKomra);

  // 6. compute gamma from M (Komrakova2013 equation (8), caution: \Delta t is unity in lattice units!)
  gama = transportCoeff / (converter.getLatticeRelaxationTime() - .5);

  // Parameter fit from Komrakova to OpenLB (Semprebon2016):
  kappa1 = 4. * aKomra;
  kappa2 = kappa1;
  alpha = xiThickness / 2.;

  // transform normalized time to physical and then to lattice:
  maxPhysTime = maxNormTime / shearRate;
  maxIter = converter.getLatticeTime( maxPhysTime );

  physVtkIter = vtkNormIter / shearRate;
  physStatIter = statNormIter / shearRate;
  vtkIter = converter.getLatticeTime( physVtkIter );
  statIter = converter.getLatticeTime( physStatIter );;


  // Prints the converter log as console output
  converter.print();

  // now print everything additionally to the unit converter (helps for double-checking)
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "shear rate             " << shearRate      << "[pu], " << converter.getLatticeVelocity(shearRate) << "[lu]" << std::endl;
  clout << "shear rate             " << shearRate      << "[pu], " << shearRate * converter.getPhysDeltaT()   << "[lu]" << std::endl;
  clout << "Re_drop                " << dropRe         << std::endl;
  clout << "Ca                     " << capillaryNr    << std::endl;
  clout << "radius                 " << radius         << "[pu], " << converter.getLatticeLength(radius)      << "[lu]" << std::endl;
  clout << "normalized t_end       " << maxNormTime    << std::endl;
  clout << "surface tension        " << surfTen        << "[pu], " << surfTenLatt << "[lu]" << std::endl;
  clout << "Cahn Nr.               " << cahnNr         << std::endl;
  clout << "xi                     " << xiThickness    << "[lu]" << std::endl;
  clout << "alpha                  " << alpha          << "[lu]" << std::endl;
  clout << "kappa1                 " << kappa1         << "[lu]" << std::endl;
  clout << "kappa2                 " << kappa2         << "[lu]" << std::endl;
  clout << "Peclet Nr.             " << pecletNr       << std::endl;
  clout << "A                      " << aKomra         << "[lu]" << std::endl;
  clout << "Mobility M             " << transportCoeff << "[lu]" << std::endl;
  clout << "Gamma                  " << gama           << "[lu]" << std::endl;
  clout << "vtkIter                " << vtkIter        << "[lu]" << std::endl;
  clout << "statIter               " << statIter       << "[lu]" << std::endl;
  clout << "maxIter                " << maxIter        << "[lu]" << std::endl;
  clout << "physInterval convCheck " << physInterval   << "[pu]" << std::endl;
  clout << "residuum convCheck     " << residuum       << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;

  // === 2nd Step: Prepare Geometry ===
  std::vector<T> extend = { nx, ny };
  std::vector<T> origin = { 0, 0 };
  IndicatorCuboid2D<T> cuboid(extend,origin);

  //changed because of segmentation fault
  #ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 1 * singleton::mpi().getSize();
  #else
  const int noOfCuboids = 1;
  #endif
  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, converter.getPhysDeltaX(), noOfCuboids);

  // set periodic boundaries to the domain
  cuboidDecomposition.setPeriodicity({ true, false });

  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  loadBalancer.print();

  // Instantiation of superGeometry
  SuperGeometry<T,2> superGeometry( cuboidDecomposition,loadBalancer );

  prepareGeometry( superGeometry, converter );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice1( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLattice2( superGeometry );

  prepareLattice( sLattice1, sLattice2, converter, superGeometry );

  // Add the lattice couplings
  // The chemical potential coupling must come before the force coupling
  SuperLatticeCoupling coupling1(
  ChemicalPotentialCoupling2D{},
  names::A{}, sLattice1,
  names::B{}, sLattice2);

  coupling1.template setParameter<ChemicalPotentialCoupling2D::ALPHA>(alpha);
  coupling1.template setParameter<ChemicalPotentialCoupling2D::KAPPA1>(kappa1);
  coupling1.template setParameter<ChemicalPotentialCoupling2D::KAPPA2>(kappa2);

  SuperLatticeCoupling coupling2(
  ForceCoupling2D{},
  names::A{}, sLattice2,
  names::B{}, sLattice1);

  coupling1.restrictTo(superGeometry.getMaterialIndicator({1}));
  coupling2.restrictTo(superGeometry.getMaterialIndicator({1}));

  // Walls
  SuperLatticeCoupling coupling3(
  InletOutletCoupling2D{},
  names::A{}, sLattice2,
  names::B{}, sLattice1);

  coupling3.restrictTo(superGeometry.getMaterialIndicator({3,4}));

  sLattice1.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());
  sLattice2.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

  {
    auto& communicator = sLattice1.getCommunicator(stage::PostCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestOverlap(sLattice1.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = sLattice2.getCommunicator(stage::PreCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestField<RhoStatistics>();
    communicator.requestOverlap(sLattice2.getOverlap());
    communicator.exchangeRequests();
  }

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  util::Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  T old_deformation = 0;
  T deformation = 0;

  for ( iT=0; iT<=maxIter; ++iT ) {

    sLattice2.setProcessingContext(ProcessingContext::Evaluation);
    sLattice1.setProcessingContext(ProcessingContext::Evaluation);
    // Measure droplet
    if ( iT%statIter==0 ) {
      old_deformation = deformation;
      clout<<"--------------- t: "<< converter.getPhysTime(iT) <<" ----------------"<<std::endl;
      deformation = measureDroplet( sLattice2, superGeometry, converter );
    }

    // Convergence check
    if ( converge.hasConverged() && util::abs(deformation - old_deformation) < defResiduum ) {
      clout << "Simulation converged." << std::endl;
      getResults( sLattice2, sLattice1, iT, superGeometry, timer, converter );
      break;
    }
    sLattice1.setProcessingContext(ProcessingContext::Simulation);
    sLattice2.setProcessingContext(ProcessingContext::Simulation);

    // Computation and output of the results
    getResults( sLattice2, sLattice1, iT, superGeometry, timer, converter );

    // Collide and stream execution
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();

    // MPI communication for lattice data
    sLattice1.communicate();
    sLattice2.communicate();

    sLattice1.executePostProcessors(stage::PreCoupling());
    sLattice2.executePostProcessors(stage::PreCoupling());

    // Execute coupling between the two lattices
    coupling1.execute();
    coupling2.execute();
    coupling3.execute();

  converge.takeValue( sLattice2.getStatistics().getAverageRho(), true );
  }

  timer.stop();
  timer.printSummary();

}

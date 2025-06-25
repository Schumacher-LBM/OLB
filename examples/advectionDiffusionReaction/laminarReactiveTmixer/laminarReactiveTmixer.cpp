/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Fedor Bukreev, Adrian Kummerländer, Mathias J. Krause
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

 /* laminarReactiveTmixer.cpp:
 * Benchmark Simulation of Laminar Reactive Micromixing Using Lattice Boltzmann Methods
 * Fedor Bukreev, Adrian Kummerländer, Julius Jeßberger, Dennis Teutscher, Stephan Simonis,
 * Dieter Bothe, and Mathias J. Krause, AIAA Journal, https://doi.org/10.2514/1.J064234
 *
 * This example is a new benchmark for laminar reactive micromixing in a T-shaped mixer
 * with A+B->C reaction and special Smagorinsky-like stabilization approach.
 */


#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;
using NSDESCRIPTOR = D3Q19<>;
using NSBulkDynamics = BGKdynamics<T,NSDESCRIPTOR>;

#define plotOverLine

using RADDESCRIPTOR = D3Q7<VELOCITY,SOURCE,OMEGA>;
using ReactionADBulkDynamics = SourcedLimitedAdvectionDiffusionBGKdynamics<T,RADDESCRIPTOR>;

// Parameters for the simulation setup
const T lx         = 1600e-6;                         // length of channel
const T ly         = 800e-6;                          // width of channel
const T lz         = 100e-6;                          // hight of channel
const T lx_inlet   = 100e-6;                          // inlet length
const T ly_outlet  = 200e-6;                          // outlet width
const int N        = 40;                              // resolution of the model
const T nu         = 1.002e-6;                        // water viscosity
const T density    = 1000.;                           // water density
const T hydrInlet  = 2.*ly_outlet*lz/(ly_outlet + lz);// hydraulic inlet diameter
const T u_mean     = 1.4;                             // mean velocity at inlets
const T u_max      = u_mean*2.1;                      // max velocity for poiseuille profile
const T residenceT = (lx + ly/2)/u_mean;              // residence time
const T maxPhysT   = 25*residenceT;                   // max. NSE simulation time in s, SI unit
const T maxPhysTAD = 0.01;                            // max. ADE simulation time in s, SI unit
const T D_A        = 1.6e-9;                          // main diffusion coefficient of component A
const T D_B        = 2.e-10;                          // diffusion coefficient of component B
const T D_P        = 2.e-10;                          // diffusion coefficient of reaction product


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,NSDESCRIPTOR> const& converter,
                      SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  T dx = converter.getConversionFactorLength();

  superGeometry.rename( 0,2 );

  Vector<T,3> extendReactor( lx, ly_outlet, lz );
  Vector<T,3> originReactor( 0., -ly_outlet/2, -lz/2 );
  IndicatorCuboid3D<T> reactor( extendReactor, originReactor );

  Vector<T,3> extendInflow( lx_inlet, ly, lz );
  Vector<T,3> originInflow( 0., -ly/2, -lz/2 );
  IndicatorCuboid3D<T> inflow( extendInflow, originInflow );

  superGeometry.rename( 2, 1, reactor );
  superGeometry.rename( 2, 1, inflow );

  Vector<T,3> extendInlet1( lx_inlet, dx, lz);
  Vector<T,3> originInlet1( 0., -ly/2, -lz/2 );
  IndicatorCuboid3D<T> inlet1( extendInlet1, originInlet1 );
  superGeometry.rename( 1, 3, inlet1 );

  Vector<T,3> extendInlet2( lx_inlet, dx, lz);
  Vector<T,3> originInlet2( 0., ly/2-dx, -lz/2 );
  IndicatorCuboid3D<T> inlet2( extendInlet2, originInlet2 );
  superGeometry.rename( 1, 4, inlet2 );

  Vector<T,3> extendOutlet( dx, ly_outlet, lz );
  Vector<T,3> originOutlet( lx-dx, -ly_outlet/2, -lz/2 );
  IndicatorCuboid3D<T> outlet( extendOutlet, originOutlet );
  superGeometry.rename( 1, 5, outlet );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLatticeNS( UnitConverter<T,NSDESCRIPTOR> const& converter,
                     SuperLattice<T,NSDESCRIPTOR>& sLatticeNS,
                     SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLatticeNS" );
  clout << "Prepare Lattice for NSE..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLatticeNS.defineDynamics<NoDynamics>(superGeometry, 0);

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inlet1 (-Y))
  // Material=4 -->bulk dynamics (inlet2 (+Y))
  // Material=5 -->bulk dynamics (outlet)
  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 3, 4, 5});
  sLatticeNS.defineDynamics<NSBulkDynamics>(bulkIndicator);

  // Material=2 -->bounce back
  sLatticeNS.defineDynamics<BounceBack>(superGeometry, 2);

  //if interpolated boundary conditions are chosen
  boundary::set<boundary::InterpolatedVelocity>(sLatticeNS, superGeometry, 3);
  boundary::set<boundary::InterpolatedVelocity>(sLatticeNS, superGeometry, 4);
  boundary::set<boundary::InterpolatedPressure>(sLatticeNS, superGeometry, 5);

  // Initial conditions
  AnalyticalConst3D<T,T> rho1( 1. );
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> u( velocity );

  //Initialize all values of distribution functions to their local equilibrium
  sLatticeNS.defineRhoU( bulkIndicator, rho1, u );
  sLatticeNS.iniEquilibrium( bulkIndicator, rho1, u );

  sLatticeNS.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLatticeNS.initialize();

  clout << "Prepare Lattice for NSE... OK" << std::endl;
}

void prepareLatticeCRAD(SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
                        SuperLattice<T, RADDESCRIPTOR>*& CRADlattice,
                        SuperGeometry<T,3>& superGeometry,
                        UnitConverter<T,NSDESCRIPTOR> const& converter,
                        T rhoA_ = 1.e-8, T rhoB_ = 1.e-8, T D = D_P)
{
  OstreamManager clout(std::cout,"prepareLatticeCRAD");
  clout << "Prepare Lattice for Coupled Reaction simualtion..." << std::endl;
  T omega = converter.getLatticeRelaxationFrequencyFromDiffusivity<RADDESCRIPTOR>( D );

  // buffer layer
  CRADlattice->defineDynamics<NoDynamics>(superGeometry, 0);

  // bulk
  CRADlattice->defineDynamics<ReactionADBulkDynamics>(superGeometry.getMaterialIndicator({1, 3, 4}));

  // Material=2 -->bounce back
  CRADlattice->defineDynamics<BounceBack>(superGeometry.getMaterialIndicator({2}));

  // Setting of the boundary conditions, inflow and outflow with Dirichlet and Neumann BCs
  boundary::set<boundary::AdvectionDiffusionDirichlet>(*CRADlattice, superGeometry.getMaterialIndicator({3}));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(*CRADlattice, superGeometry.getMaterialIndicator({4}));
  setZeroGradientBoundary<T, RADDESCRIPTOR>(*CRADlattice, superGeometry.getMaterialIndicator({5}));

  AnalyticalConst3D<T,T> rho0(1.e-8);
  AnalyticalConst3D<T,T> rhoA(rhoA_);
  AnalyticalConst3D<T,T> rhoB(rhoB_);

  T distance2Wall = converter.getConversionFactorLength()/2.;
  std::vector<T> maxVelocity1( 3,0 );
  std::vector<T> maxVelocity2( 3,0 );
  maxVelocity1[1] = converter.getLatticeVelocity(u_max);
  maxVelocity2[1] = converter.getLatticeVelocity(u_max);
  RectanglePoiseuille3D<T> poiseuilleU1( superGeometry, 3, maxVelocity1, distance2Wall, distance2Wall, distance2Wall );
  RectanglePoiseuille3D<T> poiseuilleU2( superGeometry, 4, maxVelocity2, distance2Wall, distance2Wall, distance2Wall );

  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> u( velocity );
  CRADlattice->defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({1, 2, 3, 4, 5 }), u );
  AnalyticalConst3D<T,T> omegaD( omega );
  CRADlattice->defineField<descriptors::OMEGA>(superGeometry.getMaterialIndicator({1, 2, 3, 4, 5 }), omegaD );
  CRADlattice->setParameter<descriptors::OMEGA>(omega);

  // setting actual values for the boundary
  CRADlattice->defineRho( superGeometry, 3, rhoA);
  CRADlattice->iniEquilibrium( superGeometry, 3, rhoA, poiseuilleU1 );
  CRADlattice->defineRho( superGeometry, 4, rhoB);
  CRADlattice->iniEquilibrium( superGeometry, 4, rhoB, poiseuilleU2 );
  CRADlattice->defineField<descriptors::SOURCE>(superGeometry.getMaterialIndicator({1, 2, 3, 4, 5 }), rho0);
  CRADlattice->defineRho( superGeometry.getMaterialIndicator({1, 2, 5}), rho0 );
  CRADlattice->iniEquilibrium( superGeometry.getMaterialIndicator({1, 2, 5}), rho0, u );
  CRADlattice->initialize();

  {
    auto& communicatorCR = CRADlattice->getCommunicator(stage::PreCollide());
    communicatorCR.requestOverlap(3);
    communicatorCR.exchangeRequests();
  }

  clout << "Prepare Lattice for Coupled Reaction simualtion... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( UnitConverter<T,NSDESCRIPTOR> const& converter,
                        SuperLattice<T,NSDESCRIPTOR>& sLatticeNS, int iT,
                        SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  const int iTmaxStart = converter.getLatticeTime( maxPhysT*0.02 );

  if ( iT%100 == 0 && iT <= iTmaxStart ) {

    // Smooth start curve, sinus
    SinusStartScale<T,int> startScale(iTmaxStart, T(1));

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1]= {iT};
    T frac[1]= {};
    startScale( frac,iTvec );
    std::vector<T> maxVelocity1( 3,0 );
    std::vector<T> maxVelocity2( 3,0 );
    maxVelocity1[1] = frac[0]*converter.getLatticeVelocity(u_max);
    maxVelocity2[1] = -frac[0]*converter.getLatticeVelocity(u_max);

    T distance2Wall = converter.getConversionFactorLength()/2.;
    RectanglePoiseuille3D<T> poiseuilleU1( superGeometry, 3, maxVelocity1, distance2Wall, distance2Wall, distance2Wall );
    RectanglePoiseuille3D<T> poiseuilleU2( superGeometry, 4, maxVelocity2, distance2Wall, distance2Wall, distance2Wall );
    sLatticeNS.defineU( superGeometry, 3, poiseuilleU1 );
    sLatticeNS.defineU( superGeometry, 4, poiseuilleU2 );

    sLatticeNS.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);

  }
}

// Output to console and files
void getResultsNS( SuperLattice<T,NSDESCRIPTOR>& sLatticeNS,
                 UnitConverter<T,NSDESCRIPTOR> const& converter,
                 int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer)
{
  OstreamManager clout( std::cout,"getResultsNS" );

  SuperVTMwriter3D<T> vtmWriterNS( "mixer_fluid" );
  SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity( sLatticeNS, converter );
  SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure( sLatticeNS, converter );
  vtmWriterNS.addFunctor( velocity );
  vtmWriterNS.addFunctor( pressure );

  const int  vtkIter  = converter.getLatticeTime( maxPhysT/10 );
  const int  statIter = converter.getLatticeTime( maxPhysT/100 );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid( sLatticeNS );
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank( sLatticeNS );
    vtmWriterNS.write( cuboid );
    vtmWriterNS.write( rank );
    vtmWriterNS.createMasterFile();
  }

  // Writes the ppm files
  if ( iT%vtkIter==0 ) {
    sLatticeNS.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriterNS.write( iT );
  }

  // Writes output on the console
  if ( iT%statIter==0 && iT>=0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLatticeNS.getStatistics().print(iT,converter.getPhysTime(iT));
  }

#ifdef plotOverLine
  if (iT == (converter.getLatticeTime( maxPhysT ) - 1)) {
    // Preparation
    sLatticeNS.setProcessingContext(ProcessingContext::Evaluation);
    Gnuplot<T> csvWriter("flowField");
    int dummy[4];
    T result[4];
    const int nX = superGeometry.getCuboidDecomposition().getMotherCuboid().getNx();  // number of voxels in x-direction
    const int minX = converter.getLatticeLength(lx_inlet);                       // beginning of the mixing channel
    auto mat = new SuperIndicatorIdentity3D<T> (superGeometry.getMaterialIndicator({1,5}));
    const Vector<T,3> normal( 1, 0, 0 );

    for (int iX = minX; iX < nX - 3; ++iX) {
      const T x = converter.getPhysLength(iX);
      const Vector<T,3> center( x, 0, 0 );
      auto circle = new SuperIndicatorFfromIndicatorF3D<T>(
        std::shared_ptr<IndicatorF3D<T>>(new IndicatorCircle3D<T> ( center, normal, 0.0003)), superGeometry);
      SuperIndicatorMultiplication3D<T> plane(circle, mat);

      // Average concentration: use average functor
      SuperAverage3D<T> (velocity, plane).operator()(result, dummy);
      const T avVel = result[0];
      SuperAverage3D<T> (pressure, plane).operator()(result, dummy);
      const T avPres = result[0];

      csvWriter.setData(x, {avVel, avPres});
    }
    csvWriter.writePNG();
  }
#endif
}

void getResultsCRAD( std::vector<SuperLattice<T, RADDESCRIPTOR>*>& adlattices,
                     UnitConverter<T,NSDESCRIPTOR> const& converter,
                     UnitConverter<T,RADDESCRIPTOR> const& converterAD,
                     int iT,
                     SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer)
{
  OstreamManager clout( std::cout,"getResultsCRAD" );

  SuperVTMwriter3D<T> vtmWriterCRAD( "mixer_wholeReaction" );

  // insert concentrations and velocoties for each component into vector and vtk
    std::vector<SuperLatticeDensity3D<T, RADDESCRIPTOR>*> densities;
    std::vector<SuperLatticePhysViscosity3D<T, RADDESCRIPTOR>*> diffs;
    std::vector<std::string> names = {"A","B","P"};

  for (int i = 0; i<3;i++){
      densities.emplace_back( new SuperLatticeDensity3D<T, RADDESCRIPTOR> (*adlattices[i]));
      densities[i]->getName() = "Concentration " + names[i];
      vtmWriterCRAD.addFunctor(*densities[i]);
      diffs.emplace_back( new SuperLatticePhysViscosity3D<T, RADDESCRIPTOR> (*adlattices[i], converterAD));
      diffs[i]->getName() = "Diffusisvity " + names[i];
      vtmWriterCRAD.addFunctor(*diffs[i]);
    }

  const int  vtkIter  = converter.getLatticeTime( maxPhysTAD/10 );
  const int  statIter = converter.getLatticeTime( maxPhysTAD/100 );

  if ( iT==0 ) {
    vtmWriterCRAD.createMasterFile();
  }

  // Writes the ppm files
  if ( iT%vtkIter==0 ) {
    adlattices[0]->setProcessingContext(ProcessingContext::Evaluation);
    adlattices[1]->setProcessingContext(ProcessingContext::Evaluation);
    adlattices[2]->setProcessingContext(ProcessingContext::Evaluation);
    vtmWriterCRAD.write( iT );
  }

  // Writes output on the console
  if ( iT%statIter==0 && iT>=0 ) {
    // Lattice statistics console output
    adlattices[1]->getStatistics().print(iT,converter.getPhysTime(iT));
  }

  #ifdef plotOverLine
  // Compute mixing quantities along central axis
  if (iT == (converter.getLatticeTime( maxPhysT ) - 1)) {
    adlattices[0]->setProcessingContext(ProcessingContext::Evaluation);
    adlattices[1]->setProcessingContext(ProcessingContext::Evaluation);
    adlattices[2]->setProcessingContext(ProcessingContext::Evaluation);
    // Preparation
    Gnuplot<T> csvWriter("mixing");
    int dummy[4];
    T result[4];
    const int nX = superGeometry.getCuboidDecomposition().getMotherCuboid().getNx();  // number of voxels in x-direction
    const int minX = 3;                      // beginning of the mixing channel
    auto mat = new SuperIndicatorIdentity3D<T> (superGeometry.getMaterialIndicator({1/*,5*/}));
    const Vector<T,3> normal( 1, 0, 0 );

    for (int iX = minX; iX < nX - 3; ++iX) {
      const T x = converter.getPhysLength(iX);
      const Vector<T,3> center( x, 0, 0 );
      auto circle = new SuperIndicatorFfromIndicatorF3D<T>(
        std::shared_ptr<IndicatorF3D<T>>(new IndicatorCircle3D<T> ( center, normal, 0.0003)), superGeometry);
      SuperIndicatorMultiplication3D<T> plane(circle, mat);

      // Average concentration: use average functor
      SuperAverage3D<T> (densities[2], plane).operator()(result, dummy);
      const T avConc = result[0];
      if (iX == nX - 4) {  // at outlet
        clout << "av. conc. at outlet = " << avConc << std::endl;
      }
      csvWriter.setData(x, avConc);
    }
    csvWriter.writePNG();
  }
#endif
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  UnitConverterFromResolutionAndLatticeVelocity<T,NSDESCRIPTOR> converter(
    (T)   N,               // resolution
    (T)   0.05,            // maximal lattice velocity
    (T)   hydrInlet,       // charPhysLength: reference length of simulation geometry
    (T)   u_mean,          // charPhysVelocity
    (T)   nu,              // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   density          // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("laminarReactiveTmixer");
  clout << "Maximal simulated time for NSE: " << maxPhysT << std::endl;

  UnitConverterFromResolutionAndLatticeVelocity<T,RADDESCRIPTOR> converterAD(
    (T)   N,               // resolution
    (T)   0.05,            // maximal lattice velocity
    (T)   hydrInlet,       // charPhysLength: reference length of simulation geometry
    (T)   u_mean,          // charPhysVelocity
    (T)   nu,              // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   density          // physDensity: physical density in __kg / m^3__
  );

  // === 2nd Step: Prepare Geometry ===
  T dx = converter.getConversionFactorLength();
  Vector<T,3> extendReactor( lx + 2*dx, ly_outlet + 2*dx, lz + 2*dx );
  Vector<T,3> originReactor( -dx, -(ly_outlet + 2*dx)/2, -(lz + 2*dx)/2 );
  std::shared_ptr<IndicatorCuboid3D<T> > reactor( new IndicatorCuboid3D<T> ( extendReactor, originReactor ));

  Vector<T,3> extendInflow( lx_inlet + 2*dx, ly + 2*dx, lz + 2*dx );
  Vector<T,3> originInflow( -dx, -(ly + 2*dx)/2, -(lz + 2*dx)/2 );
  std::shared_ptr<IndicatorCuboid3D<T> > inflow( new IndicatorCuboid3D<T> ( extendInflow, originInflow ));

  IndicatorIdentity3D<T> mixer(reactor + inflow);

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition3D<T> cuboidDecomposition( mixer, converter.getConversionFactorLength(), noOfCuboids, "volume" );

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

  SuperGeometry<T,3> superGeometry( cuboidDecomposition, loadBalancer, 3 );

  prepareGeometry( converter, superGeometry );

  SuperLattice<T, NSDESCRIPTOR> sLatticeNS( superGeometry );

  prepareLatticeNS( converter, sLatticeNS, superGeometry );

  SuperLattice<T, RADDESCRIPTOR> CRADlattice( superGeometry );

  // initiate reference lattice for coupling, this lattice will be the first component
  std::vector<SuperLattice<T, RADDESCRIPTOR>*> adlattices;
  adlattices.push_back(&CRADlattice);

  // initiate further lattices for the remaining components
  for (int i =1; i<3;i++){
    adlattices.emplace_back( new SuperLattice<T, RADDESCRIPTOR> (superGeometry));
  }
  // partners are used for coupling; coupled with adlattice[0]
  std::vector<SuperLattice<T, RADDESCRIPTOR>*> partners;
  for (int i =1; i<3; i++){
    partners.emplace_back(adlattices[i]);
  }

  // prepare all lattices with the respective concentration
  T rhoA = 1.e-8;
  T rhoB = 1.e-8;
  T D = D_P;
  for(int i = 0; i<3; i++){
    if ( i == 0 ) {rhoA = 1.; D = D_A;}
    if ( i == 1 ) {rhoA = 1.e-8; rhoB = 0.001; D = D_B;}
    if ( i == 2 ) {rhoA = 1.e-8; rhoB = 1.e-8; D = D_P;}
    prepareLatticeCRAD(sLatticeNS, adlattices[i], superGeometry, converter, rhoA, rhoB, D);
  }
  sLatticeNS.statisticsOff();
  adlattices[0] -> statisticsOff();
  adlattices[1] -> statisticsOff();
  adlattices[2] -> statisticsOff();
  T omegaA = converter.getLatticeRelaxationFrequencyFromDiffusivity<RADDESCRIPTOR>( D_A );
  T omegaB = converter.getLatticeRelaxationFrequencyFromDiffusivity<RADDESCRIPTOR>( D_B );
  T omegaP = converter.getLatticeRelaxationFrequencyFromDiffusivity<RADDESCRIPTOR>( D_P );

  clout << "Tau A: " << 1./omegaA << std::endl;
  clout << "Tau B: " << 1./omegaB << std::endl;
  clout << "Tau P: " << 1./omegaP << std::endl;

  /// === 3.1 Step: Prepare Coupling ===
  SuperLatticeCoupling coupling(
    LESReactionCoupling<T,3>{},
    names::NavierStokes{}, sLatticeNS,
    names::Concentration0{}, *adlattices[0],
    names::Concentration1{}, *adlattices[1],
    names::Concentration2{}, *adlattices[2]);
  coupling.setParameter<LESReactionCoupling<T,3>::LATTICE_REACTION_COEFF>(converter.getConversionFactorTime()*1.e5);
  coupling.setParameter<LESReactionCoupling<T,3>::STOCH_COEFF>({-1, -1., 1.});
  coupling.setParameter<LESReactionCoupling<T,3>::REACTION_ORDER>({1., 1., 0.});
  coupling.setParameter<LESReactionCoupling<T,3>::SMAGORINSKY_PREFACTOR>(0.2);
  coupling.setParameter<LESReactionCoupling<T,3>::SCHMIDT>({0.03, 0.03, 0.03});
  coupling.setParameter<LESReactionCoupling<T,3>::OMEGA_NSE>(1./converter.getLatticeRelaxationFrequency());
  coupling.setParameter<LESReactionCoupling<T,3>::OMEGAS_ADE>({omegaA, omegaB, omegaP});

  const int  statIterNS = converter.getLatticeTime( maxPhysT/100 );
  const int  statIterAD = converter.getLatticeTime( maxPhysTAD/100 );

  clout << "starting Coupled LES - Reactive ADE simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    if(iT%statIterNS == 0){
      sLatticeNS.statisticsOn();
    }
    if(iT == 0){
      adlattices[1] -> statisticsOn();
    }
    setBoundaryValues( converter, sLatticeNS, iT, superGeometry );
    sLatticeNS.collideAndStream();
    getResultsNS( sLatticeNS, converter, iT, superGeometry, timer );
    if(iT == 0){getResultsCRAD( adlattices, converter, converterAD, iT, superGeometry, timer );}
    if(iT >= converter.getLatticeTime( maxPhysT - maxPhysTAD )){
      if(iT%statIterAD == 0){
      adlattices[1] -> statisticsOn();
      }
      if(iT%100 == 0){
        coupling.execute();
      }
      for (int i = 0; i<3; i++){
        adlattices[i]->collideAndStream();
      }
      getResultsCRAD( adlattices, converter, converterAD, iT, superGeometry, timer );
    }
    if(iT%statIterNS == 0){
      sLatticeNS.statisticsOff();
    }
    if(iT%statIterAD == 0){
      adlattices[1] -> statisticsOff();
    }
  }

  timer.stop();
  timer.printSummary();
}

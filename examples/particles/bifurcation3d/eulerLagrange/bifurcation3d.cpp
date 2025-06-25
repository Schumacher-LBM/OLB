/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C)
 *  2023      Nicolas Hafen, Frantisek Prinz, Mathias J. Krause
 *  2011-2016 Thomas Henn, Mathias J. Krause, Marie-Luise Maier
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

/* bifurcation3d.cpp:
 * This example examines a steady particulate flow past a bifurcation. At the inlet,
 * an inflow condition with grad_n u = 0 and rho = 1 is implemented.
 * At both outlets, a Poiseuille profile is imposed on the velocity.
 * After a start time, particles are put into the bifurcation at the
 * inlet and experience a stokes drag force.
 *
 * A publication using the same geometry can be found here:
 * http://link.springer.com/chapter/10.1007/978-3-642-36961-2_5
 *  *
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace particles;
using namespace particles::subgrid;
using namespace particles::communication;
using namespace particles::dynamics;
using namespace particles::creators;
using namespace particles::io;

using T = FLOATING_POINT_TYPE;
typedef D3Q19<> DESCRIPTOR;
typedef SubgridParticle3D PARTICLETYPE;

#define BOUZIDI

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const T Re = 50;                    // Reynolds number
int N = 19;                         // resolution of the model
const T radius = 1.5e-4;            // particles radius
const T partRho = 998.2;            // particles density

const T fluidMaxPhysT = T( 5 );     // max. fluid simulation time in s, SI unit
const T particleMaxPhysT = T( 20 ); // max. particle simulation time in s, SI unit

std::size_t noOfParticles = 1000;   // total number of inserted particles

//Set capture method:
// materialCapture: based on material number
// wallCapture:     based on more accurate stl description
typedef enum {materialCapture, wallCapture} ParticleDynamicsSetup;
const ParticleDynamicsSetup particleDynamicsSetup = wallCapture;

// center of inflow and outflow regions [m]
Vector<T, 3> inletCenter( T(), T(), 0.0786395 );
Vector<T, 3> outletCenter0( -0.0235929682287551, -0.000052820468762797,
                            -0.021445708949909 );
Vector<T, 3> outletCenter1( 0.0233643529416147, 0.00000212439067050152,
                            -0.0211994104877918 );

// radii of inflow and outflow regions [m]
T inletRadius = 0.00999839;
T outletRadius0 = 0.007927;
T outletRadius1 = 0.00787134;

// normals of inflow and outflow regions
Vector<T, 3> inletNormal( T(), T(), T( -1 ) );
Vector<T, 3> outletNormal0( 0.505126, -0.04177, 0.862034 );
Vector<T, 3> outletNormal1( -0.483331, -0.0102764, 0.875377 );

//Ensure that parallel mode is used
#ifdef PARALLEL_MODE_MPI

void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      IndicatorF3D<T>& indicator, STLreader<T>& stlReader,
                      SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout, "prepareGeometry" );

  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0, 2, indicator );
  superGeometry.rename( 2, 1, stlReader );

  superGeometry.clean();

  // rename the material at the inlet
  IndicatorCircle3D<T> inletCircle( inletCenter, inletNormal,
                                    inletRadius );
  IndicatorCylinder3D<T> inlet( inletCircle,
                                2 * converter.getPhysDeltaX() );
  superGeometry.rename( 2, 3, 1, inlet );

  // rename the material at the outlet0
  IndicatorCircle3D<T> outletCircle0( outletCenter0, outletNormal0,
                                      0.95 * outletRadius0 );
  IndicatorCylinder3D<T> outlet0( outletCircle0,
                                  4 * converter.getPhysDeltaX() );
  superGeometry.rename( 2, 4, outlet0 );

  // rename the material at the outlet1
  IndicatorCircle3D<T> outletCircle1( outletCenter1, outletNormal1,
                                      0.95 * outletRadius1 );
  IndicatorCylinder3D<T> outlet1( outletCircle1,
                                  4 * converter.getPhysDeltaX() );
  superGeometry.rename( 2, 5, outlet1 );

  superGeometry.clean();
  #ifndef BOUZIDI
  superGeometry.innerClean( true );
  #endif
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

void prepareLattice( SuperLattice<T, DESCRIPTOR>& superLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     STLreader<T>& stlReader, SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout, "prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  superLattice.defineDynamics<BGKdynamics>(superGeometry, 1);

  #ifdef BOUZIDI
  setBouzidiBoundary(superLattice, superGeometry, 2, stlReader);
  #else
  boundary::set<boundary::BounceBack>(superLattice, superGeometry, 2);
  #endif

  // Material=3 -->bulk dynamics (inflow)
  superLattice.defineDynamics<BGKdynamics>(superGeometry, 3);

  // Material=4 -->bulk dynamics (outflow)
  superLattice.defineDynamics<BGKdynamics>(superGeometry, 4);
  superLattice.defineDynamics<BGKdynamics>(superGeometry, 5);

  // Setting of the boundary conditions
  boundary::set<boundary::InterpolatedPressure>(superLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedVelocity>(superLattice, superGeometry, 4);
  boundary::set<boundary::InterpolatedVelocity>(superLattice, superGeometry, 5);

  superLattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

// Generates a slowly increasing sinusoidal inflow for the first iTMax timesteps
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& superLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter, int iT, T maxPhysT,
                        SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout, "setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( 0.8*maxPhysT );
  int iTperiod = 100;  // amount of timesteps when new boundary conditions are reset

  if ( iT == 0 ) {

    AnalyticalConst3D<T, T> rhoF( 1 );
    std::vector<T> velocity( 3, T() );
    AnalyticalConst3D<T, T> uF( velocity );

    superLattice.iniEquilibrium( superGeometry, 1, rhoF, uF );
    superLattice.iniEquilibrium( superGeometry, 2, rhoF, uF );
    superLattice.iniEquilibrium( superGeometry, 3, rhoF, uF );
    superLattice.iniEquilibrium( superGeometry, 4, rhoF, uF );
    superLattice.iniEquilibrium( superGeometry, 5, rhoF, uF );

    superLattice.defineRhoU( superGeometry, 1, rhoF, uF );
    superLattice.defineRhoU( superGeometry, 2, rhoF, uF );
    superLattice.defineRhoU( superGeometry, 3, rhoF, uF );
    superLattice.defineRhoU( superGeometry, 4, rhoF, uF );
    superLattice.defineRhoU( superGeometry, 5, rhoF, uF );

    // Make the lattice ready for simulation
    superLattice.initialize();
  }

  else if ( iT <= iTmaxStart && iT % iTperiod == 0 ) {
    SinusStartScale<T, int> startScale( iTmaxStart, T( 1 ) );
    int iTvec[1] = { iT };
    T frac[1] = { T( 0 ) };
    startScale( frac, iTvec );
    T maxVelocity = frac[0] * converter.getCharLatticeVelocity() * 3. / 4.
                    * util::pow( inletRadius, 2 ) / util::pow( outletRadius0, 2 );

    CirclePoiseuille3D<T> poiseuilleU4( outletCenter0[0], outletCenter0[1],
                                        outletCenter0[2], outletNormal0[0],
                                        outletNormal0[1], outletNormal0[2],
                                        outletRadius0 * 0.95, -maxVelocity );

    CirclePoiseuille3D<T> poiseuilleU5( outletCenter1[0], outletCenter1[1],
                                        outletCenter1[2], outletNormal1[0],
                                        outletNormal1[1], outletNormal1[2],
                                        outletRadius1 * 0.95, -maxVelocity );

    superLattice.defineU( superGeometry, 4, poiseuilleU4 );
    superLattice.defineU( superGeometry, 5, poiseuilleU5 );
  }
}

//Prepare particles
void prepareParticles( UnitConverter<T,DESCRIPTOR> const& converter,
  SuperParticleSystem<T,PARTICLETYPE>& superParticleSystem,
  SolidBoundary<T,3>& wall,
  SuperGeometry<T,3>& superGeometry,
  ParticleDynamicsSetup particleDynamicsSetup,
  Randomizer<T>& randomizer )
{
  OstreamManager clout( std::cout, "prepareParticles" );
  clout << "Prepare Particles ..." << std::endl;

  //Create material indicators for particle boundaries
  std::vector<int> wallMaterials {2};
  std::shared_ptr<SuperIndicatorMaterial<T,3>> wallMaterialIndicator =
    std::make_shared<SuperIndicatorMaterial<T,3>>(superGeometry, wallMaterials);
  std::vector<int> outletMaterials {4,5};
  std::shared_ptr<SuperIndicatorMaterial<T,3>> outletMaterialIndicator =
    std::make_shared<SuperIndicatorMaterial<T,3>>(superGeometry, outletMaterials);

  //Add selected particle dynamics
  if (particleDynamicsSetup==wallCapture){
    //Create verlet dynamics with material aware wall capture
    superParticleSystem.defineDynamics<
      VerletParticleDynamicsMaterialAwareWallCaptureAndEscape<T,PARTICLETYPE>>(
        wall, wallMaterialIndicator, outletMaterialIndicator);
  } else {
    //Create verlet dynamics with material capture
    superParticleSystem.defineDynamics<
      VerletParticleDynamicsMaterialCaptureAndEscape<T,PARTICLETYPE>>(
        wallMaterialIndicator, outletMaterialIndicator);
  }

  // particles generation at inlet3
  Vector<T, 3> c( inletCenter );
  c[2] = 0.074;
  IndicatorCircle3D<T> inflowCircle( c, inletNormal, inletRadius -
                                     converter.getPhysDeltaX() * 2.5 );
  IndicatorCylinder3D<T> inletCylinder( inflowCircle, 0.01 *
                                        converter.getPhysDeltaX() );

  //Add particles
  addParticles( superParticleSystem, inletCylinder, partRho, radius, noOfParticles, randomizer );

  //Print super particle system summary
  superParticleSystem.print();
  clout << "Prepare Particles ... OK" << std::endl;
}





// Computes the pressure drop between voxels before and after the cylinder
bool getResults( SuperLattice<T, DESCRIPTOR>& superLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT, int iTperiod,
                 SuperGeometry<T,3>& superGeometry,
                 Timer<double>& fluidTimer, STLreader<T>& stlReader, bool fluidExists,
                 SuperParticleSystem<T,PARTICLETYPE>& superParticleSystem,
                 Timer<double>& particleTimer )
{

  OstreamManager clout( std::cout, "getResults" );
  SuperVTMwriter3D<T> vtmWriter( "bifurcation3d" );
  SuperVTMwriter3D<T> vtmWriterStartTime( "startingTimeBifurcation3d" );

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( superLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( superLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  vtmWriterStartTime.addFunctor( velocity );
  vtmWriterStartTime.addFunctor( pressure );

  //Create VTK writer for particles
  VTUwriter<T,PARTICLETYPE,true> superParticleWriter("particles_master",false,false);

  //Create functors
  SuperParticleGroupedFieldF<T,PARTICLETYPE,GENERAL,POSITION>
    particlePosition( superParticleSystem );
  SuperParticleGroupedFieldF<T,PARTICLETYPE,PHYSPROPERTIES,MASS>
    particleMass( superParticleSystem );
  SuperParticleGroupedFieldF<T,PARTICLETYPE,PHYSPROPERTIES,RADIUS>
    particleRadius( superParticleSystem );
  SuperParticleGroupedFieldF<T,PARTICLETYPE,MOBILITY,VELOCITY>
    particleVelocity( superParticleSystem );
  SuperParticleGroupedFieldF<T,PARTICLETYPE,FORCING,FORCE>
    particleForce( superParticleSystem );
  SuperParticleGroupedFieldF<T,PARTICLETYPE,DYNBEHAVIOUR,ACTIVE>
    particleActivity( superParticleSystem );
  superParticleWriter.addFunctor(particlePosition,"Position");
  superParticleWriter.addFunctor(particleMass,"Mass");
  superParticleWriter.addFunctor(particleRadius,"Radius");
  superParticleWriter.addFunctor(particleVelocity,"Velocity");
  superParticleWriter.addFunctor(particleForce,"Force");
  superParticleWriter.addFunctor(particleActivity,"Active");

  std::size_t fluidMaxT = converter.getLatticeTime( fluidMaxPhysT );

  if ( iT == 0 ) {
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( superLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( superLattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
    superParticleWriter.createMasterFile();
    vtmWriterStartTime.createMasterFile();


    // Print some output of the chosen simulation setup
    clout << "N=" << N <<"; maxTimeSteps(fluid)="
          << converter.getLatticeTime( fluidMaxPhysT ) << "; noOfCuboid="
          << superGeometry.getCuboidDecomposition().size() << "; Re=" << Re
          <<  "; noOfParticles=" << noOfParticles << "; maxTimeSteps(particle)="
          << converter.getLatticeTime( particleMaxPhysT )
          << "; St=" << ( 2.*partRho*radius*radius*converter.getCharPhysVelocity() ) / ( 9.*converter.getPhysViscosity()*converter.getPhysDensity()*converter.getCharPhysLength() ) << std::endl;
  }

  // Writes the .vtk and .gif files
  if ( iT % iTperiod == 0 ) {
    if ( !fluidExists && iT <= fluidMaxT ) {
      vtmWriterStartTime.write(iT);
      SuperEuklidNorm3D<T> normVel( velocity );
      BlockReduction3D2D<T> planeReduction( normVel, {0, -1, 0}, 600, BlockDataSyncMode::ReduceOnly );
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }
    if (iT > fluidMaxT) {
      // only write .vtk-files after the fluid calculation is finished
      vtmWriter.write(iT - fluidMaxT);
    }
  }

  // Writes output on the console for the fluid phase
  if (iT < converter.getLatticeTime( fluidMaxPhysT ) && iT%iTperiod == 0 ) {

    // Timer statics
    fluidTimer.update( iT );
    fluidTimer.printStep();

    // Lattice statistics
    superLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );

    // Flux at the inlet and outlet regions
    const std::vector<int> materials = { 1, 3, 4, 5 };

    IndicatorCircle3D<T> inlet(
      inletCenter + 2. * converter.getPhysDeltaX() * inletNormal,
      inletNormal, inletRadius + 2. * converter.getPhysDeltaX() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow( superLattice, converter, superGeometry, inlet, materials );
    vFluxInflow.print( "inflow", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxInflow( superLattice, converter, superGeometry, inlet, materials );
    pFluxInflow.print( "inflow", "N", "Pa" );

    IndicatorCircle3D<T> outlet0(
      outletCenter0 + 2. * converter.getPhysDeltaX() * outletNormal0,
      outletNormal0, outletRadius0 + 2. * converter.getPhysDeltaX() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow0( superLattice, converter, superGeometry, outlet0, materials );
    vFluxOutflow0.print( "outflow0", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow0( superLattice, converter, superGeometry, outlet0, materials );
    pFluxOutflow0.print( "outflow0", "N", "Pa" );

    IndicatorCircle3D<T> outlet1(
      outletCenter1 + 2. * converter.getPhysDeltaX() * outletNormal1,
      outletNormal1, outletRadius1 + 2. * converter.getPhysDeltaX() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow1( superLattice, converter, superGeometry, outlet1, materials );
    vFluxOutflow1.print( "outflow1", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow1( superLattice, converter, superGeometry, outlet1, materials );
    pFluxOutflow1.print( "outflow1", "N", "Pa" );
  }

  // Writes output on the console for the fluid phase
  if ( iT >= converter.getLatticeTime( fluidMaxPhysT ) &&
       (iT%iTperiod == 0 || iT == converter.getLatticeTime( fluidMaxPhysT )) ) {
    // advance particle timer
    particleTimer.print( iT - fluidMaxT );

    //Purge invalid particles (delete invalidated particles)
    purgeInvalidParticles<T,PARTICLETYPE>( superParticleSystem );

    //Define materials for capture rate
    std::vector<int> materialsOutout {4,5};
    SuperIndicatorMaterial<T,3> materialIndicatorOutput (superGeometry, materialsOutout);

    //Perform capture statistics
    std::size_t noActive;
    captureStatistics( superParticleSystem, materialIndicatorOutput, noActive );

    // only write .vtk-files after the fluid calculation is finished
    superParticleWriter.write( iT - fluidMaxT );

    // true as long as certain amount of active particles
    if ( noActive < 0.001 * noOfParticles
         && iT > 0.9*converter.getLatticeTime( fluidMaxPhysT + particleMaxPhysT ) ) {
      return false;
    }
    //Additional criterion added 02.02.23
    if ( noActive==0 ) {
      return false;
    }
  }
  return true;
}


#endif //PARALLEL_MODE_MPI





int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===

  initialize( &argc, &argv );

  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout, "main" );

  //Create randomizer
  Randomizer<T> randomizer;

// Ensure that parallel mode is used
#ifdef PARALLEL_MODE_MPI

  // Input treatment
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='h') {
      OstreamManager clout( std::cout,"help" );
      clout << std::endl << "Optional arguments: [N] [noP]" << std::endl;
      return 0;
    }
    N = atoi(argv[1]);
    if (N<1) {
      std::cerr << "ERROR: Resolution was set to N < 1" << std::endl;
      return 1;
    }
  }
  if (argc>2) { noOfParticles = atoi(argv[2]);}
  if (argc>3) {
    std::cerr << "ERROR: too many arguments!" << std::endl;
    return 1;
  }

  // Creat unit converter
  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> const converter(
    int {N}, // resolution: number of voxels per charPhysL
    (T)   0.557646,                    // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   inletRadius*2.,              // charPhysLength: reference length of simulation geometry
    (T)   Re*1.5e-5/( inletRadius*2 ), // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.5e-5,                      // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.225                        // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("converter");

  // === 2nd Step: Prepare Geometry ===
  STLreader<T> stlReader( "../bifurcation3d.stl", converter.getPhysDeltaX() );
  IndicatorLayer3D<T> extendedDomain( stlReader,
                                      converter.getPhysDeltaX() );

  // Create solid wall
  const unsigned latticeMaterial = 2; //Material number of wall
  const unsigned contactMaterial = 0; //Material identifier (only relevant for contact model)
  SolidBoundary<T,3> wall( std::make_unique<IndicInverse<T,3>>(stlReader),
                           latticeMaterial, contactMaterial );

  // Instantiation of an empty cuboidDecomposition
  int noOfCuboids = util::max( 16, 4 * singleton::mpi().getSize() );

  CuboidDecomposition3D<T> cuboidDecomposition( extendedDomain, converter.getPhysDeltaX(),
                                      noOfCuboids );

  // Instantiation of an empty loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

  // Default instantiation of superGeometry
  SuperGeometry<T,3> superGeometry( cuboidDecomposition, loadBalancer, 3 );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice<T, DESCRIPTOR> superLattice( superGeometry );

  //prepareLattice and setBoundaryConditions
  prepareLattice( superLattice, converter, stlReader, superGeometry );

  // === 3.1 Step: Particles ===

  // SuperParticleSystems
  SuperParticleSystem<T,PARTICLETYPE> superParticleSystem(superGeometry);

  //Create particle manager
  ParticleManager<T,DESCRIPTOR,PARTICLETYPE> particleManager(
    superParticleSystem, superGeometry, superLattice, converter);

  //Prepare particles
  prepareParticles( converter, superParticleSystem, wall, superGeometry,
    particleDynamicsSetup, randomizer );


  // === 4th Step: Main Loop with Timer ===

  Timer<double> fluidTimer( converter.getLatticeTime( fluidMaxPhysT ),
                            superGeometry.getStatistics().getNvoxel() );

  Timer<double> particleTimer( converter.getLatticeTime( particleMaxPhysT ),
                               noOfParticles );
  fluidTimer.start();

  std::size_t iT = 0;
  // amount of timesteps when getResults rewrites data
  int iTperiod = converter.getLatticeTime( .2 );

  bool fluidExists = true;

  // checks whether there is already data of the fluid from an earlier calculation
  if ( !( superLattice.load( "fluidSolution_N"+std::to_string(N) ) ) ) {

    fluidExists = false;

    // if there is no data available, it is generated
    for ( ; iT <= converter.getLatticeTime( fluidMaxPhysT ); ++iT ) {

      // during run up time boundary values are set, collide and stream step,
      // results of fluid, afterwards only particles are simulated
      setBoundaryValues( superLattice, converter, iT, fluidMaxPhysT, superGeometry );
      superLattice.collideAndStream();

      getResults( superLattice, converter, iT, iTperiod, superGeometry,
                  fluidTimer, stlReader, fluidExists,
                  superParticleSystem, particleTimer );
    }

    fluidTimer.stop();
    fluidTimer.printSummary();

    superLattice.communicate();
    // calculated results are written in a file
    superLattice.save( "fluidSolution_N"+std::to_string(N) );
  }

  // if there exists already data of the fluid from an earlier calculation, this is used
  else {

    iT = converter.getLatticeTime( fluidMaxPhysT );
    getResults( superLattice, converter, iT,
                iTperiod, superGeometry, fluidTimer, stlReader, fluidExists,
                superParticleSystem, particleTimer);

  }

  // initialize particle velocity
  initializeParticleVelocity( superLattice, superGeometry, converter, superParticleSystem );

  // after the fluid calculation, particle simulation starts
  particleTimer.start();

  for ( ; iT <= converter.getLatticeTime( fluidMaxPhysT + particleMaxPhysT ); ++iT ) {
    // particles simulation starts after run up time is over
    particleManager.execute<
      couple_lattice_to_particles<T,DESCRIPTOR,PARTICLETYPE>,
      process_dynamics<T,PARTICLETYPE>,
      update_particle_core_distribution<T,PARTICLETYPE>
    >();

    if ( !getResults( superLattice, converter, iT,
                      iTperiod, superGeometry, fluidTimer,
                      stlReader, fluidExists,
                      superParticleSystem,
                      particleTimer) ){
      break;
    }
  }
  particleTimer.stop();
  particleTimer.printSummary();
#else
  std::cerr << std::endl
            << "ERROR: Subgrid particles can only be used with MPI!"
            << std::endl << std::endl;
#endif //PARALLEL_MODE_MPI
}

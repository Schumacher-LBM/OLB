/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2014 Mathias J. Krause
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

/* aorta3d.cpp:
 * In this example the fluid flow through a bifurcation is
 * simulated. The geometry is obtained from a mesh in stl-format.
 * With Bouzidi boundary conditions the curved boundary is
 * adequately mapped and initialized fully automatically. As
 * dynamics a Smagorinsky turbulent BGK model is used to stabilize
 * the simulation for low resolutions. As output the flux at the
 * inflow and outflow region is computed. The wall stress can be
 * visualized on the stl Mesh with the Mesh.pvd file in paraview.
 * The results has been validated by comparison with other results
 * obtained with FEM and FVM.
 */

#include <olb.h>


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D3Q19<>;
using BulkDynamics = SmagorinskyBGKdynamics<T,DESCRIPTOR>;

//simulation parameters
const int N = 40;             // resolution of the model
const int M = 20;             // time discretization refinement
const bool bouzidiOn = true;  // choice of boundary condition
const T maxPhysT = 2.;        // max. simulation time in s, SI unit


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader, SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,stlReader );

  superGeometry.clean();

  // Set material number for inflow
  IndicatorCircle3D<T> inflow(  0.218125,0.249987,0.0234818, 0., 1.,0., 0.0112342 );
  IndicatorCylinder3D<T> layerInflow( inflow, 2.*converter.getPhysDeltaX() );
  superGeometry.rename( 2,3,1,layerInflow );

  // Set material number for outflow0
  //IndicatorCircle3D<T> outflow0(0.2053696,0.0900099,0.0346537,  2.5522,5.0294,-1.5237, 0.0054686 );
  IndicatorCircle3D<T> outflow0( 0.2053696,0.0900099,0.0346537, 0.,-1.,0., 0.0054686 );
  IndicatorCylinder3D<T> layerOutflow0( outflow0, 2.*converter.getPhysDeltaX() );
  superGeometry.rename( 2,4,1,layerOutflow0 );

  // Set material number for outflow1
  //IndicatorCircle3D<T> outflow1(0.2388403,0.0900099,0.0343228, -1.5129,5.1039,-2.8431, 0.0058006 );
  IndicatorCircle3D<T> outflow1( 0.2388403,0.0900099,0.0343228, 0.,-1.,0., 0.0058006 );
  IndicatorCylinder3D<T> layerOutflow1( outflow1, 2.*converter.getPhysDeltaX() );
  superGeometry.rename( 2,5,1,layerOutflow1 );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean( 3 );
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T, DESCRIPTOR>& lattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     STLreader<T>& stlReader, SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // material=1 --> bulk dynamics
  lattice.defineDynamics<BulkDynamics>(superGeometry, 1);

  if ( bouzidiOn ) {
    // material=2 --> no dynamics + bouzidi zero velocity
    setBouzidiBoundary<T,DESCRIPTOR>(lattice, superGeometry, 2, stlReader);
    // material=3 --> no dynamics + bouzidi velocity (inflow)
    setBouzidiBoundary<T,DESCRIPTOR,BouzidiVelocityPostProcessor>(lattice, superGeometry, 3, stlReader);
  }
  else {
    // material=2 --> bounceBack dynamics
    boundary::set<boundary::BounceBack>(lattice, superGeometry, 2);
    // material=3 --> bulk dynamics + velocity (inflow)
    lattice.defineDynamics<BulkDynamics>(superGeometry, 3);
    boundary::set<boundary::InterpolatedVelocity>(lattice, superGeometry, 3);
  }

  // material=4,5 --> bulk dynamics + pressure (outflow)
  lattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator({4, 5}));
  boundary::set<boundary::InterpolatedPressure>(lattice, superGeometry, 4);
  boundary::set<boundary::InterpolatedPressure>(lattice, superGeometry, 5);

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> uF( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  lattice.defineRhoU( superGeometry.getMaterialIndicator({1, 3, 4, 5}),rhoF,uF );
  lattice.iniEquilibrium( superGeometry.getMaterialIndicator({1, 3, 4, 5}),rhoF,uF );

  lattice.setParameter<descriptors::OMEGA>(omega);
  lattice.setParameter<collision::LES::SMAGORINSKY>(T(0.1));
  // Lattice initialize
  lattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing sinuidal inflow
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                        SuperGeometry<T,3>& superGeometry )
{

  // No of time steps for smooth start-up
  std::size_t iTperiod = converter.getLatticeTime( 0.5 );
  std::size_t iTupdate = 50;

  if ( iT%iTupdate == 0 ) {
    // Smooth start curve, sinus
    SinusStartScale<T,std::size_t> nSinusStartScale( iTperiod,converter.getCharLatticeVelocity() );

    // Creates and sets the Poiseuille inflow profile using functors
    std::size_t iTvec[1]= {iT};
    T maxVelocity[1]= {T()};
    nSinusStartScale( maxVelocity,iTvec );
    CirclePoiseuille3D<T> velocity( superGeometry,3,maxVelocity[0], T() );

    if ( bouzidiOn ) {
      setBouzidiVelocity(sLattice, superGeometry, 3, velocity);
      sLattice.setProcessingContext<Array<descriptors::BOUZIDI_VELOCITY>>(
        ProcessingContext::Simulation);
    }
    else {
      sLattice.defineU(superGeometry, 3, velocity);
      sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);
    }
  }
}

// Computes flux at inflow and outflow
void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                const UnitConverter<T,DESCRIPTOR>& converter,
                std::size_t iT,
                SuperGeometry<T,3>& superGeometry,
                util::Timer<T>& timer,
                STLreader<T>& stlReader,
                VTUsurfaceWriter<T>& vtuWriter)
{
  OstreamManager clout( std::cout,"getResults" );

  const std::size_t vtkIter  = converter.getLatticeTime( .1 );
  const std::size_t statIter = converter.getLatticeTime( .1 );

  if ( iT==0 ) {
    SuperVTMwriter3D<T> vtmWriter("aorta3d");
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization

    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
    vtuWriter.createMasterFile();
  }

  // Writes the vtk files
  if ( iT%vtkIter==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter("aorta3d");
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      SuperLatticePhysPressure3D pressure(sLattice, converter);
      vtmWriter.addFunctor(velocity);
      vtmWriter.addFunctor(pressure);
      task(vtmWriter, iT);
    });

    // Write interpolated wall shear stress on the STL surface
    {
      SuperLatticeDensity3D densityF(sLattice);
      AnalyticalFfromSuperF3D smoothDensityF(densityF);

      SuperLatticeStress3D stressF(sLattice);
      AnalyticalFfromSuperF3D smoothStressF(stressF);

      PhysWallShearStressOnSurface3D<T,DESCRIPTOR> interpolatedWssF(converter, smoothDensityF, smoothStressF, stlReader);
      interpolatedWssF.getName() = "interpolatedWss";

      vtuWriter.addFunctor(interpolatedWssF);

      vtuWriter.write(iT);
    }
  }

  // Writes output on the console
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Flux at the inflow and outflow region
    std::vector<int> materials = { 1, 3, 4, 5 };

    IndicatorCircle3D<T> inflow(  0.218125,0.249987-2.*converter.getPhysDeltaX(),0.0234818, 0., -1.,0., 0.0112342+2*converter.getPhysDeltaX() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow( sLattice, converter, superGeometry, inflow, materials, BlockDataReductionMode::Discrete );
    vFluxInflow.print( "inflow","ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxInflow( sLattice, converter, superGeometry, inflow, materials, BlockDataReductionMode::Discrete );
    pFluxInflow.print( "inflow","N","mmHg" );

    IndicatorCircle3D<T> outflow0( 0.2053696,0.0900099+2.*converter.getPhysDeltaX(),0.0346537, 0.,1.,0., 0.0054686+2*converter.getPhysDeltaX() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow0( sLattice, converter, superGeometry, outflow0, materials, BlockDataReductionMode::Discrete );
    vFluxOutflow0.print( "outflow0","ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow0( sLattice, converter, superGeometry, outflow0, materials, BlockDataReductionMode::Discrete );
    pFluxOutflow0.print( "outflow0","N","mmHg" );

    IndicatorCircle3D<T> outflow1( 0.2388403,0.0900099+2.*converter.getPhysDeltaX(),0.0343228, 0.,1.,0., 0.0058006+2*converter.getPhysDeltaX() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow1( sLattice, converter, superGeometry, outflow1, materials, BlockDataReductionMode::Discrete );
    vFluxOutflow1.print( "outflow1","ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow1( sLattice, converter, superGeometry, outflow1, materials, BlockDataReductionMode::Discrete );
    pFluxOutflow1.print( "outflow1","N","mmHg" );

    if ( bouzidiOn ) {
      SuperLatticeYplus3D<T, DESCRIPTOR> yPlus( sLattice, converter, superGeometry, stlReader, 3 );
      SuperMax3D<T> yPlusMaxF( yPlus, superGeometry, 1 );
      int input[4]= {};
      T yPlusMax[1];
      yPlusMaxF( yPlusMax,input );
      clout << "yPlusMax=" << yPlusMax[0] << std::endl;
    }
  }

  if ( sLattice.getStatistics().getMaxU() > 0.3 ) {
    clout << "PROBLEM uMax=" << sLattice.getStatistics().getMaxU() << std::endl;
    std::exit(0);
  }
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  const UnitConverter<T,DESCRIPTOR> converter(
    (T)   0.02246/N,     // physDeltaX: spacing between two lattice cells in __m__
    (T)   0.02246/(M*N), // physDeltaT: time step in __s__
    (T)   0.02246,       // charPhysLength: reference length of simulation geometry
    (T)   0.45,          // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.003/1055.,   // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1055           // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("aorta3d");

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  STLreader<T> stlReader( "aorta3d.stl", converter.getPhysDeltaX(), 0.001,  olb::RayMode::FastRayZ, true );
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getPhysDeltaX() );

  // Instantiation of a cuboidDecomposition with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = util::min(16*N, 8*singleton::mpi().getSize());
#else
  const int noOfCuboids = 2;
#endif
  CuboidDecomposition3D<T> cuboidDecomposition( extendedDomain, converter.getPhysDeltaX(), noOfCuboids, "volume" );
  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidDecomposition, loadBalancer );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  util::Timer<T> timer1( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer1.start();

  prepareLattice( sLattice, converter, stlReader, superGeometry );

  VTUsurfaceWriter<T> vtuWriter("surface", cuboidDecomposition, loadBalancer);
  vtuWriter.addSTL( stlReader );

  timer1.stop();
  timer1.printSummary();

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT <= converter.getLatticeTime( maxPhysT ); iT++ ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer, stlReader, vtuWriter );
  }

  timer.stop();
  timer.printSummary();
}

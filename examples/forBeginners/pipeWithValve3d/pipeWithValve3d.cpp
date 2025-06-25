/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2025 Fedor Bukreev, Adrian Kummerländer,
 *  Shota Ito, Mathias J. Krause
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

 /* pipeWithValve3d.cpp:
 * This example represents flow in a pipe with rotated valve inside.
 * Pipe and valve geometries are read from STL files.
 */

#include <olb.h>

using namespace olb;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = descriptors::D3Q19<>;

// Parameters for the simulation setup
const int N = 40;                               // Number of cells in the pipe height
const T Re = 20.;                               // Reynolds number
const T maxPhysT = 16.;                         // max. simulation time in s, SI unit
const T L = 0.05/N;                             // latticeL
const T CFL = 0.05;                             // CFL number
Vector<T,3> rotationPoint = {0.15, 0., 0.}; // coordinates of the rotation point for valve
Vector<T,3> rotationAxis = {0., 0., 1.};    // axis of rotation
T angle = 75;                                   // rotation angle in degrees

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                      STLreader<T>& pipe, STLreader<T>& valve, SuperGeometry<T,3>& superGeometry )
{
  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,pipe );
  superGeometry.clean();

  Vector<T,3> origin = superGeometry.getStatistics().getMinPhysR( 2 );
  origin[1] += converter.getPhysDeltaX()/2.;
  origin[2] += converter.getPhysDeltaX()/2.;

  Vector<T,3> extend = superGeometry.getStatistics().getMaxPhysR( 2 );
  extend[1] = extend[1]-origin[1]-converter.getPhysDeltaX()/2.;
  extend[2] = extend[2]-origin[2]-converter.getPhysDeltaX()/2.;

  // Set material number for inflow
  origin[0] = superGeometry.getStatistics().getMinPhysR( 2 )[0]-converter.getPhysDeltaX();
  extend[0] = 2*converter.getPhysDeltaX();
  IndicatorCuboid3D<T> inflow( extend,origin );
  superGeometry.rename( 2,3,inflow );

  // Set material number for outflow
  origin[0] = superGeometry.getStatistics().getMaxPhysR( 2 )[0]-converter.getPhysDeltaX();
  extend[0] = 2*converter.getPhysDeltaX();
  IndicatorCuboid3D<T> outflow( extend,origin );
  superGeometry.rename( 2,4,outflow );

  // Rotate the valve and set material number for it
  T rotationAngle = std::numbers::pi_v<T> / 180. *  angle;
  IndicatorRotate<T,3> valveRot(rotationPoint, rotationAxis, rotationAngle, valve);
  superGeometry.rename( 1,5, valveRot );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();
  superGeometry.print();
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics<BGKdynamics>(superGeometry, 1);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  // Material=3 -->fixed velocity
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);

  // Material=4 -->fixed pressure
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  // Material=5 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 5);

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  AnalyticalConst3D<T,T> uF(T(0), T(0), T(0));

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( superGeometry, 1, rhoF, uF );
  sLattice.iniEquilibrium( superGeometry, 1, rhoF, uF );

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
  int iTupdate = 30;

  if ( iT%iTupdate == 0 && iT <= iTmaxStart ) {
    PolynomialStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );
    std::vector<T> maxVelocity( 3,0 );
    maxVelocity[0] = 2.25*frac[0]*converter.getCharLatticeVelocity();
    T distance2Wall = converter.getPhysDeltaX()/2.;
    CirclePoiseuille3D<T> poiseuilleU(superGeometry, 3, maxVelocity[0], distance2Wall);
    sLattice.defineU( superGeometry, 3, poiseuilleU );

    // Update velocity on GPU
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer)
{

  OstreamManager clout( std::cout,"getResults" );

  const std::size_t vtkIter  = converter.getLatticeTime( .3 );
  const std::size_t statIter = converter.getLatticeTime( .1 );

  SuperVTMwriter3D<T> vtmWriter( "pipeWithValve3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  if ( iT==0 ) {
    vtmWriter.createMasterFile();
  }

  // Writes the vtk files
  if (iT%vtkIter == 0) {
    // Send values from GPU to CPU for evaluation
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  // Writes output on the console
  if ( iT%statIter == 0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  OstreamManager clout( std::cout,"main" );

  UnitConverter<T,DESCRIPTOR> converter(
    (T)   L,                        // physDeltaX: spacing between two lattice cells in [m]
    (T)   CFL*L/0.2,                // physDeltaT: time step in [s]
    (T)   0.5,                      // charPhysLength: reference length of simulation geometry in [m]
    (T)   0.2,                      // charPhysVelocity: highest expected velocity during simulation in [m/s]
    (T)   0.2*0.05/Re,              // physViscosity: physical kinematic viscosity in [m^2/s]
    (T)   1000.0                    // physDensity: physical density in [kg/m^3]
  );
  converter.print();


  // === 2nd Step: Prepare Geometry ===
  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  // Reading the pipe STL
  STLreader<T> pipe( "pipe.stl", converter.getPhysDeltaX(), 0.001);
  // Extending the pipe with 1 cell for outer boundaries
  IndicatorLayer3D<T> extendedDomain( pipe, converter.getPhysDeltaX() );
  // Reading the valve STL
  STLreader<T> valve( "valve.stl", converter.getPhysDeltaX(), 0.001);

  // Instantiation of a cuboidDecomposition with weights
  CuboidDecomposition3D<T> cuboidDecomposition( extendedDomain, converter.getPhysDeltaX(), singleton::mpi().getSize() );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidDecomposition, loadBalancer );

  prepareGeometry( converter, extendedDomain, pipe, valve, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and set boundaryCondition
  prepareLattice( sLattice, converter, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}
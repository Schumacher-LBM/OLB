/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mathias J. Krause, Mingliang Zhong, Stephan Simonis
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

#include <olb.h>

using namespace olb;
using namespace olb::uq;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D3Q19<>;
using BulkDynamics = ConstRhoBGKdynamics<T,DESCRIPTOR>;

extern const T vtkSave;  ///< Interval for writing VTK output (in physical seconds)
extern const T maxPhysT; ///< Maximum physical time for simulation (seconds)

void prepareGeometry( UnitConverter<T, DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator, SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid and boundary
  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,{1,1,1} );

  T eps = converter.getPhysDeltaX();
  Vector<T,3> origin( -eps, converter.getCharPhysLength() - eps, -eps );
  Vector<T,3> extend( converter.getCharPhysLength() + 2*eps, 2*eps, converter.getCharPhysLength() + 2*eps );
  IndicatorCuboid3D<T> lid( extend,origin );

  superGeometry.rename( 2,3,1,lid );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice( UnitConverter<T, DESCRIPTOR> const& converter,
                     SuperLattice<T,DESCRIPTOR>& lattice,
                     SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  lattice.defineDynamics<BulkDynamics>(superGeometry, 1);

  // Material=2,3 -->bulk dynamics, velocity boundary
  boundary::set<boundary::InterpolatedVelocity>(lattice, superGeometry, 2);
  boundary::set<boundary::InterpolatedVelocity>(lattice, superGeometry, 3);

  lattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( UnitConverter<T, DESCRIPTOR> const& converter,
                        SuperLattice<T,DESCRIPTOR>& lattice, SuperGeometry<T,3>& superGeometry, int iT )
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  if ( iT==0 ) {
    AnalyticalConst3D<T,T> rhoF( T( 1 ) );
    AnalyticalConst3D<T,T> uF( T( 0 ), T( 0 ), T( 0 ) );

    auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3});
    lattice.iniEquilibrium( bulkIndicator, rhoF, uF );
    lattice.defineRhoU( bulkIndicator, rhoF, uF );

    clout << converter.getCharLatticeVelocity() << std::endl;
    AnalyticalConst3D<T,T> uTop( converter.getCharLatticeVelocity(), T( 0 ), T( 0 ) );
    lattice.defineU( superGeometry,3,uTop );

    // Make the lattice ready for simulation
    lattice.initialize();
  }
}

void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR> const& converter, SuperGeometry<T,3>& superGeometry,
                 int iT, util::Timer<T>& timer, bool converged,
                 std::vector<int>& iTList )
{

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "cavity3d", 1, false, false );

  if ( iT==0 ) {
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    SuperLatticeDiscreteNormal3D<T, DESCRIPTOR> discreteNormal( sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3}) );
    SuperLatticeDiscreteNormalType3D<T, DESCRIPTOR> discreteNormalType( sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3}) );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.write( discreteNormal );
    vtmWriter.write( discreteNormalType );
    vtmWriter.createMasterFile();
  }

  // Writes the VTK
  if ( (iT%converter.getLatticeTime(vtkSave) == 0) || converged  ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );

    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );

    vtmWriter.write( iT );

    timer.update( iT );
    timer.printStep( 2 );
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    int rank = 0;
    #ifdef PARALLEL_MODE_MPI
      rank = singleton::mpi().getRank();
    #endif
    // Only rank 0 writes to the file
    if (rank == 0) {
      iTList.push_back(iT);
    }
  }
}


void simulateCavity3d( T physVelocity, int N, int sample )
{

  // === 1st Step: Initialization ===

//   initialize( &argc, &argv );
  OstreamManager clout( std::cout,"cylinder2d" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   0.509, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1.0,   // charPhysLength: reference length of simulation geometry
    (T)   physVelocity,   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.001, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("cavity3d");


  // === 2nd Step: Prepare Geometry ===

  // Instantiation of a unit cube by an indicator
  Vector<T,3> origin( T( 0 ) );
  Vector<T,3> extend( converter.getCharPhysLength() + 0.5 * converter.getPhysDeltaX() );
  IndicatorCuboid3D<T> cube( extend, origin );

  // Instantiation of a cuboid geometry with weights
  int noCuboids = singleton::mpi().getSize();
  CuboidDecomposition3D<T> cuboidDecomposition( cube, converter.getPhysDeltaX(), noCuboids );

  // Instantiation of a load balancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

  // Instantiation of a super geometry
  SuperGeometry<T,3> superGeometry( cuboidDecomposition, loadBalancer );

  prepareGeometry( converter, cube, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and setBoundaryConditions
  prepareLattice( converter, sLattice, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), util::pow<int>(converter.getResolution(),3) );
  int interval = converter.getLatticeTime( 1 /*config["Application"]["ConvergenceCheck"]["interval"].get<T>()*/ );
  T epsilon = 1e-3;
  util::ValueTracer<T> converge( converter.getLatticeTime(interval), epsilon );
  timer.start();
  std::vector<int> iTList;

  for ( std::size_t iT = 0; iT <= converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      getResults( sLattice, converter, superGeometry, iT, timer, converge.hasConverged(),iTList );
      converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
      break;
    }
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, superGeometry, iT );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, superGeometry, iT, timer, converge.hasConverged(),iTList );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }

  // Write iteration indices to a file (only on the main processor)
  if (singleton::mpi().isMainProcessor()) {
    std::ofstream outFile(singleton::directories().getLogOutDir() + "iteration_log.txt");
    if (outFile.is_open()) {
      for (const auto& iter : iTList) {
      outFile << iter << "\n";
      }
      outFile.close();
    } else {
      clout << "Error: Unable to open file for writing!" << std::endl;
    }
  }

  timer.stop();
  timer.printSummary();
}

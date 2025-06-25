/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mathias J. Krause, Jonas Fietz,
 *  Jonas Latt, Jonas Kratzke, Mingliang Zhong, Stephan Simonis
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

#include "olb2D.h"
#include "olb2D.hh"


using namespace olb;
using namespace olb::uq;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<>;
using BulkDynamics = ConstRhoBGKdynamics<T,DESCRIPTOR>;

extern const T vtkSave;  ///< Interval for writing VTK output (in physical seconds)
extern const T maxPhysT; ///< Maximum physical time for simulation (seconds)

void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1} );
  superGeometry.clean();

  T eps = converter.getConversionFactorLength();
  Vector<T,2> extend( T( 1 ) + 2*eps, 2*eps );
  Vector<T,2> origin( T() - eps, T( 1 ) - eps );
  IndicatorCuboid2D<T> lid( extend, origin );
  // Set material number for lid
  superGeometry.rename( 2,3,1,lid );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);

  // Material=2,3 -->bulk dynamics, velocity boundary
  boundary::set<boundary::InterpolatedVelocity<T,DESCRIPTOR,BulkDynamics>>(sLattice, superGeometry, 2);
  boundary::set<boundary::InterpolatedVelocity<T,DESCRIPTOR,BulkDynamics>>(sLattice, superGeometry, 3);

  sLattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice ... OK" << std::endl;
}


void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T, DESCRIPTOR>& sLattice,
                        int iT, SuperGeometry<T,2>& superGeometry )
{

  if ( iT==0 ) {
    // set initial values: v = [0,0]
    AnalyticalConst2D<T,T> rhoF( 1 );
    std::vector<T> velocity( 2,T() );
    AnalyticalConst2D<T,T> uF( velocity );

    auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3});
    sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );
    sLattice.defineRhoU( bulkIndicator, rhoF, uF );

    // set non-zero velocity for upper boundary cells
    velocity[0] = converter.getCharLatticeVelocity();
    AnalyticalConst2D<T,T> u( velocity );
    sLattice.defineU( superGeometry, 3, u );

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}

void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT, util::Timer<T> timer,
                 SuperGeometry<T,2>& superGeometry, bool converged,
                 std::vector<int>& iTList )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( "cavity2d", 1, false );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    SuperLatticeDiscreteNormal2D<T, DESCRIPTOR> discreteNormal( sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3}) );
    SuperLatticeDiscreteNormalType2D<T, DESCRIPTOR> discreteNormalType( sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3}) );

    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.write( discreteNormal );
    vtmWriter.write( discreteNormalType );
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if ( ( iT%converter.getLatticeTime( vtkSave )==0 && iT>0 ) || converged ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );

    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );

    timer.update(iT);
    timer.printStep(2);
    sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );

     // Only the main MPI rank logs iteration steps
     int rank = 0;
     #ifdef PARALLEL_MODE_MPI
       rank = singleton::mpi().getRank();
     #endif
     if (rank == 0) {
       iTList.push_back(static_cast<int>(iT));
     }
  }
}

void simulateCavity2d(T physVelocity, int N, int idx)
{

  // === 1st Step: Initialization ===
  OstreamManager clout( std::cout,"simulateCavity2D" );
  clout << "Start sample " << idx << std::endl;

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    N,                 // resolution
    0.5384, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    1.0,             // charPhysLength: reference length of simulation geometry
    physVelocity,               // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    0.001,                 // physViscosity: physical kinematic viscosity in __m^2 / s__
    1.0                 // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("cavity2d");

  // === 2nd Step: Prepare Geometry ===
  Vector<T,2> extend( 1,1 );
  Vector<T,2> origin( 0,0 );
  IndicatorCuboid2D<T> cuboid( extend, origin );

#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, 1.0 / N, singleton::mpi().getSize() );
#else
  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, 1.0 / N, 1 );
#endif

  cuboidDecomposition.print();

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,2> superGeometry( cuboidDecomposition, loadBalancer );
  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  prepareLattice( converter, sLattice, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  int interval = converter.getLatticeTime( 1 /*config["Application"]["ConvergenceCheck"]["interval"].get<T>()*/ );
  T epsilon = 1e-3;
  util::ValueTracer<T> converge( interval, epsilon );
  std::vector<int> iTList;

  timer.start();
  for ( std::size_t iT=0; iT <= converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      getResults( sLattice, converter, iT, timer, superGeometry, converge.hasConverged(), iTList );
      break;
    }
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, timer, superGeometry, converge.hasConverged(), iTList );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }

  // Write iteration indices to a file (only on the main processor)
  if (singleton::mpi().isMainProcessor()) {
    std::ofstream outFile(singleton::directories().getLogOutDir() + "iteration_log.txt");
    if (outFile.is_open()) {
      for (auto step : iTList) {
        outFile << step << "\n";
      }
      outFile.close();
    } else {
      clout << "Error: Unable to open iteration_log.txt for writing!" << std::endl;
    }
  }

  timer.stop();
  timer.printSummary();
}

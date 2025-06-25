/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Jonas Latt, Mathias J. Krause,
 *     Mingliang Zhong, Stephan Simonis
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

#ifndef CYLINDER_2D_H
#define CYLINDER_2D_H

#include <olb.h>

using namespace olb;
using namespace olb::uq;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<>;

#define BOUZIDI

// -------------------------------------------------------------------------
// Global parameters for the cylinder2d simulation
// -------------------------------------------------------------------------
extern const T maxPhysT;        // max. simulation time in s, SI unit
extern const T physInterval;   // interval for the convergence check in s

struct GeometryParameters {
  T L;
  T lengthX;
  T lengthY;
  T centerCylinderX;
  T centerCylinderY;
  T radiusCylinder;

  GeometryParameters(int N) {
    L = 0.1 / N;
    lengthX = 2.2;
    lengthY = 0.41 + L;
    centerCylinderX = 0.2;
    centerCylinderY = 0.2 + L / 2.0;
    radiusCylinder = 0.05;
  }
};

// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T, DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry,
                      std::shared_ptr<IndicatorF2D<T>> circle,
                      GeometryParameters geomParams)
{
  Vector<T,2> extend( geomParams.lengthX, geomParams.lengthY );
  Vector<T,2> origin;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1} );

  // Set material number for inflow
  extend[0] = 2.*geomParams.L;
  origin[0] = -geomParams.L;
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,3,1,inflow );
  // Set material number for outflow
  origin[0] = geomParams.lengthX - geomParams.L;
  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,4,1,outflow );
  // Set material number for cylinder
  superGeometry.rename( 1,5, circle );

  // Removes all not needed boundary voxels outside the surface
  const bool verbose = false;
  superGeometry.clean(verbose);
  superGeometry.innerClean(verbose);
  superGeometry.checkForErrors(verbose);
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T, DESCRIPTOR> const& converter,
                     SuperGeometry<T,2>& superGeometry,
                     std::shared_ptr<IndicatorF2D<T>> circle)
{
  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  // Setting of the boundary conditions
  //if boundary conditions are chosen to be interpolatedy, 3);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  // Material=5 -->bouzidi / bounce back
  #ifdef BOUZIDI
    setBouzidiBoundary(sLattice, superGeometry, 5, *circle);
  #else
    boundary::set<boundary::BounceBack>(sLattice, superGeometry, 5);
  #endif

  // Initial conditions
  AnalyticalConst2D<T,T> rhoF( 1 );
  std::vector<T> velocity( 2,T( 0 ) );
  AnalyticalConst2D<T,T> uF( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );
  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,2>& superGeometry,
                        T L )
{
  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( 6.4 );
  int iTupdate = 5;

  if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getCharLatticeVelocity()*3./2.*frac[0];
    T distance2Wall = L/2.;
    Poiseuille2D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall );

    sLattice.defineU( superGeometry, 3, poiseuilleU );

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
T getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,2>& superGeometry,
                 bool exportResults,
                 std::vector<int>& iTList,
                 GeometryParameters geomParams)
{
  T dragCoefficient = 0.;

  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );

  int vtkIter  = converter.getLatticeTime( 1.0 );
  int statIter = converter.getLatticeTime( 0.1 );

  if (exportResults)
  {
    SuperVTMwriter2D<T> vtmWriter( "cylinder2d", 1, false );
    SuperLatticeRefinementMetricKnudsen2D<T, DESCRIPTOR> quality( sLattice, converter);
    SuperRoundingF2D<T> roundedQuality ( quality, RoundingMode::NearestInteger);
    SuperDiscretizationF2D<T> discretization ( roundedQuality, 0., 2. );

    // vtmWriter.addFunctor( quality );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );

    T point[2] = {};
    point[0] = geomParams.centerCylinderX + 3*geomParams.radiusCylinder;
    point[1] = geomParams.centerCylinderY;
    AnalyticalFfromSuperF2D<T> intpolateP( pressure, true );
    T p;
    intpolateP( &p,point );

    if ( iT==0 ) {
      // Writes the geometry, cuboid no. and rank no. as vti file for visualization
      SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
      SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
      SuperLatticeDiscreteNormal2D<T, DESCRIPTOR> discreteNormal(
        sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3}) );
      SuperLatticeDiscreteNormalType2D<T, DESCRIPTOR> discreteNormalType(
        sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3, 4, 5}) );
      vtmWriter.write( cuboid );
      vtmWriter.write( rank );
      vtmWriter.write( discreteNormal );
      vtmWriter.write( discreteNormalType );
      vtmWriter.createMasterFile();
    }

    // Writes the vtk files
    if ( iT%vtkIter == 0 ) {
      vtmWriter.write( iT );
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

  if ( iT%statIter == 0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF2D<T> intpolatePressure( pressure, true );
    SuperLatticePhysDrag2D<T,DESCRIPTOR> drag( sLattice, superGeometry, 5, converter );
    int input[3] = {};
    T _drag[drag.getTargetDim()];
    drag( _drag,input );
    dragCoefficient = _drag[0];
  }

  return dragCoefficient;
}

T simulateCylinder( int N, T u0, bool exportResults )
{
  // === 1st Step: Initialization ===
  OstreamManager clout( std::cout,"cylinder2d" );
  GeometryParameters geomParams(N);

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int   {N},                        // resolution: number of voxels per charPhysL
    (T)   0.56,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   2.0*geomParams.radiusCylinder,       // charPhysLength: reference length of simulation geometry
    (T)   u0,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.001, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                       // physDensity: physical density in __kg / m^3__
  );

  // === 2rd Step: Prepare Geometry ===
  Vector<T,2> extend( geomParams.lengthX, geomParams.lengthY );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
  #ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = singleton::mpi().getSize();
  #else
    const int noOfCuboids = 7;
  #endif

  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, geomParams.L, noOfCuboids );

  // // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

  // // Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry( cuboidDecomposition, loadBalancer );

  Vector<T,2> center( geomParams.centerCylinderX, geomParams.centerCylinderY );
  std::shared_ptr<IndicatorF2D<T>> circle = std::make_shared<IndicatorCircle2D<T>>( center, geomParams.radiusCylinder );

  prepareGeometry( converter, superGeometry, circle, geomParams );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and set boundaryConditions
  prepareLattice( sLattice, converter, superGeometry, circle );

  // === 4th Step: Main Loop ===
  T drag = 0;

  std::vector<int> iTList;

  for ( std::size_t iT = 0; iT <= converter.getLatticeTime( maxPhysT ); ++iT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry, geomParams.L );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    if (iT % converter.getLatticeTime( .1 ) == 0) {
      drag = getResults( sLattice, converter, iT, superGeometry, exportResults, iTList, geomParams );
    }

  }

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

  return drag;
}

#endif

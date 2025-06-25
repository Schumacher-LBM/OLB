/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net
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

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>
#include <fstream>
#include <sstream>

using namespace olb;
using namespace olb::descriptors;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<descriptors::FORCE, FreeSurface::MASS, FreeSurface::EPSILON, FreeSurface::CELL_TYPE, FreeSurface::CELL_FLAGS, FreeSurface::TEMP_MASS_EXCHANGE, FreeSurface::PREVIOUS_VELOCITY>;

struct FreeSurfaceAppHelper {
  std::array<T,2> area;
  std::array<T,2> gravity_force = {{0., -9.81}};

  T char_phys_length = 1.;
  T char_phys_vel = 0.1;
  bool has_surface_tension = true;
  T surface_tension_coefficient = 0.0661;
};

template<typename T, typename DESCRIPTOR>
class FreeSurfaceBreakingDam2D final : public AnalyticalF2D<T,T> {
private:
  T lattice_size;
  std::array<T,3> cell_values;
  std::array<T, 2> area;
public:
  FreeSurfaceBreakingDam2D(T lattice_size_, const std::array<T,3>& cell_vals, const std::array<T,2>& area_):AnalyticalF2D<T,T>{1}, lattice_size{lattice_size_}, cell_values{cell_vals}, area{area_}{}

  bool operator()(T output[], const T x[]) override {
    output[0] = cell_values[0];

    if(x[1] <= area[1] * 0.6 && x[0] <= area[0] * 0.5){
      output[0] = cell_values[2];
    }else if(x[1] - lattice_size * 1.1 <= area[1] * 0.6 && x[0] - lattice_size * 1.1 <= area[0] * 0.5){
      output[0] = cell_values[1];
    }

    return true;
  }
};

void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1} );

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareBreakingDam(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry, T lattice_size, const FreeSurfaceAppHelper& helper)
{
  AnalyticalConst2D<T,T> zero( 0. );
  AnalyticalConst2D<T,T> one( 1. );
  AnalyticalConst2D<T,T> two( 2. );
  AnalyticalConst2D<T,T> four( 4. );
  FreeSurfaceBreakingDam2D<T,DESCRIPTOR> cells_analytical{ lattice_size, {0., 1., 2.}, helper.area};
  FreeSurfaceBreakingDam2D<T,DESCRIPTOR> mass_analytical{ lattice_size, {0., 0.5, 1.}, helper.area};

  AnalyticalConst2D<T,T> force_zero{0., 0.};

  for (int i: {0,1,2}) {
    sLattice.defineField<FreeSurface::MASS>(superGeometry, i, zero);
    sLattice.defineField<FreeSurface::EPSILON>(superGeometry, i, zero);
    sLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, i, zero);
    sLattice.defineField<FreeSurface::CELL_FLAGS>(superGeometry, i, zero);
    sLattice.defineField<descriptors::FORCE>(superGeometry, i, force_zero);
  }

  sLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, 1, cells_analytical);

  sLattice.defineField<FreeSurface::MASS>(superGeometry, 1, mass_analytical);
  sLattice.defineField<FreeSurface::EPSILON>(superGeometry, 1, mass_analytical);

  for (int i: {0,2}) {
    //sLattice.defineField<FreeSurface::MASS>(superGeometry, i, one);
    sLattice.defineField<FreeSurface::EPSILON>(superGeometry, i, one);
    sLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, i, four);
  }

  T force_factor = 1./ converter.getConversionFactorForce() * converter.getConversionFactorMass();
  AnalyticalConst2D<T,T> force_a{helper.gravity_force[0] * force_factor, helper.gravity_force[1] * force_factor};
  sLattice.defineField<descriptors::FORCE>(superGeometry.getMaterialIndicator({1}), force_a);

}

void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry, T lattice_size, const FreeSurfaceAppHelper& helper) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics<T,DESCRIPTOR>>(superGeometry, 0);
  // Material=1 -->bulk dynamics
  sLattice.defineDynamics<SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>>( superGeometry, 1);
  // Material=2 -->no-slip boundary
  sLattice.defineDynamics<BounceBack<T,DESCRIPTOR>>( superGeometry, 2);
  //setSlipBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 2);

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  sLattice.setParameter<collision::LES::SMAGORINSKY>(T(0.2));

  prepareBreakingDam(converter, sLattice, superGeometry, lattice_size, helper);

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,2>& sGeometry, T lattice_length, UnitConverter<T,DESCRIPTOR> const& converter){
  OstreamManager clout( std::cout,"setInitialValues" );

    AnalyticalConst2D<T,T> u{0., 0.};
    AnalyticalConst2D<T,T> one(1.);

    sLattice.defineRhoU( sGeometry.getMaterialIndicator({0,1,2}), one, u );
    for (int i: {0,1,2}) {
      sLattice.iniEquilibrium( sGeometry, i, one, u );
    }

  // Set up free surface communicator stages
  FreeSurface::initialize(sLattice);
  // Make the lattice ready for simulation
  sLattice.initialize();
}

void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer)
{
  const int vtmIter  = 100;
  const int statIter = 100;

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperVTMwriter2D<T> vtmWriter( "breakingDam2d" );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes output in ParaView
  if ( iT%vtmIter==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperVTMwriter2D<T> vtmWriter( "breakingDam2d" );
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );
    SuperLatticeExternalScalarField2D<T, DESCRIPTOR, FreeSurface::EPSILON> epsilon( sLattice );
    SuperLatticeExternalScalarField2D<T, DESCRIPTOR, FreeSurface::CELL_TYPE> cells( sLattice );
    SuperLatticeExternalScalarField2D<T, DESCRIPTOR, FreeSurface::MASS> mass( sLattice );
    epsilon.getName() = "epsilon";
    cells.getName() = "cell_type";
    mass.getName() = "mass";
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.addFunctor( epsilon );
    vtmWriter.addFunctor( cells );
    vtmWriter.addFunctor( mass );

    vtmWriter.write(iT);
  }

  // Writes output on the console
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }
}

namespace {
FreeSurfaceAppHelper free_surface_config {
  {7.31,0.42}, // area
  {0.,-9.81}, // gravity_force
  7.31, // char_phys_length
  0.1,  // char_phys_vel
  true, // has_surface_tension
  0.05 // surface_tension_coefficient
};

class FreeSurfaceConfig {
public:
  T viscosity = 1e-4;
  T density = 1e3;
  T physTime = 30.;
  T latticeRelaxationTime = .501;
  int N = 500;

  // Anti jitter value
  T transitionThreshold = 1e-3;
  // When to remove lonely cells
  T lonelyThreshold = 1.0;
};

}

int main(int argc, char **argv)
{

  initialize(&argc, &argv, false, false);

  FreeSurfaceConfig c;
  OstreamManager clerr( std::cerr, "main" );
  OstreamManager clout( std::cout, "main" );

  singleton::directories().setOutputDir("./tmp/");

  FreeSurfaceAppHelper& helper = free_surface_config;

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {c.N},     // resolution: number of voxels per charPhysL
    (T)   c.latticeRelaxationTime,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   helper.char_phys_length,     // charPhysLength: reference length of simulation geometry
    (T)   helper.char_phys_vel,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   c.viscosity, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   c.density     // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("free surface");

  T lattice_size = helper.char_phys_length / c.N;

  T force_conversion_factor = 1./converter.getConversionFactorForce()*converter.getConversionFactorMass();

  // Convert kg / s^2
  // Basically it is multiplied with s^2 / kg = s^2 * m^3 / (kg * m^2 * m) = 1. / (velocity_factor^2 * density * length_factor)
  T surface_tension_coefficient_factor = std::pow(converter.getConversionFactorTime(),2)/ (c.density * std::pow(converter.getPhysDeltaX(),3));

  clout<<"Surface: "<<surface_tension_coefficient_factor * helper.surface_tension_coefficient<<std::endl;
  clout<<"Lattice Size: "<<converter.getPhysDeltaX()<<std::endl;

  // === 2nd Step: Prepare Geometry ===
  Vector<T,2> extend( helper.area[0], helper.area[1] );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

  // Instantiation of a cuboidDecomposition with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, converter.getPhysDeltaX(), noOfCuboids );

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,2> superGeometry( cuboidDecomposition, loadBalancer, 2 );

  prepareGeometry( converter, superGeometry );

  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  clout<<"Overlap: "<<sLattice.getOverlap()<<std::endl;

  prepareLattice( converter, sLattice, superGeometry, lattice_size, helper);

  FreeSurface2DSetup<T,DESCRIPTOR> free_surface_setup{sLattice};

  free_surface_setup.addPostProcessor();

  // Set variables from freeSurfaceHelpers.h
  sLattice.setParameter<FreeSurface::DROP_ISOLATED_CELLS>(true);
  sLattice.setParameter<FreeSurface::TRANSITION>(c.transitionThreshold);
  sLattice.setParameter<FreeSurface::LONELY_THRESHOLD>(c.lonelyThreshold);
  sLattice.setParameter<FreeSurface::HAS_SURFACE_TENSION>(helper.has_surface_tension);
  sLattice.setParameter<FreeSurface::SURFACE_TENSION_PARAMETER>(surface_tension_coefficient_factor * helper.surface_tension_coefficient);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( c.physTime ), superGeometry.getStatistics().getNvoxel() );
  timer.start();
  setInitialValues(sLattice, superGeometry, lattice_size, converter);

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( c.physTime ); ++iT ) {
    getResults( sLattice, converter, iT, superGeometry, timer );
    sLattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();

  return 0;
}

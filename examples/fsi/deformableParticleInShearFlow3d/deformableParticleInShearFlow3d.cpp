/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Shota Ito, Timm Krüger
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

// Benchmark validation case for a 3D derformable particle submerged
// in a shear flow. (cf. https://doi.org/10.1016/j.camwa.2010.03.057)
// Identical simulation parameters are chosen as in the publication.

#undef PARALLEL_MODE_MPI

#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace membrane;

using T = double;
using DESCRIPTOR = descriptors::D3Q19<descriptors::LOCATION,descriptors::FORCE>;
using DESCRIPTOR_MEM = descriptors::D3Q0<fields::PHYS_R,fields::membrane::VELOCITY,
  fields::membrane::FORCE>;
using ForcedBulkDynamics = ForcedBGKdynamics<T,DESCRIPTOR>;

const int N = 50;                                 // resolution of char. length
const T Re = 0.02;                                // Reynolds number
const T G = 0.01;                                 // dimensionless shear rate
const T H = N * 1.0;                              // wall distance
const T r = H / 10.;                              // radius of sphere
const T relaxationTime = 1.0;                     // relaxation time
const T charViscosity = 1./6.;                    // kinematic viscosity
const T shearRate = Re * charViscosity / (r*r);   // shear rate
const T maxT = 120. / (shearRate / G);            // reduced time
const T charVelocity = shearRate * H / 2.;        // characteristic velocity
const T charDensity = 1.;                         // density
const T k_s = shearRate * charViscosity           // stiffness coeff membrane
  * charDensity * r / G;
const T physX = H / N;
const T physT = physX;
const std::vector<T> posCenter {H/2., H/2., H/2.};

void prepareGeometry( UnitConverter<T, DESCRIPTOR> const& converter, SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2);
  Vector<T,3> center( T(H/2.0), T(H/2.0), T(H/2.0) );
  T eps = converter.getConversionFactorLength();
  IndicatorCuboid3D<T> bulk(H, H, H - eps, center);
  superGeometry.rename(2, 1, bulk);
  center[2] = 0.0;
  IndicatorCuboid3D<T> bottomWall(H, H, eps, center);
  superGeometry.rename(2, 3, bottomWall);

  superGeometry.clean();
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
  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3});

  lattice.defineDynamics<ForcedBulkDynamics>(bulkIndicator);

  boundary::set<boundary::InterpolatedVelocity>(lattice, superGeometry, 2);
  boundary::set<boundary::InterpolatedVelocity>(lattice, superGeometry, 3);

  lattice.setParameter<descriptors::OMEGA>(omega);
  AnalyticalConst3D<T,T> rhoF( T( 1 ) );
  AnalyticalConst3D<T,T> uF( T( 0 ), T( 0 ), T( 0 ) );

  lattice.iniEquilibrium( bulkIndicator, rhoF, uF );
  lattice.defineRhoU( bulkIndicator, rhoF, uF );

  AnalyticalConst3D<T,T> uW( converter.getCharLatticeVelocity(), T( 0 ), T( 0 ) );
  AnalyticalConst3D<T,T> nuW( -1. * converter.getCharLatticeVelocity(), T( 0 ), T( 0 ) );
  lattice.defineU(superGeometry, 2, uW);
  lattice.defineU(superGeometry, 3, nuW);

  lattice.initialize();
  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR> const& converter,
                 SuperGeometry<T,3>& superGeometry,
                 MembraneParticleSystem3D<T,DESCRIPTOR_MEM>& sMembrane,
                 int iT, util::Timer<T>& timer )
{
  SuperVTMwriter3D<T> vtmWriter( "shearFlow3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  SuperLatticePhysField3D<T, DESCRIPTOR, descriptors::FORCE> force( sLattice, 1.0, "force" );

  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( force );

  const T logT   = maxT / 15.;
  const T saveT  = maxT / 15.;

  if ( iT==0 ) {
    vtmWriter.createMasterFile();
    vtmWriter.write( iT );

    sMembrane.updateMomenta();
    writeData( iT, sMembrane );
  }

  // Get statistics
  if (iT%converter.getLatticeTime(logT) == 0 && iT > 0) {
    timer.update( iT );
    timer.printStep( 2 );
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  // Writes the VTK
  if (iT%converter.getLatticeTime(saveT) == 0 && iT > 0) {
    vtmWriter.write( iT );
    sMembrane.updateMomenta();
    writeData( iT, sMembrane );
  }
}

int main( int argc, char **argv )
{
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  CLIreader args(argc, argv);

  OstreamManager clout( std::cout,"main" );
  UnitConverter<T, DESCRIPTOR> const converter(
    (T)   physX,
    (T)   physT,
    (T)   H,
    (T)   charVelocity,
    (T)   charViscosity,
    (T)   charDensity
  );
  converter.print();
  converter.write("deformableParticleInShearFlow3d");

  Vector<T,3> center( T(H/2.0), T(H/2.0), T(H/2.0) );
  IndicatorCuboid3D<T> cube( H, H, H, center );

  CuboidDecomposition3D<T> cuboidDecomposition( cube, converter.getConversionFactorLength(), singleton::mpi().getSize() );
  cuboidDecomposition.setPeriodicity({true, true, false});
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

  SuperGeometry<T,3> superGeometry( cuboidDecomposition, loadBalancer );
  prepareGeometry( converter, superGeometry );

  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );
  prepareLattice( converter, sLattice, superGeometry );

  util::Timer<T> timer( converter.getLatticeTime( maxT ), util::pow<int>(converter.getResolution(),3) );

  // Initialize membrane and read parameter files
  MembraneParticleSystem3D<T,DESCRIPTOR_MEM> sMembrane;
  sMembrane.initializeMeshesFromXML("parametersMeshes.xml");
  sMembrane.initializeParticlesFromXML("parametersPositions.xml");
  sMembrane.initializeData(loadBalancer);

  // IBM operator
  ibm::StencilWidth<4> stencil_width;
  SuperImmersedBoundaryCoupling3D ibmO(stencil_width, names::NavierStokes{}, sLattice,
    names::Points{}, sMembrane.getData());
  ibmO.setPhysDeltaX(converter.getPhysDeltaX());

  timer.start();
  for ( std::size_t iT = 0; iT <= converter.getLatticeTime( maxT ); ++iT ) {
    // Update Lagrangian particles.
    sMembrane.updateForcesFromDisplacements();

    // Communicate displacement and force
    ibmO.execute();

    // Update flow field.
    sLattice.collideAndStream();

    getResults( sLattice, converter, superGeometry, sMembrane, iT, timer );
  }
  timer.stop();
  timer.printSummary();
}

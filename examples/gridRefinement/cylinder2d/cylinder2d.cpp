/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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
using namespace olb::descriptors;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<>;

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


void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,2>& sGeometry,
                     std::shared_ptr<IndicatorF2D<T>> circle,
                     GeometryParameters geomParams)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend(geomParams.lengthX, geomParams.lengthY);
  Vector<T,2> origin;

  sGeometry.rename(0,2);
  sGeometry.rename(2,1,{1,1});

  // Set material number for inflow
  extend[0] = 2.*geomParams.L;
  origin[0] = -geomParams.L;
  IndicatorCuboid2D<T> inflow( extend, origin );
  sGeometry.rename( 2,3,1,inflow );

  // Set material number for outflow
  origin[0] = geomParams.lengthX - geomParams.L;
  IndicatorCuboid2D<T> outflow( extend, origin );
  sGeometry.rename( 2,4,1,outflow );
  // Set material number for cylinder
  sGeometry.rename( 1,5, circle );

  // Removes all not needed boundary voxels outside the surface
  sGeometry.clean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareGeometryFine(const UnitConverter<T,DESCRIPTOR>& converter,
                         SuperGeometry<T,2>& sGeometry,
                         std::shared_ptr<IndicatorF2D<T>> circle,
                         GeometryParameters geomParams)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend(geomParams.lengthX, geomParams.lengthY);
  Vector<T,2> origin;

  sGeometry.rename(0,1);
  sGeometry.rename(1,5, circle);

  sGeometry.clean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    const UnitConverter<T,DESCRIPTOR>& converter,
                    SuperGeometry<T,2>& sGeometry,
                    std::shared_ptr<IndicatorF2D<T>> circle)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = sGeometry.getMaterialIndicator({1});
  sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 2);

  // Setting of the boundary conditions
  // if boundary conditions are chosen to be local
  //boundary::set<boundary::LocalVelocity>(sLattice, sGeometry, 3);
  //boundary::set<boundary::LocalPressure>(sLattice, sGeometry, 4);

  //if boundary conditions are chosen to be interpolatedy, 3);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, sGeometry, 4);

  // if boundary conditions are chosen to be Zou He type
  //boundary::set<boundary::ZouHeVelocity>(sLattice, sGeometry, 3);
  //boundary::set<boundary::ZouHePressure>(sLattice, sGeometry, 4);

  // Material=5 -->bouzidi / bounce back
  #ifdef BOUZIDI
  setBouzidiBoundary(sLattice, sGeometry, 5, *circle);
  #else
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 5);
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

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues(SuperLattice<T,DESCRIPTOR>& sLattice,
                       const UnitConverter<T,DESCRIPTOR>& converter,
                       std::size_t iT,
                       SuperGeometry<T,2>& sGeometry)
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  std::size_t iTmaxStart = converter.getLatticeTime(6.4);
  std::size_t iTupdate = 5;

  if (iT % iTupdate == 0 && iT<= iTmaxStart) {
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getCharLatticeVelocity()*3./2.*frac[0];
    T distance2Wall = converter.getPhysDeltaX()/2.;
    Poiseuille2D<T> poiseuilleU( sGeometry, 3, maxVelocity, distance2Wall );

    sLattice.defineU( sGeometry, 3, poiseuilleU );

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void writeResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                  const UnitConverter<T,DESCRIPTOR>& converter,
                  std::size_t iT,
                  SuperGeometry<T,2>& sGeometry,
                  std::string name)
{
  OstreamManager clout(std::cout, name);

  if (iT == 0) {
    SuperVTMwriter2D<T> vtmWriter(name, 0);
    if (name == "level0") {
      SuperLatticeRank2D rank( sLattice );
      vtmWriter.write(rank);
    }
    SuperGeometryF<T,DESCRIPTOR::d> geometryF(sGeometry);
    geometryF.getName() = name + "_geometry";
    vtmWriter.write(geometryF);
    SuperLatticeCuboid2D<T,DESCRIPTOR> cuboidF(sLattice);
    cuboidF.getName() = name + "_cuboid";
    vtmWriter.write(cuboidF);
    vtmWriter.createMasterFile();
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  //sLattice.scheduleBackgroundOutputVTK([&,name,iT](auto task) {
    SuperVTMwriter2D<T> vtmWriter(name);
    SuperLatticePhysVelocity2D velocityF(sLattice, converter);
    SuperLatticePhysPressure2D pressureF(sLattice, converter);
    SuperLatticeRefinementMetricKnudsen2D qualityF(sLattice, converter);
    vtmWriter.addFunctor(qualityF);
    vtmWriter.addFunctor(velocityF);
    vtmWriter.addFunctor(pressureF);
    vtmWriter.write(iT);
    //task(vtmWriter, iT);
  //});
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  CLIreader args(argc, argv);
  const int N = args.getValueOrFallback<int>("--resolution", 11);
  const int Re = args.getValueOrFallback<int>("--reynolds", 100);
  const int maxPhysT = args.getValueOrFallback<int>("--maxPhysT", 16);

  GeometryParameters geomParams(N);

  const UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converterLevel0(
    int {N},                                   // resolution: number of voxels per charPhysL
    (T)   0.51,                                // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   2.0*geomParams.radiusCylinder,       // charPhysLength: reference length of simulation geometry
    (T)   0.2,                                 // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.2*2.*geomParams.radiusCylinder/Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                                  // physDensity: physical density in __kg / m^3__
  );
  converterLevel0.print();

  Vector<T,2> extend( geomParams.lengthX, geomParams.lengthY );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid(extend, origin);

  CuboidDecomposition2D<T> cuboidDecompositionLevel0(cuboid, geomParams.L);
  cuboidDecompositionLevel0.splitFractional(0, 0, {0.025,0.5,0.475});
  cuboidDecompositionLevel0.splitFractional(1, 1, {0.15,0.7,0.15});
  cuboidDecompositionLevel0.splitFractional(3, 0, {0.05,0.5,0.45});
  cuboidDecompositionLevel0.splitFractional(5, 1, {0.15,0.7,0.15});

  Vector<T,2> center(geomParams.centerCylinderX, geomParams.centerCylinderY);
  std::shared_ptr<IndicatorF2D<T>> circle = std::make_shared<IndicatorCircle2D<T>>(center, geomParams.radiusCylinder);

  auto cuboidDecompositionLevel1 = cuboidDecompositionLevel0;
  IndicatorCuboid2D<T> refineCuboidI({1.0, 0.25}, {0.08,0.08});
  cuboidDecompositionLevel1.remove(refineCuboidI);
  cuboidDecompositionLevel1.refine(2);
  cuboidDecompositionLevel1.print();

  auto cuboidDecompositionLevel2 = cuboidDecompositionLevel1;
  cuboidDecompositionLevel2.remove(*circle);
  cuboidDecompositionLevel2.refine(2);
  cuboidDecompositionLevel2.print();

  // Adjust weights for balancing
  for (int iC=0; iC < cuboidDecompositionLevel0.size(); ++iC) {
    auto& cuboid = cuboidDecompositionLevel0.get(iC);
    cuboid.setWeight(cuboid.getLatticeVolume());

    auto origin = cuboidDecompositionLevel0.get(iC).getOrigin();

    for (int jC=0; jC < cuboidDecompositionLevel1.size(); ++jC) {
      if (cuboidDecompositionLevel1.get(jC).getOrigin() == origin) {
        cuboid.setWeight(cuboid.getWeight() + 2*cuboidDecompositionLevel1.get(jC).getLatticeVolume());
      }
    }
    for (int jC=0; jC < cuboidDecompositionLevel2.size(); ++jC) {
      if (cuboidDecompositionLevel2.get(jC).getOrigin() == origin) {
        cuboid.setWeight(cuboid.getWeight() + 4*cuboidDecompositionLevel2.get(jC).getLatticeVolume());
      }
    }
  }

  HeuristicLoadBalancer<T> loadBalancerLevel0(cuboidDecompositionLevel0);
  SuperGeometry<T,2> sGeometryLevel0(cuboidDecompositionLevel0, loadBalancerLevel0);
  prepareGeometry(converterLevel0, sGeometryLevel0, circle, geomParams);
  SuperLattice<T,DESCRIPTOR> sLatticeLevel0(cuboidDecompositionLevel0,
                                            loadBalancerLevel0,
                                            3,
                                            converterLevel0);
  prepareLattice(sLatticeLevel0, converterLevel0, sGeometryLevel0, circle);

  RefinedLoadBalancer<T,2> loadBalancerLevel1(cuboidDecompositionLevel0,
                                              loadBalancerLevel0,
                                              cuboidDecompositionLevel1);
  auto converterLevel1 = convectivelyRefineUnitConverter(converterLevel0, 2);
  converterLevel1.print();
  SuperGeometry<T,2> sGeometryLevel1(cuboidDecompositionLevel1, loadBalancerLevel1);
  prepareGeometryFine(converterLevel1, sGeometryLevel1, circle, geomParams);
  SuperLattice<T,DESCRIPTOR> sLatticeLevel1(cuboidDecompositionLevel1,
                                            loadBalancerLevel1,
                                            3,
                                            converterLevel1);
  prepareLattice(sLatticeLevel1, converterLevel1, sGeometryLevel1, circle);

  auto coarseToFineLevel1 = refinement::lagrava::makeCoarseToFineCoupler(
    sLatticeLevel0, sGeometryLevel0,
    sLatticeLevel1, sGeometryLevel1);
  auto fineToCoarseLevel1 = refinement::lagrava::makeFineToCoarseCoupler(
    sLatticeLevel0, sGeometryLevel0,
    sLatticeLevel1, sGeometryLevel1);

  RefinedLoadBalancer<T,2> loadBalancerLevel2(cuboidDecompositionLevel1,
                                              loadBalancerLevel1,
                                              cuboidDecompositionLevel2);
  auto converterLevel2 = convectivelyRefineUnitConverter(converterLevel1, 2);
  converterLevel2.print();
  SuperGeometry<T,2> sGeometryLevel2(cuboidDecompositionLevel2, loadBalancerLevel2);
  prepareGeometryFine(converterLevel2, sGeometryLevel2, circle, geomParams);
  SuperLattice<T,DESCRIPTOR> sLatticeLevel2(cuboidDecompositionLevel2,
                                            loadBalancerLevel2,
                                            3,
                                            converterLevel2);
  prepareLattice(sLatticeLevel2, converterLevel2, sGeometryLevel2, circle);

  auto coarseToFineLevel2 = refinement::lagrava::makeCoarseToFineCoupler(
    sLatticeLevel1, sGeometryLevel1,
    sLatticeLevel2, sGeometryLevel2);
  auto fineToCoarseLevel2 = refinement::lagrava::makeFineToCoarseCoupler(
    sLatticeLevel1, sGeometryLevel1,
    sLatticeLevel2, sGeometryLevel2);

  coarseToFineLevel1->apply(meta::id<refinement::lagrava::InitializeO>{});
  coarseToFineLevel2->apply(meta::id<refinement::lagrava::InitializeO>{});

  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converterLevel0.getLatticeTime(maxPhysT),
                           sGeometryLevel0.getStatistics().getNvoxel()
                       + 2*sGeometryLevel1.getStatistics().getNvoxel()
                       + 4*sGeometryLevel2.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converterLevel0.getLatticeTime(maxPhysT); ++iT) {
    if (iT % converterLevel0.getLatticeTime(1.0) == 0) {
      writeResults(sLatticeLevel0, converterLevel0, iT, sGeometryLevel0, "level0");
      writeResults(sLatticeLevel1, converterLevel1, iT, sGeometryLevel1, "level1");
      writeResults(sLatticeLevel2, converterLevel2, iT, sGeometryLevel2, "level2");
    }

    setBoundaryValues(sLatticeLevel0, converterLevel0, iT, sGeometryLevel0);
    sLatticeLevel0.collideAndStream();

    sLatticeLevel1.collideAndStream();
    coarseToFineLevel1->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});

    {
      sLatticeLevel2.collideAndStream();
      coarseToFineLevel2->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});

      sLatticeLevel2.collideAndStream();
      coarseToFineLevel2->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});

      fineToCoarseLevel2->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
    }

    sLatticeLevel1.collideAndStream();
    coarseToFineLevel1->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});

    {
      sLatticeLevel2.collideAndStream();
      coarseToFineLevel2->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});

      sLatticeLevel2.collideAndStream();
      coarseToFineLevel2->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});

      fineToCoarseLevel2->apply(meta::id<refinement::lagrava::FineToCoarseO>{});
    }

    fineToCoarseLevel1->apply(meta::id<refinement::lagrava::FineToCoarseO>{});

    if (iT % converterLevel0.getLatticeTime(0.1) == 0) {
      timer.update(iT);
      timer.printStep();
      sLatticeLevel0.getStatistics().print(iT, converterLevel0.getPhysTime(iT));
    }
  }

  singleton::pool().wait();

  timer.stop();
  timer.printSummary();
}

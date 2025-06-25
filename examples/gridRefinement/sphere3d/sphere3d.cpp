/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
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

using T = float;
using DESCRIPTOR = D3Q19<>;

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& sGeometry)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  IndicatorCuboid3D<T> domainI(Vector<T,3>{2.5,0.5,0.5},
                               Vector<T,3>{0.0,0.0,0.0});
  IndicatorLayer3D<T> domainLayerI(domainI, converter.getPhysDeltaX());

  sGeometry.rename(0,2, domainLayerI);
  sGeometry.rename(2,1, domainI);
  sGeometry.clean();

  Vector<T,3> origin = sGeometry.getStatistics().getMinPhysR(2);
  origin[1] += converter.getPhysDeltaX()/2.;
  origin[2] += converter.getPhysDeltaX()/2.;

  Vector<T,3> extent = sGeometry.getStatistics().getMaxPhysR(2);
  extent[1] = extent[1]-origin[1]-converter.getPhysDeltaX()/2.;
  extent[2] = extent[2]-origin[2]-converter.getPhysDeltaX()/2.;

  // Set material number for inflow
  origin[0] = sGeometry.getStatistics().getMinPhysR(2)[0]-converter.getPhysDeltaX();
  extent[0] = 2*converter.getPhysDeltaX();
  IndicatorCuboid3D<T> inflow(extent, origin);
  sGeometry.rename(2,3, inflow);

  // Set material number for outflow
  origin[0] = sGeometry.getStatistics().getMaxPhysR(2)[0]-converter.getPhysDeltaX();
  extent[0] = 2*converter.getPhysDeltaX();
  IndicatorCuboid3D<T> outflow(extent, origin);
  sGeometry.rename(2,4,outflow);

  sGeometry.clean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareGeometryFine(const UnitConverter<T,DESCRIPTOR>& converter,
                         SuperGeometry<T,3>& sGeometry)
{
  OstreamManager clout( std::cout,"prepareGeometryFine");
  clout << "Prepare Geometry ..." << std::endl;

  sGeometry.rename(0,1);

  IndicatorSphere3D<T> sphereI({0.55, 0.25, 0.25}, 0.1);
  sGeometry.rename(1,2, sphereI);

  sGeometry.clean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    const UnitConverter<T,DESCRIPTOR>& converter,
                    SuperGeometry<T,3>& sGeometry)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = sGeometry.getMaterialIndicator({1});
  sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 2);

  //if interpolated boundary conditions are chosen
  boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, sGeometry, 4);

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  Vector<T,3> velocityV;
  AnalyticalConst3D<T,T> uF(velocityV);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& sLattice,
                       const UnitConverter<T,DESCRIPTOR>& converter,
                       std::size_t iT,
                       SuperGeometry<T,3>& sGeometry)
{
  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  std::size_t iTmaxStart = converter.getLatticeTime(5.0);
  std::size_t iTupdate = 30;

  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    PolynomialStartScale<T,std::size_t> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    std::size_t iTvec[1] = {iT};
    T frac[1] = {};
    StartScale(frac, iTvec);
    std::vector<T> maxVelocity( 3,0 );
    maxVelocity[0] = 2.25*frac[0]*converter.getCharLatticeVelocity();

    T distance2Wall = converter.getPhysDeltaX()/2.;
    RectanglePoiseuille3D<T> poiseuilleU( sGeometry, 3, maxVelocity, distance2Wall, distance2Wall, distance2Wall );
    sLattice.defineU( sGeometry, 3, poiseuilleU );

    clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void writeResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                  const UnitConverter<T,DESCRIPTOR>& converter,
                  std::size_t iT,
                  SuperGeometry<T,3>& sGeometry,
                  std::string name)
{
  OstreamManager clout(std::cout, name);

  if (iT == 0) {
    SuperVTMwriter3D<T> vtmWriter(name, 0);
    if (name == "level0") {
      SuperLatticeRank3D rank( sLattice );
      vtmWriter.write(rank);
    }
    SuperGeometryF<T,DESCRIPTOR::d> geometryF(sGeometry);
    geometryF.getName() = name + "_geometry";
    vtmWriter.write(geometryF);
    vtmWriter.createMasterFile();
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  sLattice.scheduleBackgroundOutputVTK([&,name,iT](auto task) {
    SuperVTMwriter3D<T> vtmWriter(name);
    SuperLatticePhysVelocity3D velocityF(sLattice, converter);
    SuperLatticePhysPressure3D pressureF(sLattice, converter);
    SuperLatticeRefinementMetricKnudsen3D qualityF(sLattice, converter);
    vtmWriter.addFunctor(qualityF);
    vtmWriter.addFunctor(velocityF);
    vtmWriter.addFunctor(pressureF);
    task(vtmWriter, iT);
  });
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  CLIreader args(argc, argv);
  const int N = args.getValueOrFallback<int>("--resolution", 11);
  const int Re = args.getValueOrFallback<int>("--reynolds", 100);
  const int maxPhysT = args.getValueOrFallback<int>("--maxPhysT", 16);

  const UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converterLevel0(
    int {N},             // resolution: number of voxels per charPhysL
    (T)   0.505,         // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   0.2,           // charPhysLength: reference length of simulation geometry
    (T)   0.2,           // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.2*0.2/Re,    // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0            // physDensity: physical density in __kg / m^3__
  );
  converterLevel0.print();

  IndicatorCuboid3D<T> domainI(Vector<T,3>{2.5,0.5,0.5},
                               Vector<T,3>{0.0,0.0,0.0});
  IndicatorLayer3D<T> domainLayerI(domainI, converterLevel0.getPhysDeltaX());
  CuboidDecomposition3D<T> cuboidDecompositionLevel0(domainLayerI, converterLevel0.getPhysDeltaX());
  cuboidDecompositionLevel0.splitFractional(0, 0, {0.15,0.5,0.35});
  cuboidDecompositionLevel0.splitFractional(1, 1, {0.1,0.8,0.1});
  cuboidDecompositionLevel0.splitFractional(3, 2, {0.1,0.8,0.1});
  cuboidDecompositionLevel0.splitFractional(5, 0, {0.5,0.5});
  cuboidDecompositionLevel0.splitFractional(7, 1, {0.5,0.5});
  cuboidDecompositionLevel0.print();

  auto cuboidDecompositionLevel1 = cuboidDecompositionLevel0;
  IndicatorCuboid3D<T> toBeRefinedI({1.0,0.1,0.1}, {0.4,0.2,0.2});
  cuboidDecompositionLevel1.remove(toBeRefinedI);
  cuboidDecompositionLevel1.refine(2);
  cuboidDecompositionLevel1.print();

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
  }

  HeuristicLoadBalancer<T> loadBalancerLevel0(cuboidDecompositionLevel0);
  SuperGeometry<T,3> sGeometryLevel0(cuboidDecompositionLevel0, loadBalancerLevel0);
  prepareGeometry(converterLevel0, sGeometryLevel0);
  SuperLattice<T,DESCRIPTOR> sLatticeLevel0(cuboidDecompositionLevel0,
                                            loadBalancerLevel0,
                                            3,
                                            converterLevel0);
  prepareLattice(sLatticeLevel0, converterLevel0, sGeometryLevel0);

  RefinedLoadBalancer<T,3> loadBalancerLevel1(cuboidDecompositionLevel0,
                                              loadBalancerLevel0,
                                              cuboidDecompositionLevel1);
  auto converterLevel1 = convectivelyRefineUnitConverter(converterLevel0, 2);
  converterLevel1.print();
  SuperGeometry<T,3> sGeometryLevel1(cuboidDecompositionLevel1, loadBalancerLevel1);
  prepareGeometryFine(converterLevel1, sGeometryLevel1);
  SuperLattice<T,DESCRIPTOR> sLatticeLevel1(cuboidDecompositionLevel1,
                                            loadBalancerLevel1,
                                            3,
                                            converterLevel1);
  prepareLattice(sLatticeLevel1, converterLevel1, sGeometryLevel1);

  auto coarseToFine = refinement::lagrava::makeCoarseToFineCoupler(
    sLatticeLevel0, sGeometryLevel0,
    sLatticeLevel1, sGeometryLevel1);
  auto fineToCoarse = refinement::lagrava::makeFineToCoarseCoupler(
    sLatticeLevel0, sGeometryLevel0,
    sLatticeLevel1, sGeometryLevel1);

  coarseToFine->apply(meta::id<refinement::lagrava::InitializeO>{});

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(converterLevel0.getLatticeTime(maxPhysT),
                           sGeometryLevel0.getStatistics().getNvoxel()
                       + 2*sGeometryLevel1.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converterLevel0.getLatticeTime(maxPhysT); ++iT) {
    if (iT % converterLevel0.getLatticeTime(1.0) == 0) {
      writeResults(sLatticeLevel0, converterLevel0, iT, sGeometryLevel0, "level0");
      writeResults(sLatticeLevel1, converterLevel1, iT, sGeometryLevel1, "level1");
    }

    setBoundaryValues(sLatticeLevel0, converterLevel0, iT, sGeometryLevel0);
    sLatticeLevel0.collideAndStream();

    sLatticeLevel1.collideAndStream();
    coarseToFine->apply(meta::id<refinement::lagrava::HalfTimeCoarseToFineO>{});

    sLatticeLevel1.collideAndStream();
    coarseToFine->apply(meta::id<refinement::lagrava::FullTimeCoarseToFineO>{});

    fineToCoarse->apply(meta::id<refinement::lagrava::FineToCoarseO>{});

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

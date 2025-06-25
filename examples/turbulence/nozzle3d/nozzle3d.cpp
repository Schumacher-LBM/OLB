/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Patrick Nathan
 *                2023-2025 Adrian Kummerlaender, Fedor Bukreev, Stephan Simonis
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

#include "fringeZone.h"

using namespace olb;
using namespace olb::descriptors;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D3Q27<>;

/// Selectable bulk models in this case
enum class BulkModel {
  RLB,
  Smagorinsky,
  LocalSmagorinsky,
  ShearSmagorinsky,
  ConsistentStrainSmagorinsky,
  Krause,
};

/// Returns BulkModel from CLI argument
BulkModel bulkModelFromString(std::string name) {
  if (name == "RLB") {
    return BulkModel::RLB;
  }
  if (name == "Smagorinsky") {
    return BulkModel::Smagorinsky;
  }
  if (name == "LocalSmagorinsky") {
    return BulkModel::LocalSmagorinsky;
  }
  if (name == "ShearSmagorinsky") {
    return BulkModel::ShearSmagorinsky;
  }
  if (name == "ConsistentStrainSmagorinsky") {
    return BulkModel::ConsistentStrainSmagorinsky;
  }
  if (name == "Krause") {
    return BulkModel::Krause;
  }
  throw std::runtime_error(name + " is not a valid BulkModel");
}

#if defined(PLATFORM_CPU_SIMD) && defined(PLATFORM_GPU_CUDA) && defined(PARALLEL_MODE_MPI)
#define HETEROGENEITY_AVAILABLE
#endif

/// Constructs indicator of inlet tube geometry
std::shared_ptr<IndicatorF3D<T>> makeInletI(const UnitConverter<T,DESCRIPTOR>& converter)
{
  Vector<T,3> origin(T(),
                      5.5*converter.getCharPhysLength()+converter.getPhysDeltaX(),
                      5.5*converter.getCharPhysLength()+converter.getPhysDeltaX());
  Vector<T,3> extend(4.*converter.getCharPhysLength()+5*converter.getPhysDeltaX(),
                      5.5*converter.getCharPhysLength()+converter.getPhysDeltaX(),
                      5.5*converter.getCharPhysLength()+converter.getPhysDeltaX());
  return std::shared_ptr<IndicatorF3D<T>>(new IndicatorCylinder3D<T>(extend, origin, converter.getCharPhysLength()));
}

/// Constructs indicator of injection tube geometry
std::shared_ptr<IndicatorF3D<T>> makeInjectionTubeI(const UnitConverter<T,DESCRIPTOR>& converter)
{
  Vector<T,3> origin(4.*converter.getCharPhysLength(),
                     5.5*converter.getCharPhysLength()+converter.getPhysDeltaX(),
                     5.5*converter.getCharPhysLength()+converter.getPhysDeltaX());
  Vector<T,3> extend(60.*converter.getCharPhysLength(),
                     5.5*converter.getCharPhysLength()+converter.getPhysDeltaX(),
                     5.5*converter.getCharPhysLength()+converter.getPhysDeltaX());
  return std::shared_ptr<IndicatorF3D<T>>(new IndicatorCylinder3D<T>(extend, origin, 5.5*converter.getCharPhysLength()));
}

std::shared_ptr<IndicatorF3D<T>> makeDomainI(const UnitConverter<T,DESCRIPTOR>& converter)
{
  auto inletCylinder = makeInletI(converter);
  auto injectionTube = makeInjectionTubeI(converter);
  return inletCylinder + injectionTube;
}

std::shared_ptr<IndicatorF3D<T>> makeBoundingI(const UnitConverter<T,DESCRIPTOR>& converter)
{
  return std::shared_ptr<IndicatorF3D<T>>(
    new IndicatorLayer3D<T>(makeDomainI(converter),
                            converter.getPhysDeltaX()));
}

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid and boundary
  superGeometry.rename(0,2);

  auto inletCylinder = makeInletI(converter);
  superGeometry.rename(2, 1, *inletCylinder);

  auto injectionTube = makeInjectionTubeI(converter);
  superGeometry.rename(2, 1, *injectionTube);

  {
    Vector<T,3> origin(converter.getPhysDeltaX(),
                       5.5*converter.getCharPhysLength()+converter.getPhysDeltaX(),
                       5.5*converter.getCharPhysLength()+converter.getPhysDeltaX());
    Vector<T,3> extend(T(),
                       5.5*converter.getCharPhysLength()+converter.getPhysDeltaX(),
                       5.5*converter.getCharPhysLength()+converter.getPhysDeltaX());

    IndicatorCylinder3D<T> cylinderIN(extend, origin, converter.getCharPhysLength());
    superGeometry.rename(1,3, cylinderIN);
  }

  {
    Vector<T,3> origin(60.*converter.getCharPhysLength()-converter.getPhysDeltaX(),
                       5.5*converter.getCharPhysLength()+converter.getPhysDeltaX(),
                       5.5*converter.getCharPhysLength()+converter.getPhysDeltaX());
    Vector<T,3> extend(60.*converter.getCharPhysLength(),
                       5.5*converter.getCharPhysLength()+converter.getPhysDeltaX(),
                       5.5*converter.getCharPhysLength()+converter.getPhysDeltaX());

    IndicatorCylinder3D<T> cylinderOUT(extend, origin, 5.5*converter.getCharPhysLength());
    superGeometry.rename(1,4, cylinderOUT);
  }

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    const UnitConverter<T,DESCRIPTOR>& converter,
                    SuperGeometry<T,3>& superGeometry,
                    BulkModel bulkModel)
{
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 3, 4});

  switch (bulkModel) {
    case BulkModel::RLB:
      sLattice.defineDynamics<RLBdynamics>(bulkIndicator);
      break;
    case BulkModel::ShearSmagorinsky:
      sLattice.defineDynamics<ShearSmagorinskyBGKdynamics>(bulkIndicator);
      sLattice.setParameter<collision::LES::SMAGORINSKY>(0.15);
      break;
    case BulkModel::Krause:
      sLattice.defineDynamics<KrauseBGKdynamics>(bulkIndicator);
      sLattice.setParameter<collision::LES::SMAGORINSKY>(0.15);
    case BulkModel::ConsistentStrainSmagorinsky:
      sLattice.defineDynamics<ConStrainSmagorinskyBGKdynamics>(bulkIndicator);
      sLattice.setParameter<collision::LES::SMAGORINSKY>(0.15);
      break;
    case BulkModel::Smagorinsky:
      sLattice.defineDynamics<SmagorinskyBGKdynamics>(bulkIndicator);
      sLattice.setParameter<collision::LES::SMAGORINSKY>(T(0.15));
      break;
    case BulkModel::LocalSmagorinsky:
    default:
      sLattice.defineDynamics<LocalSmagorinskyBGKdynamics>(bulkIndicator);
      FringeZoneSmagorinskyConstant smagorinskyFringe(converter, T{0.15});
      sLattice.defineField<collision::LES::SMAGORINSKY>(bulkIndicator, smagorinskyFringe);
      break;
  }

  // Material=2 -->bounce back
  sLattice.defineDynamics<BounceBack>(superGeometry, 2);

  const T omega = converter.getLatticeRelaxationFrequency();
  boundary::set<boundary::LocalVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedConvection>(sLattice, superGeometry, 4);

  sLattice.setParameter<descriptors::OMEGA>(omega);

  {
    auto& communicator = sLattice.getCommunicator(stage::PostCollide());
    communicator.clearRequestedCells();
    communicator.requestOverlap(1, superGeometry.getMaterialIndicator({1,2,3,4}));
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
                       SuperLattice<T,DESCRIPTOR>& sLattice,
                       SuperGeometry<T,3>& superGeometry,
                       VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>& vortex,
                       BulkModel bulkModel,
                       std::size_t iT)
{
  if (bulkModel == BulkModel::ShearSmagorinsky) {
    sLattice.setParameter<descriptors::LATTICE_TIME>(iT);
  }

  if (iT == 0) {
    AnalyticalConst3D<T,T> rhoF(1);
    Vector<T,3> velocity{};
    AnalyticalConst3D<T,T> uF(velocity);

    sLattice.defineRhoU(superGeometry.getMaterialIndicator({1, 2, 3, 4}), rhoF, uF);
    sLattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 2, 3, 4}), rhoF, uF);

    sLattice.initialize();

    superGeometry.updateStatistics();

    auto intensity = std::shared_ptr<AnalyticalF3D<T,T>>(new AnalyticalConst3D<T,T>(0.05));
    vortex.setIntensityProfile(intensity);
  }

  const T startUpFactor = 0.001 * 1.0 / converter.getCharLatticeVelocity();
  const T maxStartTPhys = (60*converter.getCharPhysLength()/converter.getCharPhysVelocity()) * startUpFactor;
  const auto maxStartT = converter.getLatticeTime(maxStartTPhys);

  auto uSol = std::shared_ptr<AnalyticalF3D<T,T>>(new CirclePowerLaw3D<T>(superGeometry, 3, converter.getCharLatticeVelocity(), 8, T()));

  if (iT <= maxStartT) {
    PolynomialStartScale<T,std::size_t> scale(maxStartT, 1);
    T frac{};
    scale(&frac, &iT);

    vortex.setVelocityProfile(converter.getConversionFactorVelocity() * frac * uSol);
    sLattice.template setProcessingContext<Array<U_PROFILE>>(ProcessingContext::Simulation);
    vortex.apply(iT);
  } else {
    vortex.apply(iT);
  }
}

void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                bool vtkOut)
{
  OstreamManager clout(std::cout,"getResults");

  if (iT==0) {
    SuperVTMwriter3D<T> vtmWriter("nozzle3d");
    SuperLatticeCuboid3D cuboid(sLattice);
    SuperLatticeRank3D rank(sLattice);
    SuperLatticePlatform platform(sLattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.write(platform);
    vtmWriter.createMasterFile();
  }

  if (iT % converter.getLatticeTime(0.1) == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    if (std::isnan(sLattice.getStatistics().getAverageRho())) {
      std::exit(-1);
    }
  }

  // Writes the vtk files
  if (vtkOut && iT % converter.getLatticeTime(1) == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtkWriter("nozzle3d");
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      SuperLatticePhysPressure3D pressure(sLattice, converter);
      vtkWriter.addFunctor(velocity);
      vtkWriter.addFunctor(pressure);
      task(vtkWriter, iT);
    });
  }
}

int main(int argc, char* argv[])
{
  CLIreader args(argc, argv);
  const int resolution = args.getValueOrFallback("--resolution", 5);
  const int maxPhysT   = args.getValueOrFallback("--max-phys-t", 150);
  const bool noResults = args.contains("--no-results");
  const BulkModel bulkModel = bulkModelFromString(args.getValueOrFallback<std::string>("--bulk-model", "LocalSmagorinsky"));
#ifdef HETEROGENEITY_AVAILABLE
  const bool heterogeneous = args.contains("--heterogeneous");
#endif

  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    int {resolution}, // resolution: number of voxels per charPhysL
    (T)   0.500018,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,          // charPhysLength: reference length of simulation geometry
    (T)   1,          // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.0002,     // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0         // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  auto boundingI = makeBoundingI(converter);

  clout << "Instantiate cuboid geometry... " << std::endl;
  std::unique_ptr<CuboidDecomposition3D<T>> cuboidDecomposition;
  if (args.contains("--load-decomposition-from")) {
    const std::string xmlPath = args.getValueOrFallback<std::string>("--load-decomposition-from", "");
    clout << "Load decomposition from " << xmlPath << std::endl;
    cuboidDecomposition = createCuboidDecomposition<T,3>(xmlPath);
    if (!cuboidDecomposition->tryRefineTo(converter.getPhysDeltaX())) {
      throw std::invalid_argument("Cuboid geometry is not refineable to reach goal deltaR");
    }
    continueMinimizeByVolume(*cuboidDecomposition, *boundingI, singleton::mpi().getSize());
  } else {
    cuboidDecomposition.reset(new CuboidDecomposition3D<T>(*boundingI, converter.getPhysDeltaX(), singleton::mpi().getSize()));
  }
  clout << "Done." << std::endl;

  clout << "Balance load... " << std::endl;
  std::unique_ptr<LoadBalancer<T>> loadBalancer;
#ifdef HETEROGENEITY_AVAILABLE
  if (heterogeneous) {
    loadBalancer.reset(new OrthogonalHeterogeneousLoadBalancer<T>(*cuboidDecomposition, 0.99));
  } else {
    loadBalancer.reset(new BlockLoadBalancer<T>(*cuboidDecomposition));
  }
#else
  loadBalancer.reset(new BlockLoadBalancer<T>(*cuboidDecomposition));
#endif
  clout << "Done." << std::endl;

  // === 2nd Step: Prepare Geometry ===
  SuperGeometry<T,3> superGeometry(*cuboidDecomposition, *loadBalancer);
  prepareGeometry(converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(sLattice, converter, superGeometry, bulkModel);

  // === Setup turbulent inlet condition ===
  Vector<T,3> originI(1./resolution, 5.5+1./resolution, 5.5+1./resolution);
  Vector<T,3> extendI(0., 5.5+1./resolution, 5.5+1./resolution);
  Vector<T,3> inflowAxis{1, 0, 0};
  IndicatorCylinder3D<T> cylinderIN(extendI, originI, 1.);
  VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR> vortex(
    superGeometry.getMaterialIndicator(3),
    cylinderIN,
    converter,
    sLattice,
    50,                                // nSeeds
    0.1,                               // nTime (s)
    converter.getCharPhysLength()*0.1, // sigma
    inflowAxis);

  // === 4th Step: Main Loop with Timer ===
  setBoundaryValues(converter, sLattice, superGeometry, vortex, bulkModel, 0);

  sLattice.writeSummary();

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());

  setBoundaryValues(converter, sLattice, superGeometry, vortex, bulkModel, 0);
  sLattice.collideAndStream();
  getResults(sLattice, converter, 0, superGeometry, timer, !noResults);

  timer.start();

  for (std::size_t iT = 1; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
    setBoundaryValues(converter, sLattice, superGeometry, vortex, bulkModel, iT);
    sLattice.collideAndStream();
    getResults(sLattice, converter, iT, superGeometry, timer, !noResults);
  }

  timer.stop();
  timer.printSummary();

  sLattice.setProcessingContext(ProcessingContext::Evaluation);
}

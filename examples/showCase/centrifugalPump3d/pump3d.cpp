/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Heiko Kipp, Fedor Bukreev, Adrian Kummerländer, Mathias J. Krause
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
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;

using DESCRIPTOR = descriptors::D3Q19<
  descriptors::POROSITY,
  descriptors::VELOCITY,

  fields::fsi::ELEMENT_TAG,
  fields::fsi::ELEMENT_TORQUE
>;

using BulkDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  typename momenta::Tuple<
    momenta::BulkDensity,
    momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>,
    momenta::BulkStress,
    momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>
>;

// Discretization Settings
const T charLatticeVelocity = 0.01;
const T rotationFreq = 24.1667; // 1/s ACHTUNG 1450rpm => 24.166Hz
const T bladeRadius = 0.15; //Durchmesser Lara 400mm
const T volume = 0.00337; // pump volume

// Time Settings
const T startT = 0.01;
const T maxPhysT = 0.66; //16 Umdrehungen
const T iTwrite = 2.3e-4; //2° Schritt

// Fluid Settings
const T physDensity = 998.21; //Wasser 20°C
const T physViscosity = 1E-5;

// Characteristic Quantities
const T charPhysLength = 0.15;
const T angularVelocity = 2. * std::numbers::pi_v<T> * rotationFreq;
const T charPhysVelocity = angularVelocity * bladeRadius;    // Assumed maximal velocity
const T turbIntensity = 0.015; //1.5% Turbulenzgrad

// Prepare geometry
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader, SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,stlReader );

  superGeometry.clean();
  superGeometry.innerClean();

  T dx = converter.getPhysDeltaX();

  Vector<T,3> inletI(0., 0., indicator.getMax()[2]-20*dx);
  Vector<T,3> inletO(0., 0., indicator.getMax()[2]+20*dx);
  IndicatorCylinder3D<T> inlet( inletI, inletO, 0.05 );
  superGeometry.rename( 2,3,1,inlet );

  Vector<T,3> outletI(-0.168, indicator.getMax()[1]-20*dx, 0.0082);
  Vector<T,3> outletO(-0.168, indicator.getMax()[1]+20*dx, 0.0082);
  IndicatorCylinder3D<T> outlet( outletI, outletO, 0.045);
  superGeometry.rename( 2,4,1,outlet );

  superGeometry.clean();
  superGeometry.innerClean();

  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}


// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>& sLattice,
                    UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  clout << "setting Velocity Boundaries ..." << std::endl;

  /// Material=0 -->do nothing
  T omega = converter.getLatticeRelaxationFrequency();
  sLattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator({1,3,4}));
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  sLattice.setParameter<descriptors::OMEGA>(omega);
  sLattice.setParameter<collision::LES::SMAGORINSKY>(0.15);

  {
    AnalyticalConst3D<T,T> porosityF(1);
    sLattice.defineField<descriptors::POROSITY>(superGeometry.getMaterialIndicator({0,1,2}), porosityF);
  }

  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}


//Set Boundary Values
void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& sLattice,
                       UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                       SuperGeometry<T,3>& superGeometry,
                       VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>& vortex)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  if (iT == 0) {
    AnalyticalConst3D<T, T> zero(0.);
    AnalyticalConst3D<T, T> one(1.);
    AnalyticalConst3D<T, T> u(0., 0., 0.);
    auto bulkIndicator = superGeometry.getMaterialIndicator({0,1,2,3,4});
    sLattice.defineField<VELOCITY>(bulkIndicator, u);

    //Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU(superGeometry.getMaterialIndicator({0,1,2,3,4}), one, u);
    sLattice.iniEquilibrium(superGeometry.getMaterialIndicator({0,1,2,3,4}), one, u);

    auto intensity = std::shared_ptr<AnalyticalF3D<T,T>>(new AnalyticalConst3D<T,T>(turbIntensity));
    vortex.setIntensityProfile(intensity);

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }

  T maxVelocity = 3.67286; //max velocity at inlet in m/s sind dann ca. 480m³/h => Q=A*c
  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime(2.68e-3);

  if (iT <= iTmaxStart && iT%10 == 0) {
    // Smooth start curve, sinus
    SinusStartScale<T,int> startScale(iTmaxStart, T(1));

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1]= {iT};
    T frac[1]= {};
    startScale( frac,iTvec );

    maxVelocity *= frac[0];

    auto poiseuilleU = std::shared_ptr<AnalyticalF3D<T,T>>(new CirclePowerLaw3D<T>(superGeometry, 3, maxVelocity, 2, T())); //2 ist power-law coeff. vortex methode s. ansys

    vortex.setVelocityProfile(poiseuilleU);
    sLattice.template setProcessingContext<Array<U_PROFILE>>(ProcessingContext::Simulation);
    vortex.apply(iT);
  }else{
    vortex.apply(iT);
  }
}


/// Computes the pressure drop between the voxels before and after the cylinder
void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter,
                int iT,
                SuperGeometry<T,3>& superGeometry,
                util::Timer<T>& timer)
{
  OstreamManager clout(std::cout, "getResults");

  if (iT == 0) {
    {
      SuperVTMwriter3D<T> vtmWriter("centrifugalPump");
      SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
      SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);
      vtmWriter.write(cuboid);
      vtmWriter.write(rank);
      vtmWriter.createMasterFile();
    }
    {
      SuperVTMwriter3D<T> vtmWriter("structure");
      vtmWriter.createMasterFile();
    }
  }

  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtkWriter("centrifugalPump");
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure(sLattice, converter);
      vtkWriter.addFunctor(velocity);
      vtkWriter.addFunctor(pressure);
      task(vtkWriter, iT);
    });
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtkWriter("structure");
      SuperLatticeField3D<T,DESCRIPTOR,descriptors::POROSITY> porosity(sLattice);
      vtkWriter.addFunctor(porosity);
      task(vtkWriter, iT);
    });

  }

  /// Writes output on the console
  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    if (std::isnan(sLattice.getStatistics().getAverageRho())) {
      std::exit(-1);
    }
  }
}

int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  initialize(&argc, &argv);

  CLIreader args(argc, argv);
  const int resolution = args.getValueOrFallback<int>("--resolution", 60);

  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR> converter(
    (int)   resolution,           //resolution
    ( T )   charLatticeVelocity,  //charLatticeVelocity
    ( T )   charPhysLength,       //charPhysLength
    ( T )   charPhysVelocity,     //charPhysVelocity
    ( T )   physViscosity,        //physViscosity
    ( T )   physDensity           //physDensity
  );
  converter.print();

  /// === 2rd Step: Prepare Geometry ===
  /// Instantiation of a cuboidDecomposition with weights
  STLreader<T> stlReader( "casing.stl", converter.getPhysDeltaX(), 0.001, RayMode::Robust);
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getPhysDeltaX() );

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 2*singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition3D<T> cuboidDecomposition( extendedDomain, converter.getConversionFactorLength(), noOfCuboids, "volume" );
  cuboidDecomposition.print();

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,3> superGeometry( cuboidDecomposition, loadBalancer, 2 );
  prepareGeometry(converter, extendedDomain, stlReader, superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);

  prepareLattice(sLattice, converter, superGeometry);

  const T dx = converter.getPhysDeltaX();
  Vector<T,3> inflowAxis{0, 0, -1};
  Vector<T,3> inletI(0., 0., stlReader.getMax()[2]);
  Vector<T,3> inletO(0., 0., stlReader.getMax()[2]+2.*dx);
  IndicatorCylinder3D<T> inflow( inletI, inletO, 0.035 );

  VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR> vortex(
    superGeometry.getMaterialIndicator(3),
    inflow,
    converter,
    sLattice,
    200,       // nSeeds
    0.0005,    // nTime (s)
    0.108*0.1, // sigma
    inflowAxis);

  SuperPorousElementEmbeddingO fsiEmbeddingO(sLattice);
  fsiEmbeddingO.registerElementType<ReferenceLatticePorosityF>(superGeometry.getMaterialIndicator(1));

  SuperPorousElementReductionO<T,DESCRIPTOR,fields::fsi::ELEMENT_TORQUE> fsiReductionO(
    sLattice,
    superGeometry.getMaterialIndicator(1));
  fsiReductionO.resize(1);
  fsiReductionO.addCollectionO(meta::id<CollectPorousBoundaryTorqueO>{});

  std::unique_ptr<BlockD<T,descriptors::fsi::REFERENCE_POROSITY_ELEMENT<DESCRIPTOR::d>>> rotorD;

  if (fsiReductionO.rankDoesFSI()) {
    clout << "Set up rotor element" << std::endl;
    const T elementDeltaX = converter.getPhysDeltaX() / 2;
    std::shared_ptr<IndicatorF3D<T>> rotorI(new STLreader<T>("rotor.stl", elementDeltaX, 0.001, RayMode::Robust));
    std::shared_ptr<IndicatorF3D<T>> rotorBoundingI(new IndicatorLayer3D<T>(*rotorI, converter.getPhysDeltaX()));
    Vector<std::size_t,DESCRIPTOR::d> rotorDim = (rotorBoundingI->getMax() - rotorBoundingI->getMin())
                                               / elementDeltaX + 1.5;
    clout << rotorI->getMin() << " " << rotorI->getMax() << std::endl;
    rotorD.reset(makeSharedBlockD<T,descriptors::fsi::REFERENCE_POROSITY_ELEMENT<DESCRIPTOR::d>>(
      loadBalancer, rotorDim).release());
    rotorD->forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
      bool inside{};
      Vector<T,DESCRIPTOR::d> physR = rotorBoundingI->getMin()
                                    + latticeR * elementDeltaX;
      rotorI->operator()(&inside, physR.data());
      if (inside) {
        rotorD->setField<descriptors::POROSITY>(latticeR, 0);
      } else {
        rotorD->setField<descriptors::POROSITY>(latticeR, 1);
      }
    });
    rotorD->setProcessingContext(ProcessingContext::Simulation);

    clout << rotorBoundingI->getMin() << std::endl;

    {
      Vector lower{rotorBoundingI->getMin()[0],
                   rotorBoundingI->getMin()[1],
                   rotorBoundingI->getMin()[2]};
      auto params = makeParametersD<T,descriptors::SPATIAL_DESCRIPTOR<3>>(
        fields::fsi::ELEMENT_TAG{},   1,
        fields::fsi::ELEMENT_LOWER{}, lower,
        fields::fsi::ELEMENT_PIVOT{}, Vector{0.,0.,0.},
        fields::fsi::ELEMENT_ROTATION{}, 0,
        fields::fsi::ELEMENT_REFERENCE_DELTA_X{}, elementDeltaX,
        fields::fsi::ELEMENT_REFERENCE_EXTENT{}, rotorDim,
        fields::fsi::ELEMENT_REFERENCE_POROSITY{}, rotorD->getField<descriptors::POROSITY>()
      );
      fsiEmbeddingO.add(params);
    }
  }

  sLattice.setParameter<fields::converter::PHYS_LENGTH>(
    1/converter.getConversionFactorLength());
  sLattice.setParameter<fields::converter::PHYS_VELOCITY>(
    1/converter.getConversionFactorVelocity());
  sLattice.setParameter<fields::converter::PHYS_DELTA_X>(
    converter.getPhysDeltaX());

  sLattice.setProcessingContext(ProcessingContext::Simulation);
  fsiEmbeddingO.initialize();

  CSV<T> torqueWriter("torque_n",
                      ';', {"time", "torque", "angle"},
                      ".csv");

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT),
                       superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
    setBoundaryValues(sLattice, converter, iT, superGeometry, vortex);
    sLattice.collideAndStream();
    getResults(sLattice, converter, iT, superGeometry, timer);

    fsiReductionO.apply();

    if (fsiReductionO.rankDoesFSI()) {
      const unsigned iElement = 1;
      auto t = converter.getPhysTorque(fsiReductionO.getField<fields::fsi::ELEMENT_TORQUE>(iElement));
      auto angle = fsiEmbeddingO.getField<fields::fsi::ELEMENT_ROTATION>(iElement-1);

      torqueWriter.writeDataFile(converter.getPhysTime(iT), {t[2], angle[2]});

      Vector angularU{0.0, 0.0, -angularVelocity};
      if (iT <= converter.getLatticeTime(startT)) {
        SinusStartScale<T,std::size_t> startScaleF(converter.getLatticeTime(startT), T{1});
        T frac{};
        startScaleF(&frac, &iT);
        angularU *= frac;
      }
      angle += angularU * converter.getPhysDeltaT();
      angle = util::fmod(angle, 2*std::numbers::pi_v<T>);

      fsiEmbeddingO.setField<fields::fsi::ELEMENT_ROTATION>(iElement-1, angle);
      fsiEmbeddingO.setField<fields::fsi::ELEMENT_U_ROTATION>(iElement-1, angularU);

      fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_ROTATION>>(ProcessingContext::Simulation);
      fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_U_ROTATION>>(ProcessingContext::Simulation);
    }

    fsiEmbeddingO.apply();
  }

  timer.stop();
  timer.printSummary();
}

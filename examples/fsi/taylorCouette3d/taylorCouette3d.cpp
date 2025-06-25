/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Fedor Bukreev, Adrian Kummerlaende
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
 *  Boston, MA  02110-1201, USA.
 */

// Taylor-Couette flow with rotating inner and static outer cylinders
// Literature with parameters and results for this example can be found here:
// https://doi.org/10.1146/annurev-fluid-122414-034353
// https://doi.org/10.1017/jfm.2014.134
// https://doi.org/10.1063/1.4863312

#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;

using T = float;
using DESCRIPTOR = descriptors::D3Q19<
  descriptors::POROSITY,
  descriptors::VELOCITY,
  fields::fsi::ELEMENT_TAG,
  fields::fsi::ELEMENT_FORCE,
  fields::fsi::ELEMENT_TORQUE,
  descriptors::AVERAGE_VELOCITY
>;

const int N = 40;
const T eta = 0.5; // \eta = r_i/r_o
const T innerR = 0.3;
const T outerR = innerR/eta;
const T gapW = outerR - innerR;
const T length = 2.*std::numbers::pi*gapW/3.;
const T Ta = util::pow(10.,10.); // Taylor number
const T viscosity = 1.e-5;
const T sigma = util::pow((innerR + outerR)/2./util::sqrt(innerR*outerR),4.);
const T rotateOmega = util::sqrt(Ta*viscosity*viscosity/(0.25*sigma*(outerR-innerR)*(outerR-innerR)*(outerR+innerR)*(outerR+innerR)));
const T charPhysU = rotateOmega*innerR;
Vector<T,3> rotateU{0.,0.,-rotateOmega};
const T turbI = 0.05;
const T maxPhysT = 1000*gapW/charPhysU;
std::default_random_engine generator(0x1337533DAAAAAAAA);

template <typename T, typename S>
class TaylorCouette3D : public AnalyticalF3D<T,S> {
protected:
  T turbulenceIntensity;
  T maxVelocity;
  T a;
  T b;

public:
  TaylorCouette3D(UnitConverter<T,DESCRIPTOR> const& converter) : AnalyticalF3D<T,S>(3)
  {
    turbulenceIntensity = turbI;
    maxVelocity = converter.getLatticeVelocity(charPhysU);
    a = -1.;
    b = 1.;
  };

  bool operator()(T output[], const S input[])
  {
    std::uniform_real_distribution<BaseType<T>> distribution(a, b);
    T nRandom1 = distribution(generator);
    T nRandom2 = distribution(generator);
    T nRandom3 = distribution(generator);

    output[0] = 0.;
    output[1] = 0.;
    output[2] = 0.;

    T radius = util::sqrt(util::pow(input[0] - outerR,2) + util::pow(input[1] - outerR,2));

    if( radius >= innerR && radius <= outerR ){

    T u_calc_tang = maxVelocity*(outerR - radius)/gapW;
    T u_calc_x = u_calc_tang*util::sin(std::acos((input[0] - outerR)/radius));
    T u_calc_y = -u_calc_tang*util::cos(std::acos((input[0] - outerR)/radius));

    if(input[1] - outerR < T(0)) u_calc_x *= T(-1);

    output[0] = turbulenceIntensity*nRandom1*u_calc_x + u_calc_x;
    output[1] = turbulenceIntensity*nRandom2*u_calc_y + u_calc_y;
    output[2] = turbulenceIntensity*nRandom3*u_calc_tang;
    }else{
    output[0] = 0.;
    output[1] = 0.;
    output[2] = 0.;
    }

    return true;
  };
};

std::shared_ptr<IndicatorF3D<T>> makeDomainI() {
  std::shared_ptr<IndicatorF3D<T>> boxI(new IndicatorCuboid3D<T>({2*outerR,2*outerR,length},{0,0,0}));
  return boxI;
}

std::shared_ptr<IndicatorF3D<T>> makeFsiRegionI() {
  return makeDomainI();
}

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  auto domainI = makeDomainI();
  auto fsiRegionI = makeFsiRegionI();

  superGeometry.rename(0, 3);
  superGeometry.rename(3, 2, {1,1,0});

  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(SuperLattice<T, DESCRIPTOR>& sLattice,
                    const UnitConverter<T, DESCRIPTOR>& converter,
                    SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  WallModelParameters<T> wallModelParameters;
  wallModelParameters.bodyForce = false;
  wallModelParameters.rhoMethod = 0;
  wallModelParameters.fNeqMethod = 0;
  wallModelParameters.samplingCellDistance = 2.5;
  wallModelParameters.interpolateSampleVelocity = true;
  wallModelParameters.useVanDriest = false;
  wallModelParameters.wallFunctionProfile = 1;
  wallModelParameters.movingWall = true;
  wallModelParameters.averageVelocity = true;
  setTurbulentWallModelDynamics(sLattice, superGeometry.getMaterialIndicator({1,2}), wallModelParameters);

  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);

  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 3);
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  AnalyticalConst3D<T,T> rhoF(1);
  AnalyticalConst3D<T,T> u0(0.,0.,0.);
  TaylorCouette3D<T,T> uCalc(converter);
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,2});
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uCalc);
  sLattice.defineRhoU(bulkIndicator, rhoF, uCalc);
  sLattice.defineField<descriptors::AVERAGE_VELOCITY>(superGeometry.getMaterialIndicator({0,1,2,3}), u0);

  {
    // set everywhere 0 porosity
    AnalyticalConst3D<T,T> porosityF(0);
    sLattice.defineField<descriptors::POROSITY>(superGeometry.getMaterialIndicator({0,1,2,3}), porosityF);
  }
  {
    // set porosity & wall-model-porosity 1 in the bulk and rotating inner cylinder
    AnalyticalConst3D<T,T> porosityF(1);
    IndicatorCylinder3D<T> cylinderI({outerR - T{0.25}*converter.getPhysDeltaX(), outerR - T{0.25}*converter.getPhysDeltaX(), -T{5.}*converter.getPhysDeltaX()},{outerR - T{0.25}*converter.getPhysDeltaX(), outerR - T{0.25}*converter.getPhysDeltaX(), length+T{10.}*converter.getPhysDeltaX()}, outerR);
    SuperIndicatorFfromIndicatorF3D<T> cylinderIndicatorF(cylinderI, superGeometry);
    sLattice.defineField<descriptors::POROSITY>(cylinderIndicatorF, porosityF);
    sLattice.defineField<descriptors::WMPOROSITY>(bulkIndicator, porosityF);
  }
  {
    // set Y1 vectors on the outer cylinder wall
    IndicatorCuboid3D<T> cuboid({2*outerR,2*outerR,length},{0,0,0});
    std::shared_ptr<IndicatorF3D<T>> cuboidLayer(new IndicatorLayer3D<T>(cuboid, T{5.}*converter.getPhysDeltaX()));
    std::shared_ptr<IndicatorF3D<T>> cylinderOutI(new IndicatorCylinder3D<T>({outerR - T{0.25}*converter.getPhysDeltaX(), outerR - T{0.25}*converter.getPhysDeltaX(), -T{5.}*converter.getPhysDeltaX()},{outerR - T{0.25}*converter.getPhysDeltaX(), outerR - T{0.25}*converter.getPhysDeltaX(), length+T{10.}*converter.getPhysDeltaX()}, outerR));
    SuperIndicatorFfromIndicatorF3D<T> cylinderIndicatorOutF(cylinderOutI, superGeometry);
    SuperIndicatorFfromIndicatorF3D<T> outerF((cuboidLayer - cylinderOutI), superGeometry);

    setWallDistance(sLattice, outerF, cylinderIndicatorOutF, cylinderOutI.get());

    sLattice.addPostProcessor(superGeometry.getMaterialIndicator(2),
                              meta::id<TurbulentWallModelPostProcessor<false,true,false,1,true>>{});
    sLattice.addPostProcessor(superGeometry.getMaterialIndicator({1,2}),
                              meta::id<TurbulentWallModelPorousFneqFDMPostProcessor<0>>{});

    sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
    sLattice.setParameter<collision::LES::SMAGORINSKY>(0.2);
    sLattice.setParameter<descriptors::SAMPLING_DISTANCE>(wallModelParameters.samplingCellDistance);
  }

  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults(SuperLattice<T,DESCRIPTOR>& sLattice,
                const UnitConverter<T,DESCRIPTOR>& converter,
                SuperGeometry<T,3>& superGeometry,
                util::Timer<T>& timer,
                std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");

  const int vtkIter  = converter.getLatticeTime(1./30);
  const int statIter = converter.getLatticeTime(1./60);

  if (iT == 0) {
    {
      SuperVTMwriter3D<T> vtmWriter("flow");
      vtmWriter.createMasterFile();
    }
  }

  if (iT % statIter == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    if (std::isnan(sLattice.getStatistics().getAverageRho())) {
      std::exit(-1);
    }
  }

  int iTstartAvg = converter.getLatticeTime(200.*gapW/charPhysU);
  if (iT == iTstartAvg) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticeVelocity3D<T,DESCRIPTOR> latticeVelocity(sLattice);
    sLattice.defineField<descriptors::AVERAGE_VELOCITY>(superGeometry.getMaterialIndicator({1,2,3,4,5,6,7}), latticeVelocity);
  }
  if (iT < iTstartAvg) {
    sLattice.setParameter<descriptors::LATTICE_TIME>(2);
  }
  else {
    sLattice.setParameter<descriptors::LATTICE_TIME>(iT - iTstartAvg + 1);
  }

  if (iT % vtkIter == 0 && iT > converter.getLatticeTime(1.)) {
    sLattice.executePostProcessors(stage::Evaluation{});
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task)
    {
      SuperVTMwriter3D<T> vtkWriter("flow");
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      SuperLatticePhysPressure3D pressure(sLattice, converter);
      SuperLatticePhysField3D<T,DESCRIPTOR,descriptors::AVERAGE_VELOCITY> velocityAv(sLattice, converter.getConversionFactorVelocity());
      vtkWriter.addFunctor(pressure);
      vtkWriter.addFunctor(velocity);
      vtkWriter.addFunctor(velocityAv);
      task(vtkWriter, iT);
    });
  }
}

int main(int argc, char* argv[]) {
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );

  OstreamManager clout(std::cout, "main");
  OstreamManager clresult(std::cout, "result");

  const UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> converter(
    (int) N,                // resolution of gap
    (T)   0.1,              // CFL
    (T)   gapW,             // charPhysLength: reference length of simulation geometry
    (T)   charPhysU,        // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   viscosity,        // physViscosity: physical kinematic viscosity in __m^3 / s__
    (T)   1.                // physDensity: physical density in __kg / m^3__
  );
  converter.print();
  clout << "maxPhysT: " << maxPhysT << std::endl;

  auto domainI = makeDomainI();
  std::shared_ptr<IndicatorF3D<T>> paddedDomainI(new IndicatorLayer3D<T>(*domainI, 3.*converter.getPhysDeltaX()));

  CuboidDecomposition3D<T> cDecomposition(*paddedDomainI, converter.getPhysDeltaX(), singleton::mpi().getSize());
  cDecomposition.setPeriodicity({false, false, true});
  BlockLoadBalancer<T> loadBalancer(cDecomposition);

  SuperGeometry<T,3> superGeometry(cDecomposition, loadBalancer, 5);
  prepareGeometry(converter, superGeometry);

  SuperLattice<T,DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(sLattice, converter, superGeometry);

  SuperPorousElementEmbeddingO fsiEmbeddingO(sLattice);
  fsiEmbeddingO.registerElementType<CylinderWithWallModelPorosityF>(superGeometry.getMaterialIndicator(2));

  Vector center{outerR- T{0.25}*converter.getPhysDeltaX(), outerR- T{0.25}*converter.getPhysDeltaX(), -T{5}*converter.getPhysDeltaX()};
  Vector<T,3> axis{0,0,1};
  auto valveParams = makeParametersD<T,descriptors::SPATIAL_DESCRIPTOR<3>>(
    fields::fsi::ELEMENT_TAG{},      1,
    fields::fsi::ELEMENT_PIVOT{},    center,
    fields::fsi::ELEMENT_ROTATION{}, 0,
    CylinderPorosityF::RADIUS{},     innerR,
    CylinderPorosityF::CENTER{},     center,
    CylinderPorosityF::LENGTH{},     (length+T{10}*converter.getPhysDeltaX()),
    CylinderPorosityF::AXIS{},       axis
  );
  fsiEmbeddingO.add(valveParams);

  SuperPorousElementReductionO<T,DESCRIPTOR,fields::fsi::ELEMENT_FORCE,fields::fsi::ELEMENT_TORQUE> fsiReductionO(
    sLattice,
    superGeometry.getMaterialIndicator(2));
  fsiReductionO.resize(1);
  fsiReductionO.addCollectionO(meta::id<CollectPorousBoundaryForceAndTorqueO>{});

  sLattice.setParameter<fields::converter::PHYS_LENGTH>(
    1/converter.getConversionFactorLength());
  sLattice.setParameter<fields::converter::PHYS_VELOCITY>(
    1/converter.getConversionFactorVelocity());
  sLattice.setParameter<fields::converter::PHYS_DELTA_X>(
    converter.getPhysDeltaX());

  sLattice.setProcessingContext(ProcessingContext::Simulation);
  fsiEmbeddingO.initialize();

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT),
                       superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
    fsiEmbeddingO.apply();

    sLattice.collideAndStream();

    sLattice.stripeOffDensityOffset(sLattice.getStatistics().getAverageRho()-(T)1);

    fsiReductionO.apply();

    if (fsiReductionO.rankDoesFSI()) {
      for (unsigned iElement=1; iElement <= fsiReductionO.getElementCount(); ++iElement) {
        fsiEmbeddingO.setField<fields::fsi::ELEMENT_ROTATION>(iElement-1, 0);
        fsiEmbeddingO.setField<fields::fsi::ELEMENT_U_ROTATION>(iElement-1, rotateU);
      }

      fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_ROTATION>>(ProcessingContext::Simulation);
      fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_U_ROTATION>>(ProcessingContext::Simulation);
    }

    getResults(sLattice, converter, superGeometry, timer, iT);
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  timer.stop();
  timer.printSummary();
}

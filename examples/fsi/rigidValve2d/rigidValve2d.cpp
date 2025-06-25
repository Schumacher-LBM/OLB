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

/**
 * This example implements the model for predicting the dynamic
 * response of a mechanical heart valve described by Stijnen et al.
 * DOI: 10.1016/j.jfluidstructs.2004.04.007. 2004.
 **/

#include <olb.h>

using namespace olb;

using T = float;
using DESCRIPTOR = descriptors::D2Q9<
  descriptors::POROSITY,
  descriptors::VELOCITY,
  fields::fsi::ELEMENT_TAG,
  fields::fsi::ELEMENT_FORCE,
  fields::fsi::ELEMENT_TORQUE
>;

using BulkDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::BulkTuple,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::fsi::HLBM
>;

const T channelH = 0.02;
const T channelL = 6*channelH;

const T cavityR  = channelH;
const T cavityX  = channelH+cavityR;
const T cavityY  = channelH;

const T valveL = 1.07*channelH;
const T valveMomentaOfInertia = 1100*util::pow(valveL,3)*(0.001)/3;

const T T_p = 2.45;
const T u_amp = 0.11;
const T u_ref = 0.04;

const T maxPhysT = T_p * 10;

std::shared_ptr<IndicatorF2D<T>> makeDomainI() {
  std::shared_ptr<IndicatorF2D<T>> channelI(new IndicatorCuboid2D<T>({channelL,channelH},{0,0}));
  std::shared_ptr<IndicatorF2D<T>> cavityI(new IndicatorCircle2D<T>({cavityX,cavityY}, cavityR));

  std::shared_ptr<IndicatorF2D<T>> cornerBoxI(new IndicatorCuboid2D<T>({0.00721568,0.00545155},{0.0592418,channelH}));
  std::shared_ptr<IndicatorF2D<T>> cornerCircleI(new IndicatorCircle2D<T>({0.0664575,0.0275},0.0075));
  auto cornerI = cornerBoxI - (cornerBoxI * cornerCircleI);

  return channelI+cavityI+cornerI;
}

std::shared_ptr<IndicatorF2D<T>> makeFsiRegionI() {
  std::shared_ptr<IndicatorF2D<T>> regionI(new IndicatorCuboid2D<T>({T{1.2}*cavityR,2*channelH},{channelH-T{0.1}*cavityR,0}));
  return regionI;
}

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,2>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  auto domainI = makeDomainI();
  auto fsiRegionI = makeFsiRegionI();

  superGeometry.rename(0, 3);
  superGeometry.rename(3, 1, domainI);
  superGeometry.rename(1, 2, fsiRegionI);

  {
    Vector<T,2> extend(channelL, channelH);
    Vector<T,2> origin;
    extend[0] = 2.*converter.getPhysDeltaX();
    origin[0] = -converter.getPhysDeltaX();
    IndicatorCuboid2D<T> inflow(extend, origin);
    superGeometry.rename(3,4,1, inflow);
  }

  // Set material number for outflow
  {
    Vector<T,2> extend(channelL, channelH);
    Vector<T,2> origin;
    origin[0] = channelL-converter.getPhysDeltaX();
    IndicatorCuboid2D<T> outflow(extend, origin);
    superGeometry.rename(3,5,1, outflow);
  }

  superGeometry.clean(true, {1,2});
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(SuperLattice<T, DESCRIPTOR>& sLattice,
                    const UnitConverter<T, DESCRIPTOR>& converter,
                    SuperGeometry<T,2>& superGeometry)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);
  sLattice.defineDynamics<BulkDynamics>(superGeometry, 2);

  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 4);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 5);

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  AnalyticalConst2D<T,T> rhoF(1);
  AnalyticalConst2D<T,T> uF(0,0);
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,2});
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);
  sLattice.defineRhoU(bulkIndicator, rhoF, uF);

  {
    AnalyticalConst2D<T,T> porosityF(1);
    sLattice.defineField<descriptors::POROSITY>(bulkIndicator, porosityF);
  }
  {
    AnalyticalConst2D<T,T> porosityF(0);
    sLattice.defineField<descriptors::POROSITY>(superGeometry.getMaterialIndicator({0,3}), porosityF);
  }

  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
                       SuperLattice<T, DESCRIPTOR>& sLattice,
                       SuperGeometry<T,2>& superGeometry,
                       std::size_t iT)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  const std::size_t iTupdate = converter.getLatticeTime(T_p/1000);
  const std::size_t iTstart = converter.getLatticeTime(T_p);

  if (iT % iTupdate == 0 && iT < iTstart) {
    PolynomialStartScale<T,int> scale(iTstart, T{1});
    int iTvec = static_cast<int>(iT);
    T frac[1] = {};
    scale(frac, &iTvec);
    const T targetVelocityX = frac[0]*converter.getLatticeVelocity(u_ref);
    AnalyticalConst2D<T,T> poiseuilleU(targetVelocityX, 0);
    sLattice.defineU(superGeometry, 4, poiseuilleU);
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
  else if (iT % iTupdate == 0 && iT >= iTstart) {
    T u { };
    const T t = std::fmod(converter.getPhysTime(iT), T_p);
    if (t / T_p <= 0.37) {
      u = u_ref + u_amp * util::sin(((2*std::numbers::pi_v<T>*t)/T_p)/0.74);
    } else if (t / T_p > 0.37 && t / T_p <= 1) {
      u = u_ref + 0.5 * u_amp * util::sin(2*std::numbers::pi_v<T>*(t/T_p+0.26)/1.26);
    }

    const T targetVelocityX = converter.getLatticeVelocity(u);
    AnalyticalConst2D<T,T> poiseuilleU(targetVelocityX, 0);
    sLattice.defineU(superGeometry, 4, poiseuilleU);
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults(SuperLattice<T,DESCRIPTOR>& superLattice,
                const UnitConverter<T,DESCRIPTOR>& converter,
                SuperGeometry<T,2>& superGeometry,
                util::Timer<T>& timer,
                std::size_t iT,
                bool vtkEnabled)
{
  OstreamManager clout(std::cout, "getResults");

  const int vtkIter  = converter.getLatticeTime(0.10);
  const int statIter = converter.getLatticeTime(0.01);

  if (iT == 0 && vtkEnabled) {
    {
      SuperVTMwriter2D<T> vtmWriter("flow");
      vtmWriter.createMasterFile();
    }
    {
      SuperVTMwriter2D<T> vtmWriter("structure");
      vtmWriter.createMasterFile();
    }
  }

  if (iT % statIter == 0) {
    timer.update(iT);
    timer.printStep();
    superLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    if (std::isnan(superLattice.getStatistics().getAverageRho())) {
      std::exit(-1);
    }
  }

  if (iT % vtkIter == 0 && vtkEnabled) {
    superLattice.executePostProcessors(stage::Evaluation{});
    superLattice.setProcessingContext(ProcessingContext::Evaluation);
    {
      SuperVTMwriter2D<T> vtkWriter("flow");
      SuperLatticePhysVelocity2D velocity(superLattice, converter);
      SuperLatticePhysPressure2D pressure(superLattice, converter);
      vtkWriter.addFunctor(pressure);
      vtkWriter.addFunctor(velocity);
      vtkWriter.write(iT);
    }
    {
      SuperVTMwriter2D<T> vtkWriter("structure");
      SuperLatticeField2D<T,DESCRIPTOR,descriptors::POROSITY> porosity(superLattice);
      SuperLatticeField2D<T,DESCRIPTOR,fields::fsi::ELEMENT_TAG> tag(superLattice);
      SuperLatticeField2D<T,DESCRIPTOR,fields::fsi::ELEMENT_FORCE> force(superLattice);
      SuperLatticeField2D<T,DESCRIPTOR,fields::fsi::ELEMENT_TORQUE> torque(superLattice);
      vtkWriter.addFunctor(porosity);
      vtkWriter.addFunctor(tag);
      vtkWriter.addFunctor(force);
      vtkWriter.addFunctor(torque);
      vtkWriter.write(iT);
    }
  }
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);
  singleton::directories().setOutputDir( "./tmp/" );

  OstreamManager clout(std::cout, "main");
  OstreamManager clresult(std::cout, "result");

  CLIreader args(argc, argv);
  const int N = args.getValueOrFallback("--resolution", 40);
  const bool vtkEnabled = !args.contains("--no-vtk");

  const T deltaX = 0.02 / N;
  const T deltaT = (9.6e-5 * (T{40}/N));

  const UnitConverter<T, DESCRIPTOR> converter(
    (T)   deltaX,           // deltaX
    (T)   deltaT,           // deltaT
    (T)   channelH,         // charPhysLength: reference length of simulation geometry
    (T)   0.04,             // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   4e-6,             // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1090              // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  auto domainI = makeDomainI();
  std::shared_ptr<IndicatorF2D<T>> paddedDomainI(new IndicatorLayer2D<T>(*domainI, converter.getPhysDeltaX()));

  CuboidDecomposition2D<T> cuboidDecomposition(*paddedDomainI, converter.getPhysDeltaX(), singleton::mpi().getSize());
  BlockLoadBalancer<T> loadBalancer(cuboidDecomposition);

  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer);
  prepareGeometry(converter, superGeometry);

  SuperLattice<T,DESCRIPTOR> superLattice(superGeometry);
  prepareLattice(superLattice, converter, superGeometry);

  SuperPorousElementEmbeddingO fsiEmbeddingO(superLattice);
  fsiEmbeddingO.registerElementType<CuboidPorosityF>(superGeometry.getMaterialIndicator(2));

  T valveW = 0.5*converter.getPhysDeltaX();

  T angle = 0;
  T angularU = 0.0;
  Vector lower{channelH, channelH+T{0.5}*valveW};

  auto valveParams = makeParametersD<T,descriptors::SPATIAL_DESCRIPTOR<2>>(
    fields::fsi::ELEMENT_TAG{},      1,
    fields::fsi::ELEMENT_LOWER{},    lower,
    fields::fsi::ELEMENT_UPPER{},    lower + Vector{valveL,valveW},
    fields::fsi::ELEMENT_PIVOT{},    lower + Vector{valveW/2,valveW/2},
    fields::fsi::ELEMENT_ROTATION{}, angle
  );
  fsiEmbeddingO.add(valveParams);

  SuperPorousElementReductionO<T,DESCRIPTOR,fields::fsi::ELEMENT_FORCE,fields::fsi::ELEMENT_TORQUE> fsiReductionO(
    superLattice,
    superGeometry.getMaterialIndicator(2));
  fsiReductionO.resize(1);
  fsiReductionO.addCollectionO(meta::id<CollectPorousBoundaryForceAndTorqueO>{});

  superLattice.setParameter<fields::converter::PHYS_LENGTH>(
    1/converter.getConversionFactorLength());
  superLattice.setParameter<fields::converter::PHYS_VELOCITY>(
    1/converter.getConversionFactorVelocity());
  superLattice.setParameter<fields::converter::PHYS_DELTA_X>(
    converter.getPhysDeltaX());

  superLattice.setProcessingContext(ProcessingContext::Simulation);
  fsiEmbeddingO.initialize();

  CSV<T> valveWriter("valve_n" + std::to_string(N), ';', {"time", "angle"}, ".csv");

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT),
                       superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
    setBoundaryValues(converter, superLattice, superGeometry, iT);
    superLattice.collideAndStream();
    getResults(superLattice, converter, superGeometry, timer, iT, vtkEnabled);

    fsiReductionO.apply();

    if (fsiReductionO.rankDoesFSI()) {
      for (unsigned iElement=1; iElement <= fsiReductionO.getElementCount(); ++iElement) {
        // Divide by deltaX due to 3D mass conversion factor choice in UnitConverter
        auto f = converter.getPhysForce(fsiReductionO.getField<fields::fsi::ELEMENT_FORCE>(iElement))
               / converter.getPhysDeltaX();
        auto t = converter.getPhysTorque(fsiReductionO.getField<fields::fsi::ELEMENT_TORQUE>(iElement))
               / converter.getPhysDeltaX();

        auto angularAccel = t / valveMomentaOfInertia;
        if (iT >= converter.getLatticeTime(T_p)) {
          angle += angularU * converter.getPhysDeltaT();
          angularU += angularAccel * converter.getPhysDeltaT();
        }

        fsiEmbeddingO.setField<fields::fsi::ELEMENT_ROTATION>(iElement-1, angle);
        fsiEmbeddingO.setField<fields::fsi::ELEMENT_U_ROTATION>(iElement-1, angularU);

        if (iT % converter.getLatticeTime(0.01) == 0) {
          clresult << "t=" << converter.getPhysTime(iT) << "; force=" << f << "; torque=" << t << "; angle=" << angle << "; angularU=" << angularU << "; angularA=" << angularAccel << std::endl;
          valveWriter.writeDataFile(converter.getPhysTime(iT), T{angle});
        }
      }

      fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_ROTATION>>(ProcessingContext::Simulation);
      fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_U_ROTATION>>(ProcessingContext::Simulation);
    }

    fsiEmbeddingO.apply();
  }

  superLattice.setProcessingContext(ProcessingContext::Evaluation);

  timer.stop();
  timer.printSummary();
}

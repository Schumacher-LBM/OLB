/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2022 Florian Raichle
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
 * micromixer3d
 * Simulating adsorption in a static mixing reactor using an Euler-Euler approach.
 * The model is based on the linear driving force model and uses advection diffusion reaction lattices for particles,
 * solute and particle loading.
 *
 * Different isotherms and mass transfer models can be used.
 */

#include <olb.h>
//#ifndef OLB_PRECOMPILED   // Unless precompiled version is used
//#endif
#include <vector>
#include <iostream>
#include <fstream>
#include <set>

//#include "../isotherms.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;

typedef D3Q19<> NSDESCRIPTOR;
typedef D3Q7<VELOCITY> ADEDESCRIPTOR;

using bulkDynamicsNS = BGKdynamics<T, NSDESCRIPTOR>;
using bulkDynamicsAD = SourcedAdvectionDiffusionBGKdynamics<T, ADEDESCRIPTOR>;

const T charU = 0.2;   // charU
const T width = 0.0133;  // hydraulic diameter = 4 * surface / perimeter
const int N = 20;         // resolution of the hydraulic diameter
const T Re = 50;          // Reynolds number
const T Sc = 100;         // Schmidt number
const T particleRadius = 5e-05;  // particles radius
const T particleConcentration = 5;
const T particleDensity = 1700;

const T isoConstA = 45;
const T isoConstB = 0.5;
const T k_f = 0.;
const T c_0 = 1.;
const T D_s = 5E-11;

const T tau = 0.6125;

const T physTotalTime = 5;   // time for fluid simulation
const T physStartTime = .6;    // time to start fluid pulsation

const int latticeOverlap = 2;

const T stlSize = 10.0;

void prepareGeometry(UnitConverter<T, NSDESCRIPTOR> const& converterNS,
                     IndicatorF3D<T>& indicator, SuperGeometry<T,3>& superGeometry)
{

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2, indicator);
  superGeometry.rename(2, 1, 1);

// Returns the minimum phys position in each direction for material 2
  Vector<T, 3> minR = superGeometry.getStatistics().getMinPhysR(2);
  Vector<T, 3> maxR = superGeometry.getStatistics().getMaxPhysR(2);
  Vector<T, 3> centerR = superGeometry.getStatistics().getCenterPhysR(2);
  Vector<T, 3> extend = superGeometry.getStatistics().getPhysExtend(2);
  extend[2] = converterNS.getPhysDeltaX();

  // sets circle of both inflows and the outflow, with direction and radius
  IndicatorCircle3D<T> inflow1(minR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / T(3));
  IndicatorCircle3D<T> inflow2(maxR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / T(3));
  IndicatorCircle3D<T> outflow(minR[0], maxR[1], centerR[2],
                               0., 1., 0., (maxR[0] - minR[0]) / T(3));

  // sets cylinder on that in-/out-flow circles with length
  IndicatorCylinder3D<T> layerInflow1(inflow1, converterNS.getPhysDeltaX());
  IndicatorCylinder3D<T> layerInflow2(inflow2, converterNS.getPhysDeltaX());
  IndicatorCylinder3D<T> layerOutflow(outflow, converterNS.getPhysDeltaX());
  // renames all boundary voxels of material fromBcMat to toBcMat if two neighbour voxel
  // in the direction of the discrete normal are fluid voxel with material fluidM in the region
  // where the indicator function is fulfilled
  superGeometry.rename(2, 3, 1, layerInflow1); // layer of inflow1 gets mat = 3
  superGeometry.rename(2, 4, 1, layerInflow2); // layer of inflow2 gets mat = 4
  superGeometry.rename(2, 5, 1, layerOutflow); // layer of outflow gets mat = 5

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLatticeNS(
  SuperGeometry<T,3>& superGeometry,
  SuperLattice<T, NSDESCRIPTOR>& sLattice,
  const T omega)

{

  OstreamManager clout(std::cout, "prepareLatticeNS");
  clout << "Prepare NSE Lattice ..." << std::endl;

  // dynamics for fluid
  sLattice.defineDynamics<bulkDynamicsNS>(superGeometry.getMaterialIndicator({1, 3, 4, 5}));
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  // boundary conditions for fluid

  // inlet
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);
  // outlet
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 5);


  // initialisation for fluid
  AnalyticalConst3D<T, T> rho1(1.);
  AnalyticalConst3D<T, T> u0(0., 0., 0.);

  sLattice.defineRhoU(superGeometry.getMaterialIndicator({1, 2, 3, 4, 5}), rho1, u0);
  sLattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 2, 3, 4, 5}), rho1, u0);

  // Make the lattice ready for simulation
  sLattice.setParameter<descriptors::OMEGA>(omega);
  sLattice.initialize();

  {
    auto &communicator = sLattice.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(sLattice.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare NSE Lattice ... OK" << std::endl;
}

void prepareLatticeAD(
  SuperGeometry<T,3>& superGeometry,
  SuperLattice<T, ADEDESCRIPTOR>*& sLatticeAD,
  const T omegaAD, int inlet, int noInlet, T rhoInlet)

{

  OstreamManager clout(std::cout, "prepareLatticeAD");
  clout << "Prepare ADE Lattice ..." << std::endl;

  // dynamics for ADE
  sLatticeAD->defineDynamics<bulkDynamicsAD>(superGeometry.getMaterialIndicator({1, 5}));
  boundary::set<boundary::BounceBack>(*sLatticeAD, superGeometry, 2);
  sLatticeAD->defineDynamics<bulkDynamicsAD>(superGeometry, inlet);
  boundary::set<boundary::BounceBack>(*sLatticeAD, superGeometry, noInlet);

  // boundary for ADE
  boundary::set<boundary::AdvectionDiffusionDirichlet>(*sLatticeAD, superGeometry, inlet);
  setZeroGradientBoundary<T,ADEDESCRIPTOR>(*sLatticeAD, superGeometry, 5);

  // initialisation for fluid
  AnalyticalConst3D<T, T> rhoI(rhoInlet);
  AnalyticalConst3D<T, T> u0(0., 0., 0.);
  AnalyticalConst3D<T, T> rhoSmall(0);

  sLatticeAD->defineRhoU(superGeometry.getMaterialIndicator({1, 2, 5}), rhoSmall, u0);
  sLatticeAD->defineRhoU(superGeometry, inlet, rhoI, u0);
  sLatticeAD->defineRhoU(superGeometry, noInlet, rhoSmall, u0);

  sLatticeAD->iniEquilibrium(superGeometry.getMaterialIndicator({1, 2, 5}), rhoSmall, u0);
  sLatticeAD->iniEquilibrium(superGeometry, inlet, rhoI, u0);
  sLatticeAD->iniEquilibrium(superGeometry, noInlet, rhoSmall, u0);

  // Make the lattice ready for simulation
  sLatticeAD->setParameter<descriptors::OMEGA>(omegaAD);
  sLatticeAD->initialize();

  clout << "Prepare ADE Lattice ... OK" << std::endl;
}

void setBoundaryValues(
  SuperGeometry<T,3>& superGeometry,
  SuperLattice<T, NSDESCRIPTOR>& sLattice,
  UnitConverter<T, NSDESCRIPTOR> const& converterNS,
  std::size_t iT)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  std::vector < T > maxVelocity(3, T());
  const T distanceToBoundary = converterNS.getPhysDeltaX() / T(2);
  const T latticeVelNS = converterNS.getLatticeVelocity(converterNS.getCharPhysVelocity());
  const size_t itStartTime = converterNS.getLatticeTime(physStartTime);

  if (iT <= itStartTime && iT % 50 == 0) {

    SinusStartScale<T, size_t> startScale(itStartTime, T(1)); //1+-amplitude
    size_t help[1] = { iT };
    T frac[3] = { T() };
    startScale(frac, help);

    // set lattice velocity on boundary
    maxVelocity[1] = latticeVelNS * frac[0];

    RectanglePoiseuille3D<T> u5(superGeometry, 5, maxVelocity,
                                distanceToBoundary, distanceToBoundary, distanceToBoundary);
    sLattice.defineU(superGeometry, 5, u5);
    sLattice.setProcessingContext<olb::Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults(
  SuperLattice<T, NSDESCRIPTOR>& sLattice,
  UnitConverter<T, NSDESCRIPTOR>& converterNS,
  std::vector<SuperLattice<T, ADEDESCRIPTOR>*> partners,
  size_t iT, SuperGeometry<T,3>& superGeometry,
  util::Timer<double>& timer)
{
  OstreamManager clout(std::cout, "getResults");

  const int itTotalTime = converterNS.getLatticeTime(physTotalTime);

  // output in beginning
  const int iTperiodConsole = itTotalTime / 100;   // Console output
  const int iTperiodVTK = itTotalTime / 100;   // Writes the vtk files

  SuperVTMwriter3D<T> vtmWriter( "microMixer3d" );

  if (iT == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank(sLattice);

    vtmWriter.write(cuboid);
    vtmWriter.write(rank);

    // have to be called before calling write(int iT=0), since it creates
    // the master pvd file, where all vti are linked!
    vtmWriter.createMasterFile();
  }


  // console output
  if (iT % iTperiodConsole == 0) {
    timer.update(iT);
    timer.printStep();

    // output for latticeStatistics
    sLattice.getStatistics().print(iT, converterNS.getPhysTime(iT));
    partners[0]->getStatistics().print(iT, converterNS.getPhysTime(iT));
  }

  // vtk and gif output
  if (iT % iTperiodVTK == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    partners[0]->setProcessingContext(ProcessingContext::Evaluation);
    partners[1]->setProcessingContext(ProcessingContext::Evaluation);
    partners[2]->setProcessingContext(ProcessingContext::Evaluation);
    SuperGeometryF<T,3> materials(superGeometry);
    SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocityNS(sLattice, converterNS);
    SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure(sLattice, converterNS);
    SuperLatticeDensity3D<T, ADEDESCRIPTOR> adsorptive(*partners[0]);
    SuperLatticeDensity3D<T, ADEDESCRIPTOR> loading(*partners[2]);
    SuperLatticeDensity3D<T, ADEDESCRIPTOR> soluteConcentration(*partners[1]);
    vtmWriter.addFunctor(velocityNS);
    vtmWriter.addFunctor(pressure);
    vtmWriter.addFunctor(adsorptive, "particle concentration");
    vtmWriter.addFunctor(loading, "loading");
    vtmWriter.addFunctor(soluteConcentration, "solute concentration");
    vtmWriter.addFunctor(materials);
    vtmWriter.write(iT);
  }
}


int main(int argc, char* argv[])
{

  /// === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  std::string configFile = "microMixer3D.xml";

  UnitConverterFromResolutionAndRelaxationTime<T, NSDESCRIPTOR> converterNS(
      (T) N,                           // resolution: number of voxels per charPhysL
      (T) tau,                          // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T) width,                     // charPhysLength: reference channelLength of simulation geometry
      (T) charU,                        // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T) charU * width /Re,            // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T) 1000                         // physDensity: physical density in __kg / m^3__
  );

  AdsorptionConverter<T, ADEDESCRIPTOR> converterADE(
      (T)   width/N,     //physDeltaX
      (T)   converterNS.getPhysDeltaT(),     //physDeltaT
      (T)   width,          //charPhysLength
      (T)   charU,          //charPhysVelocity
      (T)   charU*width / (Sc*Re),
      (T)   1.,
      (T)   1.,
      (T)   charU*width / Re
      );

  clout << "navier-stokes converter " << std::endl;
  converterNS.write("fluid");
  converterNS.print();
  converterADE.write("adsorption");
  converterADE.print();

  /// compute relaxation parameter to solve the advection-diffusion
  /// equation in the lattice Boltzmann context
  T omegaAD = converterADE.getLatticeRelaxationFrequency();
  clout << "omegaAD = " << omegaAD << std::endl;
  clout << "tauAD = " << 1./omegaAD << std::endl;
  /// === 2nd Step: Prepare Geometry ===

  // Instantiation of an empty cuboidDecomposition
  const int noOfCuboids = util::max(16, 4 * singleton::mpi().getSize());

  STLreader<T> stlReader("microMixer3d_small.stl",
                         converterNS.getPhysDeltaX(), stlSize);

  IndicatorLayer3D<T> extendedDomain( stlReader,
                                      converterNS.getPhysDeltaX() );

  CuboidDecomposition3D<T> cuboidDecomposition( extendedDomain,
                                      converterNS.getPhysDeltaX(), noOfCuboids );

  // Instantiation of an empty loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // Default instantiation of superGeometry
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, latticeOverlap);

  prepareGeometry(converterNS, stlReader, superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, NSDESCRIPTOR> sLattice(superGeometry);
  SuperLattice<T, ADEDESCRIPTOR> sLatticeAD(superGeometry);
  SuperLattice<T, ADEDESCRIPTOR> CADLattice(superGeometry);
  SuperLattice<T, ADEDESCRIPTOR> QADLattice(superGeometry);

  std::vector<SuperLattice<T, ADEDESCRIPTOR>*> partners;
  partners.emplace_back(&sLatticeAD);
  partners.emplace_back(&CADLattice);
  partners.emplace_back(&QADLattice);

  int inlet[3] = {3, 4, 3};
  int noInlet[3] = {4, 3, 4};
  T rhoInlet[3] = {1.0, 1.0, 0.};

  prepareLatticeNS(superGeometry, sLattice, converterNS.getLatticeRelaxationFrequency());

  for(int i = 0; i<3; i++){
    prepareLatticeAD(superGeometry, partners[i], omegaAD, inlet[i], noInlet[i], rhoInlet[i]);
    }

  SuperLatticeCoupling coupling(
    AdsorptionFullCoupling3D<AdsorptionReaction<Isotherm::LangmuirIsotherm>,
                             ade_forces::AdvDiffDragForce3D>{},
    names::NavierStokes{}, sLattice,
    names::Concentration0{}, sLatticeAD,
    names::Concentration1{}, CADLattice,
    names::Concentration2{}, QADLattice );

  // Setting Isotherm Parameters
  Isotherm::LangmuirIsotherm::setParameters<T>(isoConstA, isoConstB, coupling);

  // Setting Adsorption Reaction Parameters
  coupling.template setParameter<AdsorptionReaction<Isotherm::LangmuirIsotherm>::K_F>(k_f);
  coupling.template setParameter<AdsorptionReaction<Isotherm::LangmuirIsotherm>::D_S>(D_s);
  coupling.template setParameter<AdsorptionReaction<Isotherm::LangmuirIsotherm>::C_0>(c_0);
  coupling.template setParameter<AdsorptionReaction<Isotherm::LangmuirIsotherm>::R_P>(particleRadius);

  // Compute the interaction parameters
  AdsorptionReaction<Isotherm::LangmuirIsotherm>::computeParameters<T>(coupling, converterADE);

  AdsorptionReaction<Isotherm::LangmuirIsotherm>::print<T>(clout, coupling, converterADE);

  // Compute the drag force parameters
  ade_forces::AdvDiffDragForce3D::computeParametersFromRhoAndRadius<T>(particleDensity, particleRadius, coupling, converterNS);

  /// === 4th Step: Fluid Main Loop with Timer ===
  util::Timer<double> timer(converterNS.getLatticeTime(physTotalTime),
                      superGeometry.getStatistics().getNvoxel());

  timer.start();

  for (size_t iT = 0; iT < converterNS.getLatticeTime(physTotalTime); ++iT) {
    setBoundaryValues(superGeometry, sLattice, converterNS, iT);

    coupling.execute();

    getResults(sLattice, converterNS, partners, iT, superGeometry, timer);

    for (int i = 0; i<3; i++) {
      partners[i]->collideAndStream();
    }

    sLattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();
}

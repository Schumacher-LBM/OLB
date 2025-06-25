/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Adrian Kummerlaender, Yuji (Sam) Shimojima,
 *                     Jonathan Jeppener-Haltenhoff, Marc Haußmann,
 *                     Mathias J. Krause
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

/* channel3d.cpp:
 * This example examines a wall-bounded flow at high Reynolds numbers.
 * The dynamics follow the RLB collision operator of 3rd order and Smagorinsky-Lilly turbulence model.
 * The near wall region is modelled by a wall function.
 *
 * The example shows the usage of turbulent wall model with different parameters.
 * All WM paramertes are integrated into HLBM interface.
 * As showcase for stability improvement HRR collision is possible.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::util;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D3Q19<FORCE,TENSOR,VELOCITY,POROSITY,VELOCITY2>;

#define HRR_COLLISION

// Parameters for the simulation setup
const int N = 40;
const T physRefL = 1.0;     // half channel height in meters
const T lx = 2. * std::numbers::pi_v<T> * physRefL;  // streamwise length in meters
const T ly = 2. * std::numbers::pi_v<T> * physRefL/10.;  // spanwise length in meters
const T lz = 2. * physRefL;       // wall-normal length in meters
#ifdef HRR_COLLISION
const T hybridConst = 0.99;   // very strong influence of numerical diffusion. 0.99 meand 1% FDM strain rate tensor, which is already enough!
#endif

// Choose friction reynolds number ReTau
#define Case_ReTau_1000
//#define Case_ReTau_2000

/// Wallfunction parameters
//  Used method for density reconstruction
//  0: use local density
//  1: extrapolation (Guo)
//  2: constant (rho = 1.)
const int rhoMethod = 0;
//  Used method for non-equilibrium population reconstruction
//  0: extrapolation NEQ (Guo Zhaoli)
//  1: first or second order finite differnce (Malaspinas)
//  2: equilibrium scheme (no fNeq)
const int fNeqMethod = 1;
//  Used wall profile
//  0: power law profile
//  1: Spalding profile
const int wallFunctionProfile = 1;
// check if descriptor with body force is used
const bool bodyForce = true;
// interpolate sampling velocity along given normal between lattice voxels
const bool interpolateSampleVelocity = true;
// use van Driest damping function for turbulent viscosity in boundary cell
const bool useVanDriest = true;
//  distance from cell to real wall in lattice units if no geometry indicator is given as input
const T latticeWallDistance = 0.5;
//  distance from cell to velocity sampling point in lattice units
const T samplingCellDistance = 3.5;

// Reynolds number based on the friction velocity
#if defined (Case_ReTau_1000)
T ReTau = 1000.512;
#elif defined (Case_ReTau_2000)
T ReTau = 1999.756;
#endif
// Characteristic physical kinematic viscosity from DNS Data
// http://turbulence.ices.utexas.edu/channel2015/data/LM_Channel_1000_mean_prof.dat
#if defined (Case_ReTau_1000)
T charPhysNu = 5./100000.;
#elif defined (Case_ReTau_2000)
T charPhysNu = 2.3/100000.;
#endif

// physical simulated length adapted for lattice distance to boundary in meters
const T adaptedPhysSimulatedLength = 2 * physRefL / ( 1. - 2./N*(1.-latticeWallDistance) );
// Characteristic physical mean bulk velocity from Dean correlations in meters - Malaspinas and Sagaut (2014)
const T charPhysU = ( util::pow((8.0/0.073), (4.0/7.0)) * util::pow((T)ReTau, (8.0/7.0)) ) * charPhysNu / (2. * physRefL);
// Time of the simulation in seconds
const T charPhysT = physRefL / (ReTau * charPhysNu / physRefL);

const T physConvergeTime = 40. * charPhysT;  // time until until statistics sampling in seconds
const T physStatisticsTime = 150. * charPhysT ;  // statistics sampling time in seconds
const T maxPhysT = physConvergeTime + physStatisticsTime; // max. simulation time in seconds
const T statisticsSave = 1./25.;    // time between statistics samples in seconds

// seed the rng with time if SEED_WITH_TIME is set, otherwise just use a fixed seed.
#if defined SEED_WITH_TIME
#include <chrono>
auto seed = std::chrono::system_clock().now().time_since_epoch().count();
std::default_random_engine generator(seed);
#else
std::default_random_engine generator(0x1337533DAAAAAAAA);
#endif

struct CellSpatialCalculationVelocityField {
  static constexpr OperatorScope scope = OperatorScope::PerCell;
  int getPriority() const { return 2; }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    //every value is the lattice(non-dimensional) unit
    using DESCRIPTOR = typename CELL::descriptor_t;
    using V          = typename CELL::value_t;
    Vector<V, DESCRIPTOR::d> lattice_u {};
    cell.computeU(lattice_u.data());
    cell.template setField<descriptors::VELOCITY2>(lattice_u);
  }
};

template <typename T, typename S>
class Channel3D : public AnalyticalF3D<T,S> {

protected:
  T turbulenceIntensity;
  T maxVelocity;
  T distanceToWall;
  T obst_z;
  T obst_r;
  T a;
  T b;

public:
  Channel3D(UnitConverter<T,DESCRIPTOR> const& converter, T frac) : AnalyticalF3D<T,S>(3)
  {
    turbulenceIntensity = 0.05;
    maxVelocity = converter.getLatticeVelocity(converter.getCharPhysVelocity()*(8./7.)); // Centerline Velocity
    obst_r = physRefL;
    a = -1.;
    b = 1.;
  };

  bool operator()(T output[], const S input[])
  {
    std::uniform_real_distribution<BaseType<T>> distribution(a, b);
    T nRandom1 = distribution(generator);
    T nRandom2 = distribution(generator);
    T nRandom3 = distribution(generator);

    if( (util::abs(input[2] - obst_r)) < obst_r ){

    T u_calc = maxVelocity*util::pow(((obst_r-util::abs(input[2] - obst_r))/obst_r), 1./7.);

    output[0] = turbulenceIntensity*nRandom1*maxVelocity + u_calc;
    output[1] = turbulenceIntensity*nRandom2*maxVelocity;
    output[2] = turbulenceIntensity*nRandom3*maxVelocity;
    }
    return true;
  };
};

void prepareGeometry(SuperGeometry<T,3>& superGeometry,
                     IndicatorF3D<T>& indicator,
                     const UnitConverter<T,DESCRIPTOR>& converter)
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2,indicator);
  superGeometry.rename(2,1,{0,0,1});
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  olb::Vector<T, 3> PhyMax = superGeometry.getStatistics().getMaxPhysR(2);
  olb::Vector<T, 3> PhyMin = superGeometry.getStatistics().getMinPhysR(2);
  clout << "Dimension of the channel in meters: x = " << PhyMax[0] - PhyMin[0];
  clout << " ; y = " << PhyMax[1] - PhyMin[1];
  clout << " ; z = " << PhyMax[2] - PhyMin[2] << std::endl;

  clout << "Prepare Geometry ... OK" << std::endl;
}

// set up initial conditions
void setInitialConditions(SuperLattice<T, DESCRIPTOR>& sLattice,
                          const UnitConverter<T,DESCRIPTOR>& converter,
                          SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "setInitialConditions");
  clout << "Set initial conditions ..." << std::endl;

  AnalyticalConst3D<T, T> rho(1.);
  AnalyticalConst3D<T, T> rho0(0.);
  AnalyticalConst3D<T, T> u0(0.,0.,0.);
  AnalyticalConst3D<T,T> tag1(1);
  AnalyticalConst3D<T,T> tag0(0);
  Channel3D<T, T> uSol(converter, 1.);

  sLattice.defineRhoU(superGeometry, 1, rho, uSol);
  sLattice.iniEquilibrium(superGeometry, 1, rho, uSol);

  sLattice.defineRhoU(superGeometry, 2, rho, uSol);
  sLattice.iniEquilibrium(superGeometry, 2, rho, uSol);

  sLattice.defineField<reduction::TAGS_U>(superGeometry, 1, tag1);
  sLattice.defineField<reduction::TAGS_U>(superGeometry.getMaterialIndicator({0,2}), tag0);
  sLattice.defineField<descriptors::VELOCITY2>(superGeometry.getMaterialIndicator({0,1,2}), uSol);
  sLattice.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({0,1,2}), u0);
  sLattice.defineField<FORCE>(superGeometry.getMaterialIndicator({0,1,2}), u0);
  sLattice.defineField<descriptors::POROSITY>(superGeometry.getMaterialIndicator({0,2}), rho0);
  sLattice.defineField<descriptors::POROSITY>(superGeometry, 1, rho);
  sLattice.addPostProcessor<stage::PostStream>(superGeometry.getMaterialIndicator({1,2}),
                                               meta::id<CellSpatialCalculationVelocityField>{});

  clout << "Set initial conditions ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    const UnitConverter<T,DESCRIPTOR>& converter,
                    SuperGeometry<T,3>& superGeometry,
                    WallModelParameters<T>& wallModelParameters)
{
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  /// Material=1 -->bulk dynamics
  setTurbulentWallModelDynamics(sLattice, superGeometry, 1, wallModelParameters);
  sLattice.defineDynamics<BounceBack>(superGeometry, 2);
#ifdef HRR_COLLISION
  sLattice.addPostProcessor<stage::PostStream>(superGeometry.getMaterialIndicator({1}),
                                               meta::id<FDMstrainRateTensorPostProcessor>{});
  AnalyticalConst3D<T, T> hybrid(hybridConst);
  sLattice.defineField<collision::HYBRID>(superGeometry, 1, hybrid);
#endif

  /// Material = 2 --> boundary node + wallfunction

  /// === Set Initial Conditions == ///
  setTurbulentWallModel(sLattice, superGeometry, 2, wallModelParameters);
  setInitialConditions(sLattice, converter, superGeometry);

  sLattice.setParameter<descriptors::OMEGA>( converter.getLatticeRelaxationFrequency() );
  sLattice.setParameter<collision::LES::SMAGORINSKY>(T(0.12));

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Computes the pressure drop between the voxels before and after the cylinder
void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, size_t iT,
                SuperGeometry<T,3>& superGeometry, Timer<double>& timer,
                SuperLatticeTimeAveragedF3D<T>& sAveragedVel,
                SuperLatticeTimeAveragedF3D<T>& sAveragedPress)
{
  OstreamManager clout(std::cout, "getResults");

  const T checkstatistics = converter.getPhysTime(100);//(T)maxPhysT/200.;

  SuperVTMwriter3D<T> vtmWriter("channel3d");
  SuperLatticePhysVelocity3D velocity(sLattice, converter);
  SuperLatticePhysPressure3D pressure(sLattice, converter);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);
  vtmWriter.addFunctor(sAveragedVel);
  vtmWriter.addFunctor(sAveragedPress);


  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D cuboid(sLattice);
    SuperLatticeRank3D rank(sLattice);

    vtmWriter.write(cuboid);
    vtmWriter.write(rank);

    vtmWriter.createMasterFile();
  }

  // Writes output on the console
  if (iT % converter.getLatticeTime(checkstatistics) == 0) {
    // Timer console output
    timer.update(iT);
    timer.printStep(2);
    // Lattice statistics console output
    sLattice.getStatistics().print(iT, iT * converter.getPhysDeltaT());
    clout << "Max. physical velocity(m/s): " << converter.getPhysVelocity(sLattice.getStatistics().getMaxU()) << std::endl;
    clout << "Max u+:" << converter.getPhysVelocity(sLattice.getStatistics().getMaxU())
          / (ReTau * charPhysNu / (converter.getCharPhysLength()/2.)) << std::endl;
  }

  if (iT%converter.getLatticeTime(statisticsSave) == 0 && iT > converter.getLatticeTime(physConvergeTime)) {
    // Add ensemble to temporal averaged velocity
    sLattice.communicate();
    sAveragedVel.addEnsemble();
    sAveragedPress.addEnsemble();
  }

  if (iT%converter.getLatticeTime(checkstatistics)/*converter.getLatticeTime(maxPhysT/20)*/==0 || iT==converter.getLatticeTime(maxPhysT)-1) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    // Writes the vtk files
    vtmWriter.write(iT);
  }
}

int main(int argc, char* argv[])
{

  /// === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> const converter(
    int {N},                  // resolution: number of voxels per charPhysL
    (T)   0.50025,            // relaxation time
    (T)   adaptedPhysSimulatedLength, // charPhysLength: reference length of simulation geometry
    (T)   1.0,                // charPhysVelocity: mean bulk velocity in __m / s__
    (T)   charPhysNu,           // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                 // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "Converge time(s): " << physConvergeTime << std::endl;
  clout << "Lattice converge time: " << converter.getLatticeTime(physConvergeTime) << std::endl;
  clout << "Max. Phys. simulation time(s): " << maxPhysT << std::endl;
  clout << "Max. Lattice simulation time: " << converter.getLatticeTime(maxPhysT) << std::endl;
  clout << "Frequency Statistics Save(Hz): " << 1./statisticsSave << std::endl;
  clout << "Statistics save period(s): " << statisticsSave << std::endl;
  clout << "Lattice statistics save period: " << converter.getLatticeTime(statisticsSave) << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "Channel height(m): " << adaptedPhysSimulatedLength << std::endl;
  clout << "y+ value: " << (ReTau * converter.getPhysViscosity() / (physRefL)) * ((2. / T(N + 2 * latticeWallDistance)) * latticeWallDistance)/ converter.getPhysViscosity() << std::endl;
  clout << "y+ value spacing: " << (ReTau * converter.getPhysViscosity() / (physRefL)) * (converter.getPhysDeltaX()) / converter.getPhysViscosity() << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif

  T dx = converter.getPhysDeltaX();
  Vector<T,3> extend(lx, ly, adaptedPhysSimulatedLength);
  Vector<T,3> origin(0., 0., -(1.-latticeWallDistance)*dx);
  IndicatorCuboid3D<T> cuboid(extend, origin );

  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, converter.getPhysDeltaX(), noOfCuboids );

  cuboidDecomposition.setPeriodicity({true, true, false});

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, 4);

  prepareGeometry(superGeometry, cuboid, converter);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);


  WallModelParameters<T> wallModelParameters;
  wallModelParameters.bodyForce = bodyForce;
  wallModelParameters.rhoMethod = rhoMethod;
  wallModelParameters.fNeqMethod = fNeqMethod;
  wallModelParameters.samplingCellDistance = samplingCellDistance;
  wallModelParameters.interpolateSampleVelocity = interpolateSampleVelocity;
  wallModelParameters.useVanDriest = useVanDriest;
  wallModelParameters.wallFunctionProfile = wallFunctionProfile;
  wallModelParameters.latticeWallDistance = latticeWallDistance;

  prepareLattice(sLattice, converter, superGeometry, wallModelParameters);

  //forcing of the channel
  SuperLatticeFieldReductionO<T, DESCRIPTOR, descriptors::VELOCITY2, reduction::SumO,
                              reduction::checkBulkTag<reduction::TAGS_U>> sumProcessor(sLattice);
  SuperLatticeCoupling coupling(
    TurbulentChannelForce<T>{},
    names::NavierStokes{}, sLattice);
  coupling.setParameter<TurbulentChannelForce<T>::CHAR_LATTICE_U>(converter.getCharLatticeVelocity());
  coupling.setParameter<TurbulentChannelForce<T>::CHAR_LATTICE_L>(converter.getLatticeLength(adaptedPhysSimulatedLength));
  coupling.setParameter<TurbulentChannelForce<T>::LATTICE_UTAU>((T)converter.getLatticeVelocity(ReTau * converter.getPhysViscosity()/physRefL));
  coupling.setParameter<TurbulentChannelForce<T>::LATTICE_CHANNEL_VOLUME>(converter.getLatticeLength(lz) *converter.getLatticeLength(ly) *converter.getLatticeLength(lx));
  coupling.setParameter<TurbulentChannelForce<T>::LATTICE_U_SUM>({converter.getCharLatticeVelocity(),0.,0.});
  coupling.restrictTo(superGeometry.getMaterialIndicator(1));
  coupling.execute();

  /// === 5th Step: Definition of turbulent Statistics Objects ===

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> sVel(sLattice, converter);
  SuperLatticeTimeAveragedF3D<T> sAveragedVel(sVel);

  SuperLatticePhysPressure3D<T, DESCRIPTOR> sPress(sLattice, converter);
  SuperLatticeTimeAveragedF3D<T> sAveragedPress(sPress);

  /// === 4th Step: Main Loop with Timer ===
  Timer<double> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (size_t iT=0; iT<converter.getLatticeTime(maxPhysT); ++iT) {
    /// === 6th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer, sAveragedVel, sAveragedPress);
    coupling.setParameter<TurbulentChannelForce<T>::LATTICE_U_SUM>(sumProcessor.compute());
    coupling.execute();
    /// === 7th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
  }
  timer.stop();
  timer.printSummary();

  return 0;
}

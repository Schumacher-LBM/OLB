/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2023 Stephan Simonis, Mathias J. Krause, Patrick Nathan,
 *                     Alejandro C. Barreto, Marc Haußmann, Liam Sauterleute
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

/* tgv3d.cpp:
 * The Taylor-Green-Vortex (TGV) is one of the simplest configuration,
 * where you can investigate the generation of small structures and the resulting turbulence.
 * The 2pi periodic box domain and the single mode initial conditions contribute to the simplicity.
 * In consequence, the TGV is a common benchmark case for
 * Direct Numerical Simulations (DNS) and Large Eddy Simulations (LES).
 *
 * This example shows the usage and the effects of different subgrid scale turbulence models.
 * The molecular dissipation rate, the eddy dissipation rate and
 * the effective dissipation rate are calculated and plotted over the simulation time.
 * This results can be compared with a published DNS solution, e.g.
 * Brachet, Marc E., et al. "Small-scale structure of the Taylor–Green vortex."
 * Journal of Fluid Mechanics 130 (1983): 411-452.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;

// Choose turbulence model or collision scheme
//#define RLB
#define _SMAGORINSKY
//#define WALE
//#define ConsistentStrainSmagorinsky
//#define ShearSmagorinsky
//#define Krause
//#define DNS
//#define KBC
//#define K_EPSILON

#define finiteDiff //for N<256

#ifdef ShearSmagorinsky
typedef D3Q19<AV_SHEAR> DESCRIPTOR;
#elif defined (WALE)
typedef D3Q19<EFFECTIVE_OMEGA,VELO_GRAD> DESCRIPTOR;
#elif defined(KBC)
typedef D3Q19<GAMMA> DESCRIPTOR;
#elif defined(K_EPSILON)
using DESCRIPTOR = D3Q19<OMEGA>;
using KEDESCRIPTOR = D3Q7<VELOCITY, SOURCE, OMEGA, K, EPSILON>;
#else
typedef D3Q19<> DESCRIPTOR;
#endif

#if defined(RLB)
using BulkDynamics = RLBdynamics<T,DESCRIPTOR>;
#elif defined(DNS)
using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;
#elif defined(WALE)
using BulkDynamics = WALEBGKdynamics<T,DESCRIPTOR>;
#elif defined(ShearSmagorinsky)
using BulkDynamics = ShearSmagorinskyBGKdynamics<T,DESCRIPTOR>;
#elif defined(Krause)
using BulkDynamics = KrauseBGKdynamics<T,DESCRIPTOR>;
#elif defined(ConsistentStrainSmagorinsky)
using BulkDynamics = ConStrainSmagorinskyBGKdynamics<T,DESCRIPTOR>;
#elif defined(KBC)
using BulkDynamics = KBCdynamics<T, DESCRIPTOR>;
#elif defined(K_EPSILON)
using BulkDynamics = dynamics::Tuple<T,DESCRIPTOR,momenta::BulkTuple,equilibria::SecondOrder,collision::ParameterFromCell<descriptors::OMEGA,collision::BGK>>;
using KEBulkDynamics = SourcedLimitedAdvectionDiffusionBGKdynamics<T,KEDESCRIPTOR>;
#else
using BulkDynamics = SmagorinskyBGKdynamics<T,DESCRIPTOR>;
#endif

// Global constants
const T pi = 4.0 * util::atan(1.0);
const T volume = util::pow(2. * pi, 3.); // volume of the 2pi periodic box


// Parameters for the simulation setup
const T maxPhysT = 10;    // max. simulation time in s, SI unit
int N = 64;              // resolution of the model
T Re = 800;               // defined as 1/kinematic viscosity
T smagoConst = 0.1;       // Smagorisky Constant, for ConsistentStrainSmagorinsky smagoConst = 0.033
T vtkSave = 0.25;         // time interval in s for vtk and gnuplot output
const T intensity = 0.04; // turbulent intensity for K_EPSILON dynamics


bool plotDNS = true;      //available for Re=800, Re=1600, Re=3000 (maxPhysT<=10)
std::vector<std::vector<T>> values_DNS;

template <typename T, typename _DESCRIPTOR>
class Tgv3D : public AnalyticalF3D<T,T> {

protected:
  T u0;

// initial solution of the TGV
public:
  Tgv3D(UnitConverter<T,_DESCRIPTOR> const& converter, T frac) : AnalyticalF3D<T,T>(3)
  {
    u0 = converter.getCharLatticeVelocity();
  };

  bool operator()(T output[], const T input[]) override
  {
    const T x = input[0];
    const T y = input[1];
    const T z = input[2];

    output[0] = u0 * util::sin(x) * util::cos(y) * util::cos(z);
    output[1] = -u0 * util::cos(x) * util::sin(y) * util::cos(z);
    output[2] = 0;

    return true;
  };
};

#ifdef K_EPSILON
template <typename T, typename _DESCRIPTOR>
class TgvK3D : public AnalyticalF3D<T,T> {

protected:
  T u0;

// initial values of k for the the TGV
public:
  TgvK3D(UnitConverter<T,_DESCRIPTOR> const& converter, T frac) : AnalyticalF3D<T,T>(1)
  {
    u0 = converter.getCharLatticeVelocity();
  };

  bool operator()(T output[], const T input[]) override
  {
    output[0] = 3./2. * intensity * intensity * u0 * u0;
    return true;
  };
};

template <typename T, typename _DESCRIPTOR>
class TgvE3D : public AnalyticalF3D<T,T> {

protected:
  T u0, length, dX;

// initial values of epsilon for the TGV
public:
  TgvE3D(UnitConverter<T,_DESCRIPTOR> const& converter, T frac) : AnalyticalF3D<T,T>(1)
  {
    u0 = converter.getCharLatticeVelocity();
    length = converter.getCharPhysLength();
    dX = converter.getPhysDeltaX();
  };

  bool operator()(T output[], const T input[]) override
  {
    T k = 3./2. * intensity * intensity * u0 * u0;
    output[0] = util::pow(0.09, 3./4.)*util::pow(k, 3./2.)/(length/dX)/0.04;

    return true;
  };
};
#endif



void prepareGeometry(SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,1);
  superGeometry.communicate();

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  /// Material=1 -->bulk dynamics
  const T omega = converter.getLatticeRelaxationFrequency();
  sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);
  sLattice.setParameter<descriptors::OMEGA>(omega);
  #if !defined(RLB) && !defined(DNS) && !defined(KBC) && !defined(K_EPSILON)
  sLattice.setParameter<collision::LES::SMAGORINSKY>(smagoConst);
  #endif

  sLattice.initialize();
  clout << "Prepare Lattice ... OK" << std::endl;
}

#ifdef K_EPSILON
void prepareLatticeKE(SuperLattice<T,KEDESCRIPTOR>& sLatticeKE,
                    UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice KE ...";

  /// Material=1 -->bulk dynamics
  sLatticeKE.defineDynamics<KEBulkDynamics>(superGeometry, 1);

  clout << " OK" << std::endl;
}

void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& sLatticeRANS,
                       SuperLattice<T, KEDESCRIPTOR>& sLatticeK,
                       SuperLattice<T, KEDESCRIPTOR>& sLatticeE,
                       UnitConverter<T,DESCRIPTOR> const& converter,
                       SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout,"setBoundaryValues");

  const T omega = converter.getLatticeRelaxationFrequency();
  const T dX = converter.getPhysDeltaX();
  const T dT = converter.getPhysDeltaT();
  T convK = dX*dX/dT/dT;
  T convE = dX*dX/dT/dT/dT;
  T lMix = converter.getCharPhysLength()*0.04;
  T turbKinEnergy = util::pow(converter.getPhysViscosity()/0.1/(lMix),2.);

  AnalyticalConst3D<T,T> omega0(omega);
  AnalyticalConst3D<T,T> rho0(1e-8);
  AnalyticalConst3D<T,T> rhoK(turbKinEnergy/convK);
  AnalyticalConst3D<T,T> rhoE(0.09*util::pow(turbKinEnergy, 1.5)/lMix/0.1/convE);
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> u( velocity );

  AnalyticalConst3D<T,T> rho(1.);
  Tgv3D<T,DESCRIPTOR> uSol(converter, 1);
  TgvE3D<T,DESCRIPTOR> eSol(converter, 1);
  TgvK3D<T,DESCRIPTOR> kSol(converter, 1);

  sLatticeRANS.defineRhoU(superGeometry, 1, rho, uSol);
  sLatticeRANS.iniEquilibrium(superGeometry, 1, rho, uSol);
  sLatticeRANS.defineField<descriptors::OMEGA>(superGeometry.getMaterialIndicator(1), omega0);
  sLatticeRANS.defineField<descriptors::SCALAR>(superGeometry.getMaterialIndicator(1), rho0);

  sLatticeK.defineRhoU(superGeometry, 1, kSol, uSol);
  sLatticeK.iniEquilibrium(superGeometry, 1, kSol, uSol);
  sLatticeE.defineRhoU(superGeometry, 1, eSol, uSol);
  sLatticeE.iniEquilibrium(superGeometry, 1, eSol, uSol);

  sLatticeK.defineField<descriptors::SOURCE>(superGeometry.getMaterialIndicator(1), rho0);
  sLatticeK.defineField<descriptors::K>(superGeometry.getMaterialIndicator(1), kSol);
  sLatticeK.defineField<descriptors::OMEGA>(superGeometry.getMaterialIndicator(1), omega0);
  sLatticeK.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator(1), uSol);
  sLatticeK.setParameter<descriptors::OMEGA>(omega);

  sLatticeE.defineField<descriptors::SOURCE>(superGeometry.getMaterialIndicator(1), rho0);
  sLatticeE.defineField<descriptors::EPSILON>(superGeometry.getMaterialIndicator(1), eSol);
  sLatticeE.defineField<descriptors::OMEGA>(superGeometry.getMaterialIndicator(1), omega0);
  sLatticeE.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator(1), uSol);
  sLatticeE.setParameter<descriptors::OMEGA>(omega);

  sLatticeRANS.initialize();
  sLatticeK.initialize();
  sLatticeE.initialize();

  {
    auto& communicatorK = sLatticeK.getCommunicator(stage::PostCoupling());
    communicatorK.requestField<descriptors::VELOCITY>();
    communicatorK.requestField<descriptors::K>();
    communicatorK.requestOverlap(2);
    communicatorK.exchangeRequests();
  }

  {
    auto& communicatorE = sLatticeE.getCommunicator(stage::PostCoupling());
    communicatorE.requestField<descriptors::VELOCITY>();
    communicatorE.requestField<descriptors::EPSILON>();
    communicatorE.requestOverlap(2);
    communicatorE.exchangeRequests();
  }
}

#else
void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& sLattice,
                       UnitConverter<T,DESCRIPTOR> const& converter,
                       SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout,"setBoundaryValues");

  AnalyticalConst3D<T,T> rho(1.);
  Tgv3D<T,DESCRIPTOR> uSol(converter, 1);

  sLattice.defineRhoU(superGeometry, 1, rho, uSol);
  sLattice.iniEquilibrium(superGeometry, 1, rho, uSol);

  sLattice.initialize();
}
#endif

// Interpolate the data points to the output interval
void getDNSValues()
{
  std::string file_name;
  //Brachet, Marc E., et al. "Small-scale structure of the Taylor–Green vortex." Journal of Fluid Mechanics 130 (1983): 411-452; Figure 7
  if (util::abs(Re - 800.0) < std::numeric_limits<T>::epsilon() && maxPhysT <= 10.0 + std::numeric_limits<T>::epsilon()) {
    file_name= "Re800_Brachet.inp";
  }
  else if (util::abs(Re - 1600.0) < std::numeric_limits<T>::epsilon() && maxPhysT <= 10.0 + std::numeric_limits<T>::epsilon()) {
    file_name = "Re1600_Brachet.inp";
  }
  else if (util::abs(Re - 3000.0) < std::numeric_limits<T>::epsilon() && maxPhysT <= 10.0 + std::numeric_limits<T>::epsilon()) {
    file_name = "Re3000_Brachet.inp";
  }
  else {
    std::cout<<"Reynolds number not supported or maxPhysT>10: DNS plot will be disabled"<<std::endl;
    plotDNS = false;
    return;
  }
  std::ifstream data(file_name);
  std::string line;
  std::vector<std::vector<T>> parsedDat;
  while (std::getline(data, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    std::vector<T> parsedRow;
    while (std::getline(lineStream, cell, ' ')) {
      parsedRow.push_back(atof(cell.c_str()));

    }
    if (parsedDat.size() > 0 && parsedRow.size() > 1) {
      parsedRow.push_back((parsedRow[1] - parsedDat[parsedDat.size() - 1][1]) / (parsedRow[0] - parsedDat[parsedDat.size()-1][0]));
      parsedRow.push_back(parsedDat[parsedDat.size()-1][1] - parsedRow[2] * parsedDat[parsedDat.size()-1][0]);
    }
    parsedDat.push_back(parsedRow);
  }

  int steps = maxPhysT / vtkSave + 1.5;
  for (int i=0; i < steps; i++) {
    std::vector<T> inValues_temp;
    inValues_temp.push_back(i * vtkSave);
    if (inValues_temp[0] < parsedDat[0][0]) {
      inValues_temp.push_back(parsedDat[1][2] * inValues_temp[0] + parsedDat[1][3]);
    }
    else if (inValues_temp[0] > parsedDat[parsedDat.size()-1][0]) {
      inValues_temp.push_back(parsedDat[parsedDat.size()-1][2] * inValues_temp[0] +
                              parsedDat[parsedDat.size()-1][3]);
    }
    else {
      for (size_t j=0; j < parsedDat.size()-1; j++)  {
        if (inValues_temp[0] > parsedDat[j][0] && inValues_temp[0] < parsedDat[j+1][0]) {
          inValues_temp.push_back(parsedDat[j+1][2] * inValues_temp[0] + parsedDat[j+1][3]);
        }
      }
    }
    values_DNS.push_back(inValues_temp);
  }
}



void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
#ifdef K_EPSILON
                SuperLattice<T, KEDESCRIPTOR>& sLatticeK,
                SuperLattice<T, KEDESCRIPTOR>& sLatticeE,
#endif
                UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                SuperGeometry<T,3>& superGeometry, util::Timer<double>& timer)
{
  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter3D<T> vtmWriter("tgv3d");

  T omega = converter.getLatticeRelaxationFrequency();

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);
    superGeometry.rename(0,2);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
    if (plotDNS==true) {
      getDNSValues();
    }
  }

  static Gnuplot<T> gplot("Turbulence_Dissipation_Rate");

  if (iT % converter.getLatticeTime(vtkSave) == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);

#ifdef K_EPSILON
    SuperLatticeDensity3D<T, KEDESCRIPTOR> K(sLatticeK);
    SuperLatticeDensity3D<T, KEDESCRIPTOR> E(sLatticeE);

    K.getName() = "k";
    E.getName() = "e";

    vtmWriter.addFunctor( K );
    vtmWriter.addFunctor( E );
#endif
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write(iT);

    // write output of velocity as JPEG
    SuperEuklidNorm3D<T> normVel( velocity );
    BlockReduction3D2D<T> planeReduction( normVel, {0, 0, 1} );
    heatmap::write(planeReduction, iT);

    timer.update(iT);
    timer.printStep(2);
    sLattice.getStatistics().print(iT,converter.getPhysTime( iT ));

    int input[3];
    T output[1];

    ParametersOfOperatorD<T,DESCRIPTOR,BulkDynamics> bulkDynamicsParams{};
    bulkDynamicsParams.set<descriptors::OMEGA>(omega);
    if constexpr (BulkDynamics::parameters::contains<descriptors::LATTICE_TIME>()) {
      bulkDynamicsParams.set<descriptors::LATTICE_TIME>(iT);
    }

#if defined (finiteDiff)
    std::list<int> matNumber;
    matNumber.push_back(1);
    SuperLatticePhysDissipationFD3D<T, DESCRIPTOR> diss(superGeometry, sLattice, matNumber, converter);
#if !defined (DNS) && !defined(KBC) && !defined(RLB) && !defined(K_EPSILON)
    bulkDynamicsParams.set<collision::LES::SMAGORINSKY>(smagoConst);
    SuperLatticePhysEffectiveDissipationFD3D<T, DESCRIPTOR> effectiveDiss(
      superGeometry, sLattice, matNumber, converter,
      [&](Cell<T,DESCRIPTOR>& cell) -> double {
        return BulkDynamics::CollisionO().computeEffectiveOmega(cell, bulkDynamicsParams);
      });
#endif
#else
    SuperLatticePhysDissipation3D<T, DESCRIPTOR> diss(sLattice, converter);
#if !defined (DNS) && !defined(KBC) && !defined(RLB) && !defined(K_EPSILON)
    bulkDynamicsParams.set<collision::LES::SMAGORINSKY>(smagoConst);
    SuperLatticePhysEffectiveDissipation3D<T, DESCRIPTOR> effectiveDiss(
      sLattice, converter, smagoConst,
      [&](Cell<T,DESCRIPTOR>& cell) -> double {
        return BulkDynamics::CollisionO().computeEffectiveOmega(cell, bulkDynamicsParams);
      });
#endif
#endif
    SuperIntegral3D<T> integralDiss(diss, superGeometry, 1);
    integralDiss(output, input);
    T diss_mol = output[0];
    diss_mol /= volume;
    T diss_eff = diss_mol;

#if !defined (DNS) && !defined(KBC) && !defined(RLB) && !defined(K_EPSILON)
    SuperIntegral3D<T> integralEffectiveDiss(effectiveDiss, superGeometry, 1);
    integralEffectiveDiss(output, input);
    diss_eff = output[0];
    diss_eff /= volume;
#elif defined(K_EPSILON)
    const T dX = converter.getPhysDeltaX();
    const T dT = converter.getPhysDeltaT();
    T convE = dX*dX/dT/dT/dT;
    SuperLatticeDensity3D<T, KEDESCRIPTOR> effectiveDiss(sLatticeE);
    SuperIntegral3D<T> integralEffectiveDiss(effectiveDiss, superGeometry, 1);
    integralEffectiveDiss(output, input);
    diss_eff = output[0]*convE;
    diss_eff /= volume;
#endif

    T diss_eddy = diss_eff - diss_mol;
    if (plotDNS==true) {
      int step = converter.getPhysTime(iT) / vtkSave + 0.5;
      gplot.setData(converter.getPhysTime(iT), {diss_mol, diss_eddy, diss_eff, values_DNS[step][1]}, {"molecular dissipation rate", "eddy dissipation rate", "effective dissipation rate","Brachet et al."}, "bottom right");
    }
    else {
      gplot.setData(converter.getPhysTime(iT), {diss_mol, diss_eddy, diss_eff}, {"molecular dissipation rate", "eddy dissipation rate", "effective dissipation rate"}, "bottom right");
    }
    gplot.writePNG();
  }

  /// write pdf at last time step
  if (iT == converter.getLatticeTime(maxPhysT)-1) {
    gplot.writePDF();
  }
}

int main(int argc, char* argv[])
{

  /// === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout( std::cout,"main" );

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    int {int(std::nearbyint(N/(2*pi)))},        // resolution: number of voxels per charPhysL
    (T)   0.507639, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,        // charPhysLength: reference length of simulation geometry
    (T)   1,        // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./Re,    // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0       // physDensity: physical density in __kg / m^3__
  );


  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("tgv3d");

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif

  CuboidDecomposition3D<T> cuboidDecomposition(0, converter.getPhysDeltaX(), N, noOfCuboids);

  cuboidDecomposition.setPeriodicity({true, true, true});

  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // === 2nd Step: Prepare Geometry ===
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, 3);
  prepareGeometry(superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(sLattice, converter, superGeometry);

#ifdef K_EPSILON
  SuperLattice<T, KEDESCRIPTOR> sLatticeK(superGeometry);
  prepareLatticeKE(sLatticeK, converter, superGeometry);

  SuperLattice<T, KEDESCRIPTOR> sLatticeE(superGeometry);
  prepareLatticeKE(sLatticeE, converter, superGeometry);
#endif

#if defined(WALE)
  std::list<int> mat;
  mat.push_back(1);
  std::unique_ptr<SuperLatticeF3D<T, DESCRIPTOR>> functor(new SuperLatticeVelocityGradientFD3D<T, DESCRIPTOR>(superGeometry, sLattice, mat));
#endif
#if defined(K_EPSILON)
  SuperLatticeCoupling coupling(
    RANSKE<T>{},
    names::NavierStokes{}, sLattice,
    names::TurbKineticEnergy{}, sLatticeK,
    names::DissipationRate{}, sLatticeE);
    coupling.setParameter<RANSKE<T>::VISC>(converter.getLatticeViscosity());
    coupling.setParameter<RANSKE<T>::RNG>(T{1});
#endif

  /// === 4th Step: Main Loop with Timer ===
  util::Timer<double> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  // === 5th Step: Definition of Initial and Boundary Conditions ===
#if !defined(K_EPSILON)
  setBoundaryValues(sLattice, converter, superGeometry);
#else
  setBoundaryValues(sLattice, sLatticeK, sLatticeE, converter, superGeometry);
#endif

  for (std::size_t iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
    sLattice.setParameter<descriptors::LATTICE_TIME>(iT);
#if defined(WALE)
    sLattice.defineField<descriptors::VELO_GRAD>(superGeometry, 1, *functor);
#endif
#if defined(ShearSmagorinsky)
    sLattice.setParameter<descriptors::LATTICE_TIME>(iT);
#endif

    /// === 6th Step: Computation and Output of the Results ===
#if !defined(K_EPSILON)
    getResults(sLattice, converter, iT, superGeometry, timer);
#else
    getResults(sLattice, sLatticeK, sLatticeE, converter, iT, superGeometry, timer);
#endif

    /// === 7th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
#if defined(K_EPSILON)
    coupling.execute();
    sLatticeK.collideAndStream();
    sLatticeE.collideAndStream();
#endif
  }
  timer.stop();
  timer.printSummary();

  return 0;
}

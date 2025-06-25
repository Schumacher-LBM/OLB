/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Maximilian Gaedtke, Larissa Dietz
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

/* The solution for the melting problem (solid-liquid phase change)
coupled with natural convection is found using the lattice Boltzmann
method after Rongzong Huang and Huiying Wu (2015)[1]. The equilibrium
distribution function for the temperature is modified in order to deal
with the latent-heat source term. That way, iteration steps or solving
a group of linear equations is avoided, which results in enhanced efficiency.
The phase interface is located by the current total enthalpy, and
its movement is considered by the immersed moving boundary scheme after
Noble and Torczynski (1998)[2]. This method was validated by comparison
with experimental values (e.g. Gau and Viskanta (1986) [3]).

[1] Rongzong Huang, Huiying Wu, Phase interface effects in the total enthalpy-based lattice
Boltzmann model for solid–liquid phase change, Journal of Computational Physics 294 (2015) 346–362.

[2] D. Noble, J. Torczynski, A lattice-Boltzmann method for partially saturated
computational cells, Int. J. Modern Phys. C 9 (8) (1998) 1189–1202.

[3] C. Gau, R. Viskanta, Melting and Solidification of a Pure Metal on a
Vertikal Wall, Journal of Heat Transfer  108(1) (1986): 174–181.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;

using NSDESCRIPTOR = D2Q9<POROSITY,VELOCITY_SOLID,FORCE,OMEGA>;
using TDESCRIPTOR  = D2Q5<VELOCITY,TEMPERATURE>;

using TotalEnthalpyAdvectionDiffusionDynamics = TotalEnthalpyAdvectionDiffusionTRTdynamics<T,TDESCRIPTOR>;

// Parameters for the simulation setup
const T lx  = 88.9e-3;      // length of the channel
const T ly  = 63.5e-3;      // height of the channel
int     N = 128;            // resolution of the model
T       tau = 0.51;         // relaxation time
const T Ra = 2e6;           // Rayleigh number
const T Pr = 0.0216;        // Prandtl number
const T Ste = 0.039;        // Stephan number
const T maxPhysT = 1140.;   // simulation time

const T Tcold = 0.5;
const T Tmelt = (302.8 - 301.3)/(311.0 - 301.3) + 0.5;
const T Thot = 1.5;

const T lambda_s = 33.5; // W / m K
const T lambda_l = 32.0; // W / m K
const T R_lambda =  lambda_s/lambda_l;

const T cp_s = 1.0; // J / kg K
const T cp_l = 1.0; // J / kg K
const T R_cp = cp_s/cp_l;

//for this case, the harmonic mean (cp_ref) is applicable
const T cp_ref = 2.0 * cp_s * cp_l / (cp_s + cp_l); // J / kg K

const T R_alpha = lambda_s / lambda_l * cp_l / cp_s;
const T density = 1.; // kg / m^3 //?
const T L = cp_l * (Thot - Tmelt) / Ste; // J / kg

T lattice_Hcold, lattice_Hhot;
T physDeltaX, physDeltaT;

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry<T,2>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,4);

  std::vector<T> extend(2,T());
  extend[0] = lx;
  extend[1] = ly;
  std::vector<T> origin(2,T());
  origin[0] = converter.getPhysLength(1);
  origin[1] = 0.5*converter.getPhysLength(1);
  IndicatorCuboid2D<T> cuboid2(extend, origin);

  superGeometry.rename(4,1,cuboid2);

  std::vector<T> extendwallleft(2,T(0));
  extendwallleft[0] = converter.getPhysLength(1);
  extendwallleft[1] = ly;
  std::vector<T> originwallleft(2,T(0));
  originwallleft[0] = 0.0;
  originwallleft[1] = 0.0;
  IndicatorCuboid2D<T> wallleft(extendwallleft, originwallleft);

  std::vector<T> extendwallright(2,T(0));
  extendwallright[0] = converter.getPhysLength(1);
  extendwallright[1] = ly;
  std::vector<T> originwallright(2,T(0));
  originwallright[0] = lx+converter.getPhysLength(1);
  originwallright[1] = 0.0;
  IndicatorCuboid2D<T> wallright(extendwallright, originwallright);

  superGeometry.rename(4,2,1,wallleft);
  superGeometry.rename(4,3,1,wallright);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;

}

template <typename SuperLatticeCoupling>
void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                     SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice<T, TDESCRIPTOR>& ADlattice,
                     SuperLatticeCoupling& coupling,
                     SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  T omega  = converter.getLatticeRelaxationFrequency();
  T Tomega = converter.getLatticeThermalRelaxationFrequency();

  NSlattice.defineDynamics<ForcedPSMBGKdynamics>(superGeometry.getMaterialIndicator({1, 2, 3, 4}));
  ADlattice.defineDynamics<TotalEnthalpyAdvectionDiffusionDynamics>(superGeometry.getMaterialIndicator({1, 2, 3}));
  boundary::set<boundary::BounceBack>(ADlattice, superGeometry, 4);

  /// sets boundary
  boundary::set<boundary::RegularizedTemperature>(ADlattice, superGeometry.getMaterialIndicator({2,3}));
  boundary::set<boundary::InterpolatedVelocity>(NSlattice, superGeometry.getMaterialIndicator({2,3,4}));

  /// Compute pre factor
  std::vector<T> dir{0.0, 1.0};
  std::vector<T> forcePrefactor{0, 0};
  T boussinesqForcePrefactor = Ra / util::pow(T(N),3) * Pr * util::pow(cp_ref / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5), 2);
  clout << "boussinesq " << Ra / util::pow(T(N), 3) * Pr * lambda_l * lambda_l << std::endl;

  // we normalize the direction of force vector
  T normDir = T();
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    normDir += dir[iD]*dir[iD];
  }
  normDir = util::sqrt(normDir);
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    dir[iD] /= normDir;
  }

  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    forcePrefactor[iD] = boussinesqForcePrefactor * dir[iD];
  }

  NSlattice.setParameter<descriptors::OMEGA>(omega);
  ADlattice.setParameter<descriptors::OMEGA>(Tomega);

  ADlattice.setParameter<collision::TRT::MAGIC>(T(0.25));

  ADlattice.setParameter<TotalEnthalpy::T_S>(Tmelt);
  ADlattice.setParameter<TotalEnthalpy::T_L>(Tmelt);
  ADlattice.setParameter<TotalEnthalpy::CP_S>(cp_s);
  ADlattice.setParameter<TotalEnthalpy::CP_L>(cp_l);
  ADlattice.setParameter<TotalEnthalpy::LAMBDA_S>(cp_ref / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5) * R_lambda);
  ADlattice.setParameter<TotalEnthalpy::LAMBDA_L>(cp_ref / descriptors::invCs2<T,TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5));
  ADlattice.setParameter<TotalEnthalpy::L>(L);

  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::T_S>(Tmelt);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::T_L>(Tmelt);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::CP_S>(cp_s);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::CP_L>(cp_l);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::L>(L);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::FORCE_PREFACTOR>(forcePrefactor);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::T_COLD>(Tcold);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::DELTA_T>(T(1.));

  /// define initial conditions
  AnalyticalConst2D<T,T> rho(1.);
  AnalyticalConst2D<T,T> u0(0.0, 0.0);
  AnalyticalConst2D<T,T> T_cold(lattice_Hcold);
  AnalyticalConst2D<T,T> T_hot(lattice_Hhot);
  AnalyticalConst2D<T,T> omegaField(omega);
  NSlattice.defineField<descriptors::OMEGA>(superGeometry.getMaterialIndicator({1, 2, 3, 4}), omegaField);

  /// for each material set Rho, U and the Equilibrium
  NSlattice.defineRhoU(superGeometry.getMaterialIndicator({1, 2, 3, 4}), rho, u0);
  NSlattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 2, 3, 4}), rho, u0);

  ADlattice.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({1, 2, 3}), u0);
  ADlattice.defineRho(superGeometry.getMaterialIndicator({1, 3}), T_cold);
  ADlattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 3}), T_cold, u0);
  ADlattice.defineRho(superGeometry, 2, T_hot);
  ADlattice.iniEquilibrium(superGeometry, 2, T_hot, u0);

  /// Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                        SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                        SuperLattice<T, TDESCRIPTOR>& ADlattice,
                        int iT, SuperGeometry<T,2>& superGeometry)
{

// nothing to do here

}

void getResults( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                 SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice<T, TDESCRIPTOR>& ADlattice, int iT,
                 SuperGeometry<T,2>& superGeometry,
                 util::Timer<T>& timer,
                 bool converged)
{

  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter2D<T> vtkWriter("galliumMelting2d");
  SuperLatticeField2D<T, TDESCRIPTOR, VELOCITY> velocity(ADlattice);
  SuperLatticePhysPressure2D<T, NSDESCRIPTOR> pressure(NSlattice, converter);

  SuperLatticeDensity2D<T, TDESCRIPTOR> enthalpy(ADlattice);
  enthalpy.getName() = "enthalpy";
  SuperLatticeField2D<T, NSDESCRIPTOR, POROSITY> liquid_frac(NSlattice);
  liquid_frac.getName() = "liquid fraction";
  SuperLatticeField2D<T, TDESCRIPTOR, TEMPERATURE> temperature(ADlattice);
  temperature.getName() = "temperature";
  SuperLatticeField2D<T, NSDESCRIPTOR, FORCE> force(NSlattice);
  force.getName() = "force";
  vtkWriter.addFunctor( pressure );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( enthalpy );
  vtkWriter.addFunctor( liquid_frac );
  vtkWriter.addFunctor( temperature );
  vtkWriter.addFunctor( force );

  const int vtkIter = converter.getLatticeTime(0.5);

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files
  if (iT % vtkIter == 0 || converged) {

    timer.update(iT);
    timer.printStep();

    /// NSlattice statistics console output
    NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));

    /// ADlattice statistics console output
    ADlattice.getStatistics().print(iT,converter.getPhysTime(iT));

    /*if ( NSlattice.getStatistics().getAverageRho() != NSlattice.getStatistics().getAverageRho() or ADlattice.getStatistics().getAverageRho() != ADlattice.getStatistics().getAverageRho() ) {
      clout << "simulation diverged! stopping now." << std::endl;
      exit(1);
    }*/

    vtkWriter.write(iT);
  }
}

int main(int argc, char *argv[])
{

  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  T char_lattice_u = 0.2;

  if (argc >= 2) {
    N = atoi(argv[1]);
  }
  if (argc >= 3) {
    tau = atof(argv[2]);
  }
  if (argc >= 4) {
    char_lattice_u = atof(argv[3]);
  }

  const T char_u = util::sqrt( 9.81 * 1.2e-4 * (311. - 302.8) * 6093. );
  const T conversion_u = char_u / char_lattice_u;

  physDeltaX = lx / N;
  physDeltaT = physDeltaX / conversion_u;
  physDeltaT = 6093. / 1.81e-3 / descriptors::invCs2<T,NSDESCRIPTOR>() * (tau - 0.5) * physDeltaX * physDeltaX;

  lattice_Hcold = cp_s * Tcold;
  lattice_Hhot = cp_l * Thot;

  clout << "H_cold " << lattice_Hcold << " H_hot " << lattice_Hhot << std::endl;

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const converter(
    (T) physDeltaX, // physDeltaX
    (T) physDeltaT, // physDeltaT
    (T) lx, // charPhysLength
    (T) char_u, // charPhysVelocity
    (T) 1.81e-3 / 6093., // physViscosity
    (T) 6093., // physDensity
    (T) 32., // physThermalConductivity
    (T) 381., // physSpecificHeatCapacity
    (T) 1.2e-4, // physThermalExpansionCoefficient
    (T) Tcold, // charPhysLowTemperature
    (T) Thot // charPhysHighTemperature
  );
  converter.print();
  clout << "lattice cp " << converter.getLatticeSpecificHeatCapacity(cp_l) << std::endl;

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  extend[0] = lx + 2*converter.getPhysLength(1);
  extend[1] = ly + converter.getPhysLength(1);
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of an empty cuboidDecomposition
  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  /// Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer);

  prepareGeometry(superGeometry, converter);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice<T, TDESCRIPTOR> ADlattice(superGeometry);
  SuperLattice<T, NSDESCRIPTOR> NSlattice(superGeometry);

  SuperLatticeCoupling coupling(
      TotalEnthalpyPhaseChangeCoupling{},
      names::NavierStokes{}, NSlattice,
      names::Temperature{}, ADlattice);
  coupling.restrictTo(superGeometry.getMaterialIndicator({1}));

  // prepareLattice and setBoundaryConditions
  prepareLattice(converter, NSlattice, ADlattice, coupling, superGeometry);

  /// === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT)+1; ++iT) {

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, NSlattice, ADlattice, iT, superGeometry);

    /// === 6th Step: Collide and Stream Execution ===
    coupling.execute();
    NSlattice.collideAndStream();
    ADlattice.collideAndStream();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, false);
  }

  timer.stop();
  timer.printSummary();

}

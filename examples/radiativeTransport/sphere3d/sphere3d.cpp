/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014-2018 Albert Mink
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

/* sphere3d.cpp:
 * A 3D implementation of radiation being emitted from a sphere and
 * spreading in a bigger sphere through highly scattering media. The macroscopic
 * method is based on the diffusion approximation.
 * The theoretical background and validation of said method are detailed in
 * [A. Mink, G. Thäter, H. Nirschl, and M. J. Krause. “A 3D Lattice Boltzmann method
 * for light simulation in participating media”. In: Journal of Computational Science
 * 17.Part 2 (2016). Discrete Simulation of Fluid Dynamics 2015, pp. 431–437. DOI:
 * 10.1016/j.jocs.2016.03.014.]
 */

#include "olb3D.h"
#include "olb3D.hh"     // include full template code

#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;

typedef double T;

#define DESCRIPTOR descriptors::D3Q7<descriptors::VELOCITY> // descriptors::D3Q27DescriptorLebedev

// PostProcessor for updating the source term in each cell
struct ComputeSourceTerm {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct C_T : public descriptors::FIELD_BASE<1> {};
  struct MU_EFF : public descriptors::FIELD_BASE<1> {};

  using parameters = meta::list<C_T, MU_EFF>;

  int getPriority() const { return 0; }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V  = typename CELL::value_t;
    V phi    = cell.computeRho();
    V mu_eff = parameters.template get<MU_EFF>();
    V source = -mu_eff * mu_eff * phi * parameters.template get<C_T>();
    cell.template setField<descriptors::SOURCE>(source);
  }
};

// Compute the analytical solution and error
void error(SuperGeometry<T, 3>& superGeometry, SuperLattice<T, DESCRIPTOR>& lattice,
           UnitConverter<T, DESCRIPTOR> const& converter, T absorption, T scattering)
{
  OstreamManager clout(std::cout, "error");

  int input[1];
  T   normAnaSol[1], absErr[1], relErr[1];
  // analytical solution
  PLSsolution3D<T, T, DESCRIPTOR>               dSol(absorption, scattering);
  SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> dSolLattice(dSol, lattice);
  // approximated solution
  SuperLatticeDensity3D<T, DESCRIPTOR> d(lattice);

  SuperL2Norm3D<T> dL2Norm(dSolLattice - d, superGeometry, 1);
  SuperL2Norm3D<T> dSolL2Norm(dSolLattice, superGeometry, 1);
  dL2Norm(absErr, input);
  dSolL2Norm(normAnaSol, input);
  relErr[0] = absErr[0] / normAnaSol[0];
  clout << "density-L2-error(util::abs)=" << absErr[0] << "; density-L2-error(rel)=" << relErr[0] << std::endl;

  SuperLinfNorm3D<T> dLinfNorm(dSolLattice - d, superGeometry, 1);
  dLinfNorm(absErr, input);
  clout << "density-Linf-error(util::abs)=" << absErr[0] << std::endl;
}

// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter, SuperGeometry<T, 3>& superGeometry,
                     IndicatorSphere3D<T>& sphereAll, IndicatorSphere3D<T>& sphereEm)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  // Material numbers from 0 to 3 defined by indicators
  // the domain, a sphere of radius 1
  superGeometry.rename(0, 2, sphereAll);
  superGeometry.rename(2, 1, {1, 1, 1});
  // the emitter sphere of radius 0.1
  superGeometry.rename(1, 3, sphereEm);

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

// Set up the geometry of the simulation
void prepareLattice(UnitConverter<T, DESCRIPTOR> const& converter, SuperLattice<T, DESCRIPTOR>& sLattice,
                    SuperGeometry<T, 3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // the simulation parameters
  const T omega = converter.getLatticeRelaxationFrequency();
  // T latticeSink = 3.*converter.getLatticeAbsorption()*(converter.getLatticeAbsorption()+converter.getLatticeScattering()) / 8.;

  // define dynamics
  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics<T, DESCRIPTOR>>(superGeometry, 0);

  // Material=1 -->sourced advection diffusion dynamics OR poisson dynamics with sink term
  // when using poisson dynamics, setting the sink parameter is necessary
  // when using sourced advection diffusion dynamics, the sink term is handeled by the post processor
  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  // sLattice.defineDynamics<PoissonDynamics<T,DESCRIPTOR>>(bulkIndicator);
  sLattice.defineDynamics<SourcedAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>>(bulkIndicator);

  // Material=2,3 -->first order equilibrium
  sLattice.defineDynamics<EquilibriumBoundaryFirstOrder>(superGeometry, 2);
  sLattice.defineDynamics<EquilibriumBoundaryFirstOrder>(superGeometry, 3);

  // set the parameters for the dynamics
  sLattice.setParameter<descriptors::OMEGA>(omega);
  // sLattice.setParameter<collision::Poisson::SINK>(latticeSink);

  AnalyticalConst3D<T, T> velocity(0.0, 0.0, 0.0);
  sLattice.defineField<descriptors::VELOCITY>(bulkIndicator, velocity);

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

// Write data to termimal and file system
void getResults(SuperLattice<T, DESCRIPTOR>& sLattice, int iT, UnitConverter<T, DESCRIPTOR> const& converter,
                SuperGeometry<T, 3>& superGeometry, util::Timer<T>& timer, int saveIter, T absorption, T scattering)
{
  OstreamManager                       clout(std::cout, "getResults");
  SuperVTMwriter3D<T>                  vtmWriter("sphere3d");
  SuperLatticeDensity3D<T, DESCRIPTOR> density(sLattice);
  vtmWriter.addFunctor(density);

  if (iT == 0) {
    // Writes geometry, cuboid no. and rank no. to file system
    SuperLatticeCuboid3D cuboid(sLattice);
    SuperLatticeRank3D   rank(sLattice);

    vtmWriter.write(cuboid);
    vtmWriter.write(rank);

    vtmWriter.createMasterFile();
  }

  if (iT % saveIter == 0) {
    std::cout << "Write added functor at time " << iT << "." << std::endl;
    vtmWriter.write(iT);
  }

  // Writes output on the console
  if (iT % saveIter == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  if (iT % 500 == 0) {
    SuperLatticeDensity3D<T, DESCRIPTOR> density(sLattice);
    AnalyticalFfromSuperF3D<T>           analytDensity(density, true, 1);
    PLSsolution3D<T, T, DESCRIPTOR>      dSol(absorption, scattering);
    for (int nY = -25; nY < 25; ++nY) {
      double position[3] = {0, nY / 25., 0};
      // double postition[3] = {0, 0, 0};
      double densityEval[1] = {0};
      analytDensity(densityEval, position);
      double densitySol[1] = {0};
      dSol(densitySol, position);
      // printf("%1.3f\t%1.3f\t%1.3f\n", position[1], densityEval[0], densitySol[0] );
      // fout << postition[1] << ", " << densityEval[0] << ", " << densitySol[0]  << std::endl;
    }
    // fout.close();
  }
}

void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& lattice, SuperGeometry<T, 3>& superGeometry, T absorption,
                       T scattering)
{
  // since adv-diffusion model is used, the velocity it set to 0
  AnalyticalConst3D<T, T>         u0(0.0, 0.0, 0.0);            // 3D -> 3D
  AnalyticalConst3D<T, T>         rho0(0.0);                    // 3D -> 3D
  PLSsolution3D<T, T, DESCRIPTOR> dSol(absorption, scattering); // analytical Solution

  // initialize media with density from analytical solution
  // at iT=0 the error is given by the maschinen genauigkeit
  lattice.iniEquilibrium(superGeometry, 1, rho0, u0);
  lattice.iniEquilibrium(superGeometry, 2, dSol, u0);
  lattice.iniEquilibrium(superGeometry, 3, dSol, u0);
  lattice.defineRho(superGeometry, 2, dSol);
  lattice.defineRho(superGeometry, 3, dSol);
  lattice.initialize();
}

int main(int argc, char* argv[])
{
  // ===== 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  std::string fName("sphere3d.xml");
  XMLreader   config(fName);

  int N, maxPhysT, saveIter;
  T   writeSec, ABSORPTION, SCATTERING;
  config["Application"]["DiscParam"]["resolution"].read(N);
  config["Application"]["PhysParam"]["maxTime"].read(maxPhysT);
  config["Application"]["absorption"].read(ABSORPTION);
  config["Application"]["scattering"].read(SCATTERING);
  config["Output"]["saveSec"].read(writeSec);

  T latticeRelaxationTime = 1.0;
  T mu_eff                = util::sqrt(3 * ABSORPTION * (ABSORPTION + SCATTERING));
  // the radiative unit converter specifically deals with values associated with radiation, necessary when using poisson dynamics
  // const RadiativeUnitConverter<T,DESCRIPTOR> converter( N, 1, ABSORPTION, SCATTERING);
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
      (int)N,                   // resolution: number of voxels per charPhysL
      (T)latticeRelaxationTime, // latticeRelaxationTime: relaxation time, has to be greater than 0.5!
      (T)1.0,                   // charPhysLength: reference length of simulation geometry
      (T)0.0, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__  || none for radiation
      (T)1.0, // physViscosity: physical kinematic viscosity in __m^2 / s__  || radiative diffusion coefficient D
      (T)1.0  // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  saveIter = converter.getLatticeTime(writeSec / 32.);

  // ===== 2nd Step: Prepare Geometry =====
  Vector<T, 3>         origin;
  IndicatorSphere3D<T> sphereEm(origin, 0.1);
  IndicatorSphere3D<T> sphereAll(origin, 1);

  // Instantiation of a cuboid
  CuboidDecomposition<T, 3>* cuboid =
      new CuboidDecomposition<T, 3>(sphereAll, converter.getConversionFactorLength(), 2 * singleton::mpi().getSize());
  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T>* loadBalancer = new HeuristicLoadBalancer<T>(*cuboid);
  // Instatiation of a superGeometry
  SuperGeometry<T, 3> superGeometry(*cuboid, *loadBalancer, 2);

  prepareGeometry(converter, superGeometry, sphereAll, sphereEm);

  // ===== 3rd Step: Prepare Lattice =====
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(converter, sLattice, superGeometry);

  // Add PostProcessor and set parameters
  sLattice.addPostProcessor<stage::PreCollide>(meta::id<ComputeSourceTerm> {});
  sLattice.setParameter<ComputeSourceTerm::C_T>(converter.getConversionFactorTime());
  sLattice.setParameter<ComputeSourceTerm::MU_EFF>(mu_eff);

  // === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());
  timer.start();

  setBoundaryValues(sLattice, superGeometry, ABSORPTION, SCATTERING);

  util::ValueTracer<T> converge(N, 1e-6);
  int                  iT;
  for (iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    converge.takeValue(sLattice.getStatistics().getAverageRho(), true);
    if (converge.hasConverged()) {
      clout << "Simulation converged. -- " << iT << std::endl;
      break;
    }
    if (iT % saveIter == 0) {
      sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
      timer.print(iT);
      error(superGeometry, sLattice, converter, ABSORPTION, SCATTERING);
    }

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, iT, converter, superGeometry, timer, saveIter, ABSORPTION, SCATTERING);
  }
  // Write and save the results in a text file
  std::string   name = std::string("line") + std::to_string(iT) + std::string(".txt");
  std::ofstream fout(name, std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << name << std::endl;
  }
  SuperLatticeDensity3D<T, DESCRIPTOR> density(sLattice);
  AnalyticalFfromSuperF3D<T>           analytDensity(density, true, 1);
  PLSsolution3D<T, T, DESCRIPTOR>      dSol(ABSORPTION, SCATTERING);
  int                                  n = converter.getResolution() / 2;
  for (int nY = 0; nY < n; ++nY) {
    double position[3]    = {0, double(nY) / n, 0};
    double densityEval[1] = {0};
    analytDensity(densityEval, position);
    double densitySol[1] = {0};
    dSol(densitySol, position);
    if (singleton::mpi().getRank() == 0) {
      fout << position[1] << ", " << densityEval[0] << ", " << densitySol[0] << ", "
           << util::abs(densityEval[0] - densitySol[0]) / densitySol[0] << std::endl;
      printf("%1.3f\t%1.3f\t%1.3f\t%1.4f\n", position[1], densityEval[0], densitySol[0],
             util::abs(densityEval[0] - densitySol[0]) / densitySol[0]);
    }
  }
  fout.close();

  timer.stop();
  timer.printSummary();

  return 0;
}

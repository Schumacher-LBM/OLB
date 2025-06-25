/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

using namespace olb;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = descriptors::DualForcedD3Q19Descriptor;
using DYNAMICS = ForcedBGKdynamics<T,DESCRIPTOR>;

const int resolution = 11;
const T latticeVelocity = 0.07;
const T physLength = 1;
const T physVelocity = 1;
const T physViscosity = 0.1;
const T physDensity = 1;
const T physMaxT = 6.0;
const T physStartT = 4.0;

void prepareGeometry(SuperGeometry<T,3>& superGeometry) {
  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, {1,1,1});
  superGeometry.print();
}

void prepareLattice(UnitConverter<T,DESCRIPTOR>& converter,
                    SuperGeometry<T,3>& superGeometry,
                    SuperLattice<T,DESCRIPTOR>& sLattice) {
  // Define physics
  sLattice.template defineDynamics<DYNAMICS>(superGeometry.getMaterialIndicator({1}));
  boundary::set<boundary::LocalVelocity>(sLattice, superGeometry.getMaterialIndicator({2}));
  sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  // Set force field
  ForceTestFlow3D<T,T,DESCRIPTOR> forceSol(converter);
  const T latticeScaling(converter.getConversionFactorMass() / converter.getConversionFactorForce());
  AnalyticalScaled3D<T,T> force(forceSol, latticeScaling);  // conversion to lattice units
  sLattice.template defineField<descriptors::FORCE>(superGeometry.getMaterialIndicator({1}), force);

  // Set initial values
  AnalyticalConst3D<T,T> rhoF(1);
  Vector<T,3> u{0,0,0};
  AnalyticalConst3D<T,T> uF(u);
  sLattice.defineRhoU(superGeometry.getMaterialIndicator({1,2}), rhoF, uF);
  sLattice.iniEquilibrium(superGeometry.getMaterialIndicator({1,2}), rhoF, uF);
  sLattice.stripeOffDensityOffset(sLattice.getStatistics().getAverageRho() - T{1});

  sLattice.initialize();
}

void setBoundaryValues(UnitConverter<T,DESCRIPTOR>& converter,
                       SuperLattice<T,DESCRIPTOR>& sLattice,
                       std::size_t iT,
                       SuperGeometry<T,3>& superGeometry) {
  const std::size_t itStart = converter.getLatticeTime(physStartT);
  if (iT <= itStart) {
    // compute scaling factor
    PolynomialStartScale<T,T> startScale(itStart, T(1));
    T iTvec[1] = {T(iT)};
    T frac[1] = {};
    startScale(frac, iTvec);

    // analytical solution
    VelocityTestFlow3D<T,T,DESCRIPTOR> velocity(converter);
    AnalyticalScaled3D<T,T> uBoundaryStart(velocity, frac[0] / converter.getConversionFactorVelocity());
    sLattice.defineU(superGeometry, 2, uBoundaryStart);
  }
}

void getResults(UnitConverter<T,DESCRIPTOR>& converter,
                SuperLattice<T,DESCRIPTOR>& sLattice,
                LatticeResults<T,DESCRIPTOR>& results,
                std::size_t iT, util::Timer<T>& timer) {
  const std::size_t iTvtk = converter.getLatticeTime(physMaxT/5);
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/5);

  // writes the vtk files
  if (iT%iTvtk == 0) {
    results.write(converter, sLattice, iT);
  }

  // Get statistics
  if (iT%iTlog == 0) {
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

auto build(std::string name) {
  // === 1st Step: Initialization ===
  OstreamManager clout(std::cout, name);
  clout << "Building simulation setup..." << std::endl;
  auto lData = std::make_unique<LatticeData<T,DESCRIPTOR>>(name);

  // Provide the unit converter the characteristic entities
  auto& converter = lData->create<UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR>>(
    resolution,        // resolution
    latticeVelocity,   // charLatticeVelocity
    physLength,        // charPhysLength: reference length of simulation geometry in [m]
    physVelocity,      // charPhysVelocity: highest expected velocity during simulation in [m/s]
    physViscosity,     // physViscosity: physical kinematic viscosity in [m^2/s]
    physDensity        // physDensity: physical density [kg/m^3]
  );
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  Vector extend{physLength, physLength, physLength};
  Vector origin{0, 0, 0};
  IndicatorCuboid3D<T> cuboid(extend, origin);
  auto& cuboidDecomposition = lData->create<CuboidDecomposition3D<T>>(cuboid, converter.getPhysDeltaX(),singleton::mpi().getSize());

  auto& loadBalancer = lData->create<HeuristicLoadBalancer<T>>(cuboidDecomposition);

  auto& superGeometry = lData->create<SuperGeometry<T,3>>(cuboidDecomposition, loadBalancer);
  prepareGeometry(superGeometry);

  lData->create<SuperLattice<T,DESCRIPTOR>>(superGeometry);
  clout << "Done." << std::endl;
  return lData;
}

void simulate(std::unique_ptr<LatticeData<T,DESCRIPTOR>>& lData) {
  auto& converter = lData->getUnitConverter();
  auto& superGeometry = lData->getSuperGeometry();
  auto& sLattice = lData->getSuperLattice();
  auto& results = lData->getLatticeResults();

  OstreamManager clout(std::cout, lData->getName());
  clout << "Starting simulation..." << std::endl;

  // === 4th Step: Main Loop with Timer ===
  const std::size_t iTmax = converter.getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, sLattice, iT, superGeometry);
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults(converter, sLattice, results, iT, timer);
  }

  timer.stop();
  timer.printSummary();
  clout << "Done." << std::endl;
}

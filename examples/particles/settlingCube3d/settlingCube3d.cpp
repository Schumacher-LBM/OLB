/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2021 Nicolas Hafen, Mathias J. Krause
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

/* settlingCube3d.cpp:
 * The case examines the settling of a cubical silica particle
 * under the influence of gravity.
 * The object is surrounded by water in a rectangular domain
 * limited by no-slip boundary conditions.
 * For the calculation of forces an DNS approach is chosen
 * which also leads to a back-coupling of the particle on the fluid,
 * inducing a flow.
 *
 * The simulation is based on the homogenised lattice Boltzmann approach
 * (HLBM) introduced in "Particle flow simulations with homogenised
 * lattice Boltzmann methods" by Krause et al.
 * and extended in "Towards the simulation of arbitrarily shaped 3D particles
 * using a homogenised lattice Boltzmann method" by Trunk et al.
 * for the simulation of 3D particles.
 *
 * This example demonstrates the usage of HLBM in the OpenLB framework.
 * To improve parallel performance, the particle decomposition scheme
 * described in 10.48550/arXiv.2312.14172 is used.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace olb::particles;
using namespace olb::particles::dynamics;

using T = FLOATING_POINT_TYPE;

// Define lattice type
typedef PorousParticleD3Q19Descriptor DESCRIPTOR;

// Define particleType
#ifdef PARALLEL_MODE_MPI
// Particle decomposition improves parallel performance
typedef ResolvedDecomposedParticle3D PARTICLETYPE;
#else
typedef ResolvedParticle3D PARTICLETYPE;
#endif

#define WriteVTK

// Discretization Settings
int res = 30;
T const charLatticeVelocity = 0.01;

// Time Settings
T const maxPhysT = 0.5;       // max. simulation time in s
T const iTwrite = 0.02;       // write out intervall in s

// Domain Settings
T const lengthX = 0.01;
T const lengthY = 0.01;
T const lengthZ = 0.05;

// Fluid Settings
T const physDensity = 1000;
T const physViscosity = 1E-5;

//Particle Settings
T centerX = lengthX*.5;
T centerY = lengthY*.5;
T centerZ = lengthZ*.9;
T const cubeDensity = 2500;
T const cubeEdgeLength = 0.0025;
Vector<T,3> cubeCenter = {centerX,centerY,centerZ};
Vector<T,3> cubeOrientation = {0.,15.,0.};
Vector<T,3> cubeVelocity = {0.,0.,0.};
Vector<T,3> externalAcceleration = {.0, .0, -T(9.81) * (T(1) - physDensity / cubeDensity)};

// Characteristic Quantities
T const charPhysLength = lengthX;
T const charPhysVelocity = 0.15;    // Assumed maximal velocity


// Prepare geometry
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, {1, 1, 1});

  superGeometry.clean();
  superGeometry.innerClean();

  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}


// Set up the geometry of the simulation
void prepareLattice(
  SuperLattice<T, DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR> const& converter,
  SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  clout << "setting Velocity Boundaries ..." << std::endl;

  /// Material=0 -->do nothing
  sLattice.defineDynamics<PorousParticleBGKdynamics>(superGeometry, 1);
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  {
    auto& communicator = sLattice.getCommunicator(stage::PostPostProcess());
    communicator.requestFields<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR>();
    communicator.requestOverlap(sLattice.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}


//Set Boundary Values
void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& sLattice,
                       UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                       SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  if (iT == 0) {
    AnalyticalConst3D<T, T> zero(0.);
    AnalyticalConst3D<T, T> one(1.);
    sLattice.defineField<descriptors::POROSITY>(superGeometry.getMaterialIndicator({0,1,2}), one);
    // Set initial condition
    AnalyticalConst3D<T, T> ux(0.);
    AnalyticalConst3D<T, T> uy(0.);
    AnalyticalConst3D<T, T> uz(0.);
    AnalyticalConst3D<T, T> rho(1.);
    AnalyticalComposed3D<T, T> u(ux, uy, uz);

    //Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU(superGeometry, 1, rho, u);
    sLattice.iniEquilibrium(superGeometry, 1, rho, u);

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}


/// Computes the pressure drop between the voxels before and after the cylinder
void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T,3>& superGeometry, Timer<T>& timer,
                XParticleSystem<T, PARTICLETYPE>& xParticleSystem )
{
  OstreamManager clout(std::cout, "getResults");

#ifdef WriteVTK
  SuperVTMwriter3D<T> vtkWriter("sedimentation");
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticePhysExternalPorosity3D<T, DESCRIPTOR> externalPor(sLattice, converter);
  vtkWriter.addFunctor(velocity);
  vtkWriter.addFunctor(pressure);
  vtkWriter.addFunctor(externalPor);

  if (iT == 0) {
    /// Writes the converter log file
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    vtkWriter.write(iT);
  }
#endif

  /// Writes output on the console
  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
#ifdef PARALLEL_MODE_MPI
    communication::forParticlesInSuperParticleSystem<
        T, PARTICLETYPE, conditions::valid_particle_centres>(
        xParticleSystem,
        [&](Particle<T, PARTICLETYPE>&       particle,
            ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
          io::printResolvedParticleInfo(particle);
        });
#else
    for (std::size_t iP=0; iP<xParticleSystem.size(); ++iP) {
      auto particle = xParticleSystem.get(iP);
      io::printResolvedParticleInfo(particle);
    }
#endif
  }
}

int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR> converter(
    (int)   res,                  //resolution
    ( T )   charLatticeVelocity,  //charLatticeVelocity
    ( T )   charPhysLength,       //charPhysLength
    ( T )   charPhysVelocity,     //charPhysVelocity
    ( T )   physViscosity,        //physViscosity
    ( T )   physDensity           //physDensity
  );
  converter.print();

  /// === 2rd Step: Prepare Geometry ===
  /// Instantiation of a cuboidDecomposition with weights
  Vector<T,3> origin( 0. );
  Vector<T,3> extend( lengthX, lengthY, lengthZ );
  IndicatorCuboid3D<T> cuboid(extend, origin);

#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition3D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());
#else
  CuboidDecomposition3D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), 7);
#endif
  cuboidDecomposition.print();

  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, 2);
  prepareGeometry(converter, superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);

  // Prepare lattice
  prepareLattice(sLattice, converter, superGeometry);

  // Create ParticleSystem
#ifdef PARALLEL_MODE_MPI
  SuperParticleSystem<T,PARTICLETYPE> particleSystem(superGeometry);
#else
  ParticleSystem<T,PARTICLETYPE> particleSystem;
#endif

  //Create particle manager handling coupling, gravity and particle dynamics
  ParticleManager<T,DESCRIPTOR,PARTICLETYPE> particleManager(
    particleSystem, superGeometry, sLattice, converter, externalAcceleration);

  // Create and assign resolved particle dynamics
  particleSystem.defineDynamics<
    VerletParticleDynamics<T,PARTICLETYPE>>();

  // Calculate particle quantities
  T epsilon = 0.5*converter.getPhysDeltaX();
  Vector<T,3> cubeExtend( cubeEdgeLength );

  // Create Particle 1
  creators::addResolvedCuboid3D( particleSystem, cubeCenter,
                                 cubeExtend, epsilon, cubeDensity, cubeOrientation );

  // Create Particle 2
  cubeCenter = {centerX,lengthY*T(0.51),lengthZ*T(.7)};
  cubeOrientation = {0.,0.,15.};
  creators::addResolvedCuboid3D( particleSystem, cubeCenter,
                                 cubeExtend, epsilon, cubeDensity, cubeOrientation );

  // Check ParticleSystem
  particleSystem.checkForErrors();

  /// === 4th Step: Main Loop with Timer ===
  Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());
  timer.start();

  /// === 5th Step: Definition of Initial and Boundary Conditions ===
  setBoundaryValues(sLattice, converter, 0, superGeometry);

  clout << "MaxIT: " << converter.getLatticeTime(maxPhysT) << std::endl;

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT)+10; ++iT) {

    // Execute particle manager
    particleManager.execute<
      couple_lattice_to_particles<T,DESCRIPTOR,PARTICLETYPE>,
#ifdef PARALLEL_MODE_MPI
      communicate_surface_force<T,PARTICLETYPE>,
#endif
      apply_gravity<T,PARTICLETYPE>,
      process_dynamics<T,PARTICLETYPE>,
#ifdef PARALLEL_MODE_MPI
      update_particle_core_distribution<T, PARTICLETYPE>,
#endif
      couple_particles_to_lattice<T,DESCRIPTOR,PARTICLETYPE>
    >();

    // Get Results
    getResults(sLattice, converter, iT, superGeometry, timer, particleSystem );

    // Collide and stream
    sLattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();
}

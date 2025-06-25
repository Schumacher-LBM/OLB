/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2024 Pascal Sitter, Adrian Kummerlaender
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

/* cylinder2d.cpp:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 */

#include <olb.h>

#include <emscripten.h>
#include <emscripten/bind.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace emscripten;

using T          = double;
using DESCRIPTOR = D2Q9<>;

#define BOUZIDI

struct GeometryParameters {
  T L;
  T lengthX;
  T lengthY;
  T centerCylinderX;
  T centerCylinderY;
  T radiusCylinder;

  GeometryParameters(int N) {
    L = 0.1 / N;
    lengthX = 2.2;
    lengthY = 0.41 + L;
    centerCylinderX = 0.2 + L / 2.0;
    centerCylinderY = 0.2 + L / 2.0;
    radiusCylinder = 0.05;
  }
};

// Values needed to draw Image in JS
T*    values;
int   valuesLength;
int   numberCellsX;
int   numberCellsY;
float charPhysVelocity;
float physMaxU;
float avgRho;
float physTime = 0;
float mlups    = 0;
// Parameters for the simulation setup
bool    valuesIsVelocity = true;
int     N                = 10;      // resolution of the model
T       Re               = 20.;     // Reynolds number
T       maxPhysT         = 16;      // max. simulation time in s, SI unit
GeometryParameters geometryParameters(N);

//Global variables
std::unique_ptr<SuperLattice<T, DESCRIPTOR>> sLattice;
std::unique_ptr<SuperGeometry<T, 2>>         superGeometry;
std::unique_ptr<CuboidDecomposition<T,2>>    cuboidDecomposition;
std::unique_ptr<HeuristicLoadBalancer<T>>    loadBalancer;
std::unique_ptr<UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR>>
                                converter;
std::unique_ptr<util::Timer<T>> timer;

//------------------------------------------------------------------------------
// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter,
                     SuperGeometry<T, 2>&                superGeometry,
                     std::shared_ptr<IndicatorF2D<T>>    circle,
                     GeometryParameters geomParams)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T, 2> extend(geomParams.lengthX, geomParams.lengthY);
  Vector<T, 2> origin;

  superGeometry.rename(0, 2);

  superGeometry.rename(2, 1, {1, 1});

  // Set material number for inflow
  extend[0] = 2. * geomParams.L;
  origin[0] = -geomParams.L;
  IndicatorCuboid2D<T> inflow(extend, origin);
  superGeometry.rename(2, 3, 1, inflow);
  // Set material number for outflow
  origin[0] = geomParams.lengthX - geomParams.L;
  IndicatorCuboid2D<T> outflow(extend, origin);
  superGeometry.rename(2, 4, 1, outflow);
  // Set material number for cylinder
  superGeometry.rename(1, 5, circle);

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}
//------------------------------------------------------------------------------
// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>&        sLattice,
                    UnitConverter<T, DESCRIPTOR> const& converter,
                    SuperGeometry<T, 2>&                superGeometry,
                    std::shared_ptr<IndicatorF2D<T>>    circle)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator(1);
  sLattice.defineDynamics<SmagorinskyBGKdynamics>(bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  // Setting of the boundary conditions
  //if boundary conditions are chosen to be interpolated
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);


// Material=5 -->bouzidi / bounce back
#ifdef BOUZIDI
  setBouzidiBoundary(sLattice, superGeometry, 5, *circle);
#else
//  setBounceBackBoundary(sLattice, superGeometry, 5);
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 5);
#endif

  // Initial conditions
  AnalyticalConst2D<T, T> rhoF(1);
  std::vector<T>          velocity(2, T(0));
  AnalyticalConst2D<T, T> uF(velocity);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(bulkIndicator, rhoF, uF);
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}
//------------------------------------------------------------------------------
// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues(int iT, T L)
{

  OstreamManager clout(std::cout, "setBoundaryValues");

  // No of time steps for smooth start-up
  int iTmaxStart = (*converter).getLatticeTime(maxPhysT * 0.4);
  int iTupdate   = 5;

  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    // Smooth start curve, polynomial
    PolynomialStartScale<T, T> StartScale(iTmaxStart, T(1));

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T(iT)};
    T frac[1]  = {};
    StartScale(frac, iTvec);
    T maxVelocity   = (*converter).getCharLatticeVelocity() * 3. / 2. * frac[0];
    T distance2Wall = L / 2.;
    Poiseuille2D<T> poiseuilleU(*superGeometry, 3, maxVelocity, distance2Wall);

    (*sLattice).defineU(*superGeometry, 3, poiseuilleU);

    (*sLattice)
        .setProcessingContext<
            Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
            ProcessingContext::Simulation);
  }
}
//------------------------------------------------------------------------------
// Function to compute velocity values at core spatial locations
void computeVelocityValues(SuperLatticePhysVelocity2D<T, DESCRIPTOR>& velocity)
{
  (*sLattice).getBlock(0).forCoreSpatialLocations([&](LatticeR<2> latticeR) {
    T   v[DESCRIPTOR::d] {};
    int pos[3] = {0, latticeR[0], latticeR[1]};
    velocity(v, pos);
    // 2D to 1D Mapping of velocityVector Norms
    values[latticeR[0] + latticeR[1] * (*sLattice).getBlock(0).getNx()] =
        util::norm<2>(v);
  });
}
//------------------------------------------------------------------------------
// Function to compute pressure values at core spatial locations
void computePressureValues(SuperLatticePhysPressure2D<T, DESCRIPTOR>& pressure)
{
  (*sLattice).getBlock(0).forCoreSpatialLocations([&](LatticeR<2> latticeR) {
    T   p[DESCRIPTOR::d] {};
    int pos[3] = {0, latticeR[0], latticeR[1]};
    pressure(p, pos);
    // 2D to 1D Mapping of velocityVector Norms
    values[latticeR[0] + latticeR[1] * (*sLattice).getBlock(0).getNx()] =
        util::norm<2>(p);
  });
}
//------------------------------------------------------------------------------
// Computes the pressure drop between the voxels before and after the cylinder , util::Timer<T>& timer
void getResults(std::size_t iT, GeometryParameters geomParams)
{
  OstreamManager clout(std::cout, "getResults");

  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(*sLattice, *converter);
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(*sLattice, *converter);

  //
  // Two different methods so only one if statement is needed and not in the loop
  // Otherwise, they are the same apart from their functor.
  //
  if (valuesIsVelocity) {
    computeVelocityValues(velocity);
  }
  else {
    computePressureValues(pressure);
  }


  physTime = (*converter).getPhysTime(iT);
  physMaxU =
      (*converter).getPhysVelocity((*sLattice).getStatistics().getMaxU());
  // May be needed for convergence check
  avgRho = (*sLattice).getStatistics().getAverageRho();

  const int statIter = (*converter).getLatticeTime(.1);

  T point[2] = {};
  point[0]   = geomParams.centerCylinderX + 3 * geomParams.radiusCylinder;
  point[1]   = geomParams.centerCylinderY;
  AnalyticalFfromSuperF2D<T> intpolateP(pressure, true);
  T                          p;
  intpolateP(&p, point);

  if (iT % statIter == 0) {
    (*sLattice).setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    (*timer).update(iT);
    mlups = (*timer).getMLUPs();
    // timer.printStep();

    // Lattice statistics console output
    (*sLattice).getStatistics().print(iT, (*converter).getPhysTime(iT));

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF2D<T>            intpolatePressure(pressure, true);
    SuperLatticePhysDrag2D<T, DESCRIPTOR> drag(*sLattice, *superGeometry, 5,
                                               *converter);

    T point1[2] = {};
    T point2[2] = {};

    point1[0] = geomParams.centerCylinderX - geomParams.radiusCylinder;
    point1[1] = geomParams.centerCylinderY;

    point2[0] = geomParams.centerCylinderX + geomParams.radiusCylinder;
    point2[1] = geomParams.centerCylinderY;

    T p1, p2;
    intpolatePressure(&p1, point1);
    intpolatePressure(&p2, point2);

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1 - p2;
    clout << "; pressureDrop=" << pressureDrop;

    int input[3] = {};
    T   _drag[drag.getTargetDim()];
    drag(_drag, input);
    clout << "; drag=" << _drag[0] << "; lift=" << _drag[1] << std::endl;
  }
}
//------------------------------------------------------------------------------
extern "C" {
EMSCRIPTEN_KEEPALIVE
size_t setupCylinder2d()
{
  // === 1st Step: Initialization ===
  geometryParameters = GeometryParameters(N);

  converter.reset(new UnitConverterFromResolutionAndLatticeVelocity<T,
                                                                    DESCRIPTOR>(
      int {N}, // resolution: number of voxels per charPhysL
      (T)0.05, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T)2.0 *
          geometryParameters.radiusCylinder, // charPhysLength: reference length of simulation geometry
      (T)0.2, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)0.2 * 2. * geometryParameters.radiusCylinder /
          Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)1.0  // physDensity: physical density in __kg / m^3__
      ));

  // Prints the converter log as console output
  (*converter).print();

  // === 2rd Step: Prepare Geometry ===
  Vector<T, 2>         extend(geometryParameters.lengthX, geometryParameters.lengthY);
  Vector<T, 2>         origin;
  IndicatorCuboid2D<T> cuboid(extend, origin);

  // Instantiation of a cuboidDecomposition with weights
  cuboidDecomposition.reset(new CuboidDecomposition<T,2>(cuboid, geometryParameters.L, 1));

  // Instantiation of a loadBalancer
  loadBalancer.reset(new HeuristicLoadBalancer<T>(*cuboidDecomposition));

  // Instantiation of a superGeometry
  superGeometry.reset(new SuperGeometry<T, 2>(*cuboidDecomposition, *loadBalancer));

  Vector<T, 2>                     center(geometryParameters.centerCylinderX, geometryParameters.centerCylinderY);
  std::shared_ptr<IndicatorF2D<T>> circle =
      std::make_shared<IndicatorCircle2D<T>>(center, geometryParameters.radiusCylinder);

  prepareGeometry(*converter, *superGeometry, circle, geometryParameters);

  // === 3rd Step: Prepare Lattice ===
  sLattice.reset(new SuperLattice<T, DESCRIPTOR>(*superGeometry));

  //prepareLattice and set boundaryConditions
  prepareLattice(*sLattice, *converter, *superGeometry, circle);

  //initialization of velocityValues
  valuesLength = (*sLattice).getBlock(0).getNcells();
  numberCellsX = (*sLattice).getBlock(0).getNx();
  numberCellsY = (*sLattice).getBlock(0).getNy();
  values       = new T[valuesLength];

  charPhysVelocity = (*converter).getCharPhysVelocity();
  physMaxU         = 0;
  avgRho           = 0;
  physTime         = 0;
  mlups            = 0;

  timer.reset(new util::Timer<T>((*converter).getLatticeTime(maxPhysT),
                                 (*superGeometry).getStatistics().getNvoxel()));

  return (*converter).getLatticeTime(maxPhysT);
}
}
//------------------------------------------------------------------------------
extern "C" {
EMSCRIPTEN_KEEPALIVE
void runCylinder2d(size_t timeFrom, size_t timeTo)
{
  // === 4th Step: Main Loop with Timer ===
  for (std::size_t iT = timeFrom; iT < timeTo; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(iT, geometryParameters.L);

    // === 6th Step: Collide and Stream Execution ===
    (*sLattice).collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults(iT, geometryParameters);
  }
}
}
//------------------------------------------------------------------------------
extern "C" {
EMSCRIPTEN_KEEPALIVE
int getNumberCellsX() { return numberCellsX; }
}
extern "C" {
EMSCRIPTEN_KEEPALIVE
int getNumberCellsY() { return numberCellsY; }
}
extern "C" {
EMSCRIPTEN_KEEPALIVE
int getValuesCells() { return valuesLength; }
}
extern "C" {
EMSCRIPTEN_KEEPALIVE
float getCharPhysVelocity() { return charPhysVelocity; }
}
extern "C" {
EMSCRIPTEN_KEEPALIVE
float getPhysMaxU() { return physMaxU; }
}
extern "C" {
EMSCRIPTEN_KEEPALIVE
float getAvgRho() { return avgRho; }
}
extern "C" {
EMSCRIPTEN_KEEPALIVE
float getPhysTime() { return physTime; }
}
extern "C" {
EMSCRIPTEN_KEEPALIVE
float getMLUPs() { return mlups; }
}
//------------------------------------------------------------------------------
extern "C" {
EMSCRIPTEN_KEEPALIVE T* getValuesPointer() { return values; }
EMSCRIPTEN_KEEPALIVE void freeMemory()
{
  delete[] values;
  values = nullptr;
}
}
//------------------------------------------------------------------------------
extern "C" {
EMSCRIPTEN_KEEPALIVE void startTimer() { (*timer).start(); }
}
extern "C" {
EMSCRIPTEN_KEEPALIVE void stopTimer() { (*timer).stop(); }
}
//------------------------------------------------------------------------------
extern "C" {
EMSCRIPTEN_KEEPALIVE
void setSimulationParameters(T Re_, int N_)
{
  Re              = Re_;
  N               = N_;
  geometryParameters = GeometryParameters(N);
}
}
extern "C" {
EMSCRIPTEN_KEEPALIVE
void setOutputType(bool isVelocity) { valuesIsVelocity = isVelocity; }
}
//------------------------------------------------------------------------------
extern "C" {
EMSCRIPTEN_KEEPALIVE
int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  setupCylinder2d();
  OstreamManager clout(std::cout, "WASM");
  clout << "OpenLB Web Widget loaded." << std::endl;
}
}

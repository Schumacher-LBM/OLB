/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2006-2025 Fedor Bukreev, Shota Ito,
 *  Jonas Latt, Mathias J. Krause, Vojtech Cvrcek,
 *  Peter Weisbrod, Adrian Kummerländer
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

/* cylinder2d.h:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 */

#include <olb.h>

using namespace olb;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = descriptors::D2Q9<>;

// Parameters for the simulation setup
const int N = 10;       // resolution of the model
const T CFL = 0.05;     // characteristic CFL number
const T Re = 20.;       // Reynolds number
const T maxPhysT = 16;  // max. simulation time in s, SI unit
const T L = 0.1/N;      // latticeL
const T lengthX = 2.2;
const T lengthY = .41+L;
const T centerCylinderX = 0.2;
const T centerCylinderY = 0.2+L/2.;
const T radiusCylinder = 0.05;

// Stores geometry information in form of material numbers
void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,2>& superGeometry,
                     IndicatorF2D<T>& circle)
{
  Vector<T,2> extend(lengthX, lengthY);
  Vector<T,2> origin;
  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, {1,1});

  // Set material number for inflow
  extend[0] = 2.*L;
  origin[0] = -L;
  IndicatorCuboid2D<T> inflow(extend, origin);
  superGeometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = lengthX - L;
  IndicatorCuboid2D<T> outflow(extend, origin);
  superGeometry.rename(2, 4, 1, outflow);

  // Set material number for cylinder
  superGeometry.rename(1, 5, circle);

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();
  superGeometry.print();
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    const UnitConverter<T, DESCRIPTOR>& converter,
                    SuperGeometry<T,2>& superGeometry,
                    IndicatorF2D<T>& circle)
{
  // Material=1 -->bulk dynamics
  sLattice.defineDynamics<BGKdynamics>(superGeometry, 1);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  // Material=3 -->fixed velocity
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);

  // Material=4 -->fixed pressure
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  // Material=5 -->bouzidi
  setBouzidiBoundary(sLattice, superGeometry, 5, circle);

  // Initial conditions
  AnalyticalConst2D<T,T> rhoF(1);
  AnalyticalConst2D<T,T> uF(0, 0);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry, 1, rhoF, uF);
  sLattice.iniEquilibrium(superGeometry, 1, rhoF, uF);
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  sLattice.initialize();
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& sLattice,
                       const UnitConverter<T, DESCRIPTOR>& converter, std::size_t iT,
                       SuperGeometry<T,2>& superGeometry)
{
  // Number of time steps for smooth start-up
  const std::size_t iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
  const std::size_t iTupdate = iTmaxStart/1000;

  if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
    // Smooth start curve, polynomial
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getCharLatticeVelocity()*3./2.*frac[0];
    T distance2Wall = L/2.;
    Poiseuille2D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall );
    sLattice.defineU( superGeometry, 3, poiseuilleU );

    // Update velocity on GPU
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                const UnitConverter<T, DESCRIPTOR>& converter, std::size_t iT,
                SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer)
{
  OstreamManager clout( std::cout,"getResults" );
  const std::size_t vtkIter  = converter.getLatticeTime(0.3);
  const std::size_t statIter = converter.getLatticeTime(0.8);

  SuperVTMwriter2D<T> vtmWriter("cylinder2d");
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure(sLattice, converter);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the vtk files
  if (iT%vtkIter == 0 && iT > 0) {
    // Send values from GPU to CPU for evaluation
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  // Writes the console log
  if (iT%statIter == 0) {
    // Send values from GPU to CPU for evaluation
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT,converter.getPhysTime(iT));

    // Pressure drop
    AnalyticalFfromSuperF2D<T> intpolatePressure( pressure, true );

    T point1[2] = {};
    T point2[2] = {};

    point1[0] = centerCylinderX - radiusCylinder;
    point1[1] = centerCylinderY;

    point2[0] = centerCylinderX + radiusCylinder;
    point2[1] = centerCylinderY;

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop << std::endl;
  }

}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  // Set up the unit converter
  const UnitConverter<T,DESCRIPTOR> converter(
    (T)   L,                        // physDeltaX: spacing between two lattice cells in [m]
    (T)   CFL*L/0.2,                // physDeltaT: time step in [s]
    (T)   2.0*radiusCylinder,       // charPhysLength: reference length of simulation geometry in [m]
    (T)   0.2,                      // charPhysVelocity: highest expected velocity during simulation in [m/s]
    (T)   0.2*2.*radiusCylinder/Re, // physViscosity: physical kinematic viscosity in [m^2/s]
    (T)   1.0                       // physDensity: physical density in [kg/m^3]
  );
  converter.print();

  // === 2rd Step: Prepare Geometry ===
  Vector extend{lengthX, lengthY};
  Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extend, origin);
  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, L, singleton::mpi().getSize());

  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer);
  Vector center{centerCylinderX, centerCylinderY};
  IndicatorCircle2D<T> circle(center, radiusCylinder);
  prepareGeometry(converter, superGeometry, circle);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(sLattice, converter, superGeometry, circle);

  // === 4th Step: Main Loop with Timer ===
  std::size_t iTmax = converter.getLatticeTime(maxPhysT);
  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(sLattice, converter, iT, superGeometry);
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer);
  }

  timer.stop();
  timer.printSummary();
}

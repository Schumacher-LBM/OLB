/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006 - 2025 Mathias J. Krause, Jonas Fietz,
 *                            Jonas Latt, Jonas Kratzke, Shota Ito
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

/* cavity2d.cpp:
 * This example illustrates a minimal working example for a
 * fluid simulation with OpenLB; a flow in a cuboid, the lid-
 * driven cavity.
 */

#include <olb.h>
#include </home/clara/OpenLB/release-1.8.0/examples/noise/noiseauxiliary/acousticPulse.h>

using namespace olb;

using T = FLOATING_POINT_TYPE;

using DESCRIPTOR = descriptors::D2Q9<>;

using BulkDynamics = ConstRhoBGKdynamics<T,DESCRIPTOR>;

const T physDeltaX        = 0.0078125;   // grid spacing [m]
const T physDeltaT        = 0.00078125;  // temporal spacing [s]
const T physLength        = 1.0;         // length of the squared cuboid [m]
const T physLidVelocity   = 1.0;         // velocity imposed on lid [m/s]
const T physViscosity     = 0.001;        // kinetic viscosity of fluid [m*m/s]
const T physDensity       = 1.0;         // fluid density [kg/(m*m*m)]
const T physMaxT          = 0.5;        // maximal simulation time [s]

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,2>& superGeometry)
{
  // Set material numbers to assign physics to lattice nodes
  superGeometry.rename(0, 1);
  //superGeometry.rename(2, 1, {1,1});

  T dx = converter.getPhysDeltaX();
  Vector extend{physLength + 2*dx, 2*dx};
  Vector origin{-dx, physLength - dx};
  //IndicatorCuboid2D lid(extend, origin);
 // superGeometry.rename(2,3,1,lid);

  superGeometry.getStatistics().print();
}

void prepareLattice(const UnitConverter<T,DESCRIPTOR>& converter,
                    SuperLattice<T, DESCRIPTOR>& sLattice,
                    SuperGeometry<T,2>& superGeometry)
{
  // Material=1 -->bulk dynamics
  sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}



//Original
// void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
//                        SuperLattice<T, DESCRIPTOR>& sLattice,
//                        std::size_t iT, SuperGeometry<T,2>& superGeometry)
// {
//   if (iT == 0) {
//     auto domain = superGeometry.getMaterialIndicator({1,2,3});
//    //     AnalyticalConst2D<T,T> rhoF(1);  // Definiert rho auf 1
//     AnalyticalConst2D<T,T> uWall(0, 0);   // Definiert die Geschwindigkeit überall auf 0 in X und Y Richtung
//     AnalyticalConst2D<T,T> uLid(converter.getCharLatticeVelocity(), 0); // Definiert die Geschwindigkeit an der lid auf eine Geschwindigkeit in X Richtung und 0 in Y Richtung

//     // Initialize populations to equilibrium state
//     sLattice.iniEquilibrium(domain, rhoF, uWall);
//     sLattice.defineRhoU(domain, rhoF, uWall);

//     // Assign the non-zero velocity to the lid
//     sLattice.defineU(superGeometry, 3, uLid);
//     sLattice.initialize();
//   }
// }


// Hier werden die Startbedingungen für eine Gausschen Puls definiert
void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
  SuperLattice<T, DESCRIPTOR>& sLattice,
  std::size_t iT, SuperGeometry<T,2>& superGeometry)
{
if (iT == 0) {   // Der Puls wird zum Zeitpunkt t=0 definiert

auto domain = superGeometry.getMaterialIndicator({1,2,3}); // Die Materialien werden Fluid (1), gleitfreie Grenzen (2) und Geschiwindigkeit gesetzt (3) Quelle: Handbuch

// Dichtefeld mit höherer Dichte im Zentrum (z. B. als Gauß-Glocke)
Vector<T,2> center(0.5, 0.5); // Zentrum der Welle

//T radius = 0.05; // Wird momentan noch nicht genutzt
// Gaußförmige Dichteverteilung
AcousticPulse<2,T> densityWave(
2., 1e-3, 120., center);           // Rho an einer Stelle auf 1.5 definiert
AnalyticalConst2D<T,T> rhoF(1);       // Definiert rho auf 1
AnalyticalConst2D<T,T> uWall(0, 0);   // Definiert die Geschwindigkeit überall auf 0 in X und Y Richtung
AnalyticalConst2D<T,T> uLid(0, 0);    // Definiert die Geschwindigkeit an der lid auf eine Geschwindigkeit in X Richtung und 0 in Y Richtung

// Initialize populations to equilibrium state
sLattice.iniEquilibrium(domain, densityWave, uWall);
sLattice.defineRhoU(domain, densityWave, uLid);


// Assign the non-zero velocity to the lid
// sLattice.defineRho(superGeometry,densityWave,uWall);
sLattice.defineU(superGeometry, 3, uLid);
sLattice.initialize();
}
}


void getResults(const UnitConverter<T,DESCRIPTOR>& converter,
                SuperLattice<T, DESCRIPTOR>& sLattice,
                std::size_t iT, util::Timer<T> timer)
{
  const std::size_t iTvtk = converter.getLatticeTime(physMaxT/100.);  // Ausgabe von 100 Zeitschritten. physMaxT ist eine reale Zeiteinheit; getLatticeTime konvertiert in die Lattice-Zeiteinheit
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/20.);

  SuperVTMwriter2D<T> vtmWriter("cavity2d");
  SuperLatticePhysVelocity2D velocity(sLattice, converter);
  SuperLatticePhysPressure2D pressure(sLattice, converter);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT%iTvtk == 0 && iT >= 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  // Get statistics
  if (iT%iTlog == 0) {
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout( std::cout,"main" ); // writing all output first in a userdefined Buffer of type OMBuf. On a flush it spits out at first the userdefined text in squared brackets and afterwards everything from the buffer

  // Provide the unit converter the characteristic entities
  const UnitConverter<T,DESCRIPTOR> converter (
    physDeltaX,        // physDeltaX: spacing between two lattice cells in [m]
    physDeltaT,        // physDeltaT: time step in [s]
    physLength,        // charPhysLength: reference length of simulation geometry in [m]
    physLidVelocity,   // charPhysVelocity: highest expected velocity during simulation in [m/s]
    physViscosity,     // physViscosity: physical kinematic viscosity in [m^2/s]
    physDensity        // physDensity: physical density [kg/m^3]
  );
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  Vector extend{physLength, physLength};
  Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extend, origin);
  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());
  cuboidDecomposition.setPeriodicity({true,true});
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer);
  prepareGeometry(converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(converter, sLattice, superGeometry);

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
    getResults(converter, sLattice, iT, timer);
  }

  timer.stop();
  timer.printSummary();
}

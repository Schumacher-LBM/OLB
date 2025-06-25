/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender
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

#include <olb.h>

using namespace olb;

int main(int argc, char **argv) {
  initialize(argc, argv);

  using T = double;
  using DESCRIPTOR = descriptors::D3Q19<>;

  const T deltaX = 1e-2;
  const T deltaT = 1e-3;
  UnitConverter<T,DESCRIPTOR> converter(deltaX, deltaT, 1, 1, 1e-2, 1);
  converter.print();

  IndicatorCuboid3D<T> domainI({1,1,1}, {0,0,0});
  // Decompose periodic domain into cuboids
  CuboidDecomposition<T,3> cDecomposition(domainI, converter.getPhysDeltaX(), singleton::mpi().getSize());
  cDecomposition.setPeriodicity(true);

  // Assign cuboids to available parallel processing units
  HeuristicLoadBalancer balancer(cDecomposition);

  // Set up lattice structure
  SuperLattice<T,DESCRIPTOR> sLattice(cDecomposition, balancer);

  // Use BGK collision step everywhere (default second order equilibrium)
  sLattice.defineDynamics<BGKdynamics>();
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  // Perform a single collide and stream cycle
  sLattice.collideAndStream();

  // Output macroscopic velocity as VTK
  SuperVTMwriter3D<T> writer("minimal3d");
  SuperLatticePhysVelocity3D velocityF(sLattice, converter);
  writer.write(velocityF);
}

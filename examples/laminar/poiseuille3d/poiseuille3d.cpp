/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Marc Haußmann, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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

/* poiseuille3d.cpp:
 * This example examines a 3D Poseuille flow
 * It illustrates the usage of different flow setups, boundary conditions
 * and computation of error norms.
 * The following simulation options can be combined freely:
 * - if the compiler flag ENABLE_MRT is set, mrt collision operators are used
 * - forced/ nonForced flow
 * - different boundary conditions
 * - simulation only or eoc convergence analysis
 */

// the main code of the simulation is in poiseuille3d.h as it is also used by the
// example ../../pdeSolverEoc/poiseuille3dEoc

// set flag in order to use mrt collision operators instead of bgk
//#define ENABLE_MRT

#include "poiseuille3d.h"

//Initialize Gnuplot
static Gnuplot<T> gplot("centerVelocity");

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  int N = 21;
  bool eoc = false;

  int status = readParameters(argc, argv, N, flowType, boundaryType);
  if (status != 0) {
    return status;
  }
  simulatePoiseuille(N, gplot, eoc);

}

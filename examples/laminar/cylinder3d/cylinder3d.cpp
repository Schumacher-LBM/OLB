/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2013 Mathias J. Krause, Thomas Henn, Tim Dornieden
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

/* cylinder3d.cpp:
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
 * It also shows the usage of the STL-reader and explains how
 * to set boundary conditions automatically.
 */


#include "cylinder3d.h"

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  /// Simulation Parameter
  int N = 10;

  bool eoc = false;

  static Gnuplot<T> gplot( "drag" );

  simulateCylinder(N, gplot, eoc);

  return 0;
}

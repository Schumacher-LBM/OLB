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
// example ../../laminar/poiseuille3d

// set flag in order to use mrt collision operators instead of bgk
//#define ENABLE_MRT

#include "../../laminar/poiseuille3d/poiseuille3d.h"

/// Initialize gnuplot
  static Gnuplot<T> gplot(
    "Velocity_and_StrainRate_eoc",
    false,
    "set terminal png size 720, 720 font 'Arial,10'",
    Gnuplot<T>::LOGLOG,
    Gnuplot<T>::LINREG);

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  int startN = 21;
  int maxN = startN + 31;
  bool eoc = true;

  int status = readParameters(argc, argv, startN, flowType, boundaryType);
  if (status != 0) {
    return status;
  }

  if ((boundaryType == freeSlip) || (boundaryType == partialSlip)) {
    throw std::invalid_argument(
      "eoc computation is currently not supported for slip boundary conditions");
  }

  // set the labels for the plot
  gplot.setLabel("Resolution test", "average Error");

  // loop over the different simulations
  for(int simuN = startN; simuN < maxN; simuN += 10){

    /// Run the simulations
    clout << "Starting next simulation with N = " << simuN << std::endl;
    simulatePoiseuille(simuN, gplot, eoc);
  }

  gplot.writePNG();

  return 0;
}

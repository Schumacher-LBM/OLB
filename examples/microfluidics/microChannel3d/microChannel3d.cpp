/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Adrian Kummerländer,
 *  Mathias J. Krause
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
// This is an 3D example for a laminar flow in a round pipe
// at higher Knudsen numbers with Knudsen slip condition and
// thermal creep along the wall

#include "microChannel3d.h"
/// Initialize gnuplot
static Gnuplot<T> gplot(
  "Velocity_eoc",
  false,
  "set terminal png size 720, 720 font 'Arial,10'",
  Gnuplot<T>::LOGLOG,
  Gnuplot<T>::LINREG);

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  CLIreader args(argc, argv);
  const T relTime = args.getValueOrFallback("--relTime", 1.4);

  /// Simulation Parameter
  int startN = 41;
  int maxN = 162;

  // set the labels for the plot
  gplot.setLabel("Resolution test", "average Error");

  // loop over the different simulations
  for(int simuN = startN; simuN < maxN; simuN = simuN*2-1){

    /// Run the simulations
    clout << "Starting next simulation with N = " << simuN << std::endl;
    simulateMicroChannel3D(simuN, gplot, relTime);
  }

  gplot.writePNG();

  return 0;
}

/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Philipp Spelten
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

/* linePlot.h:
 * Function creating and extracting a line plot of a given
 * AnalyticalF3D along x-axis, y-axis, x-y-diagonal, or
 * x-y-z-diagonal.
*/

#ifndef LINE_PLOT_H
#define LINE_PLOT_H

namespace olb {

typedef enum { horizontal, vertical, diagonal2d, diagonal3d } SamplingDirection;

template <unsigned int ndim, typename T>
void linePlot(AnalyticalF3D<T, T>& data, size_t ndatapoints, T dist, std::string title, std::string ylabel,
              SamplingDirection direction, bool halfDomain = true, bool setRange = false, T ymin = 0, T ymax = 0)
{
  Gnuplot<T> gplot(title);
  gplot.setLabel("distance [m]", ylabel);
  int nmin = 0;
  if (!halfDomain)
    nmin = -int(ndatapoints / 2);
  for (int n = nmin; n <= int(ndatapoints / 2); n++) {
    T input[ndim] = {0, 0, 0};
    T distance    = 0;
    switch (direction) {
    case horizontal:
      input[0] = n * dist;
      distance = n * dist;
      break;
    case vertical:
      input[1] = n * dist;
      distance = n * dist;
      break;
    case diagonal2d:
      input[0] = n * dist;
      input[1] = n * dist;
      distance = n * dist * std::sqrt(2);
      break;
    case diagonal3d:
      input[0] = n * dist;
      input[1] = n * dist;
      input[2] = n * dist;
      distance = n * dist * std::sqrt(3);
      break;
    }
    T output[1];
    data(output, input);
    gplot.setData(distance, output[0]);
  }
  if (setRange)
    gplot.setYrange(ymin, ymax);
  gplot.writePNG(-1, -1, title);
};

} // namespace olb

#endif

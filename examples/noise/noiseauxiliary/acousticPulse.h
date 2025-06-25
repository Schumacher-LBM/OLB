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

/* acousticPulse.h:
 * Function returning a gaussian pressure field described by
 * rho(x) = rho0 + amplitude*e^(-alpha*||x-x0||)
*/

#ifndef ACOUSTIC_PULSE_H
#define ACOUSTIC_PULSE_H

namespace olb {

template <unsigned ndim, typename T>
class AcousticPulse : public AnalyticalF<ndim, T, T> {
protected:
  T               rho0;
  T               amplitude;
  T               alpha;
  Vector<T, ndim> x0;

public:
  AcousticPulse(T rho0, T amplitude, T alpha, Vector<T, ndim> x0 = Vector<T, ndim>(0.))
      : AnalyticalF<ndim, T, T>(1)
      , rho0(rho0)
      , amplitude(amplitude)
      , alpha(alpha)
      , x0(x0) {};

  bool operator()(T output[], const T input[]) override
  {
    T distance = 0;
    for (unsigned d = 0; d < ndim; d++) {
      T distance_i = input[d] - x0[d];
      distance += distance_i * distance_i; // actually, distance would be square root, but would be squared in exponent
    }
    output[0] = rho0 + amplitude * util::exp(-alpha * distance);
    return true;
  };
};

} // namespace olb

#endif
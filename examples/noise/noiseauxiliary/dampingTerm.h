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

/* dampingTerm.h:
 * The newtonian cooling term suggested in M. Israeli, S.A. Orszag,
 * Approximation of radiation boundary conditions,
 * J. Comput. Phys. 41 (1981) 115–135.
*/

#ifndef DAMPING_TERM_H
#define DAMPING_TERM_H

namespace olb {

template <unsigned ndim, typename T, typename DESCRIPTOR>
class DampingTerm : public AnalyticalF<ndim, T, T> {
protected:
  Vector<T, ndim> x_0, Lx;
  int             A = 4;
  int             n = 4;
  T               _bdPU;
  T               _ds;

public:
  DampingTerm(T boundaryDepthPU, Vector<T, ndim> domain_lengths, T dampingStrength = 1.)
      : AnalyticalF<ndim, T, T>(1)
      , _bdPU(boundaryDepthPU)
      , _ds(dampingStrength)
  {
    for (size_t d = 0; d < ndim; d++) {
      Lx[d]  = domain_lengths[d] / 2; // Lx is usually half the domain
      x_0[d] = Lx[d] - _bdPU;         // subtract boundary depth from domain length
      // x_0[d] /= Lx[d];  // normalize x_0 (Lx is still needed in operator() to normalize x; Lx will be replaced by 1)
    }
  };

  bool operator()(T output[], const T input[]) override
  {
    // normalize x to [0,1]
    Vector<T, ndim> x, distance_from_border;
    bool            is_boundary = false;
    for (size_t d = 0; d < ndim; d++) {
      x[d] = (std::abs(input[d]) - x_0[d]) / _bdPU; // x_0 becomes 0
      x[d] = std::max(x[d], T(0)); // if x[d]<0, ignore for X; if x[d]>1, set to one to avoid negative values of sigma
      if (x[d] > 0) {
        is_boundary = true;
      }
      distance_from_border[d] = Lx[d] - std::abs(input[d]);
    }
    if (!is_boundary) {
      output[0] = 0;
      return true;
    }

    T X = 0.;
    for (size_t d = 0; d < ndim; d++) {
      if (distance_from_border[d] ==
          *std::min_element(std::begin(distance_from_border), std::end(distance_from_border))) {
        X = x[d];
      }
    }
    T sigma;
    sigma = A * (((std::pow(X, n)) * (1 - X) * (n + 1) * (n + 2)) / std::pow(_bdPU, n + 2)); // x_0 left out
    sigma /= 9.8304 / std::pow(_bdPU, 6); // calculated as max of function (at X=0.8) for A=n=4
    sigma = std::max(sigma, T(0));
    sigma = std::min(sigma, T(1)); // hard set 0<sigma<1
    sigma *= _ds;
    output[0] = sigma;
    return true;
  };
};

} // namespace olb

#endif

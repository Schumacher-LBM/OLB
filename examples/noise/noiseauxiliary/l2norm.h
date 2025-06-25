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

/* l2norm.h:
 * Function returning the pressure L2-norm of a given SuperLattice
*/

#ifndef L2_NORM_H
#define L2_NORM_H

namespace olb {

template <unsigned int ndim, typename T, typename DESCRIPTOR>
T L2Norm(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T, ndim>& sGeometry,
         UnitConverter<T, DESCRIPTOR> const& converter, int fluidMaterial)
{
  OstreamManager                            clout(std::cout, "L2Norm");
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressurenD(sLattice, converter);
  AnalyticalConst<ndim, T, T>               rho0nD(0.);
  SuperAbsoluteErrorL2Norm3D<T>             absPressureErrorL2Norm(pressurenD, rho0nD,
                                                                   sGeometry.getMaterialIndicator(fluidMaterial));
  T                                         result[ndim];
  int                                       tmp[] = {int()};
  absPressureErrorL2Norm(result, tmp);
  return result[0];
};

} // namespace olb

#endif
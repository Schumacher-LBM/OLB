/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Stephan Simonis
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

#pragma once

namespace olb {

template <typename T, typename _DESCRIPTOR>
class FringeZoneSmagorinskyConstant : public AnalyticalF3D<T,T> {
private:
  const T csBulk;
  const T charPhysLength;
  const T physDeltaX;

  T s(T x) {
    if ( x < std::numeric_limits<T>::epsilon() ) {
      return 0.;
    }
    else if (x > T(1)-std::numeric_limits<T>::epsilon() ) {
      return T(1);
    }
    else {
      return T(1)/(T(1) + util::exp( T(1)/(x - T(1) ) + T(1)/x ) );
    }
  }

public:
  FringeZoneSmagorinskyConstant(const UnitConverter<T,_DESCRIPTOR>& converter,
                                T bulkSmagorinskyConstant)
    : AnalyticalF3D<T,T>(1)
    , csBulk{bulkSmagorinskyConstant}
    , charPhysLength{converter.getCharPhysLength()}
    , physDeltaX{converter.getPhysDeltaX()}
  { }

  bool operator()(T output[], const T input[]) override {
    const T x = input[0];

    T empiricalK = 200.;                       // increase if required Re=10^4 -> 200, Re 2*10^4 - 300
    T xStart = 45.*charPhysLength-physDeltaX;  // T xStart = 55.*charPhysLength-physDeltaX;
    T xEnd = 70.*charPhysLength;               // 60.*charPhysLength-2.*physDeltaX; // if inside domain
    T xRise = 15.*charPhysLength;              // T xRise = 5.*charPhysLength;
    T xFall = 2.*charPhysLength;

    // OPTION A:
    /* Adaption of fringe zone [Lundbladh et al. 1999, url: http://www.fluidosol.se/thesismod/paper9.pdf]
     *
     *               _____________ ................. empiricalK*csBulk
     *              /             \
     *             /               \
     *            /                 \
     *           /                   \
     *          /                     \
     * ________/                       \________ ... csBulk
     *         |xStart                 |xEnd
     *         |-----|           |-----|
     *          xRise             xFall
     *
     * Fringe function:
     *  empiricalK*csBulk   (maximal Smagorinsky constant)
     *  xStart              (begin of the fringe zone)
     *  xEnd                (end of the fringe zone)
     *  bRise               (rise distance)
     *  bFall               (fall distance)
     *
     * S is a smooth step function:
     *   S(x)=0,                                         for x<=0,
     *   S(x)=1/( 1 + util::exp( (1/(x-1)) + (1/x) ) ),  for 0<x<1,
     *   S(x)=1,                                         for x>=1.
     */
    output[0] = csBulk*(empiricalK*(s((x - xStart)/xRise ) - s((x - xEnd)/xFall + T(1.)) ) + 1.0);

    // OPTION B:
    /* Adaption of sponge zone [Xue et al. 2022, doi: 10.1063/5.0090641]
     *
     *  \nu_{sponge} = \nu_{eff} [ K ( (x - x_{start})/(x_{end}-x_{start}) )^{p} + 1 ],
     * where:
     *  K = 1000  (empirical)
     *  p = 3     (empirical)
     *  \nu_{eff} (as always, effective viscosity of LBM LES scheme)
     *  x_{start} (starting point of sponge layer)
     *  x_{end}   (end point of sponge layer)
     *  x         (current grid point of evaluation)
     *
     * Sponge Smagorinsky constant as function of spatial location in x.
     * The constant effectuates a change of viscosity by a factor of order O(K(x-...)^{p}+1).
     */
    // output[0] = csBulk * ( empiricalK * util::pow((x - xStart)/(xEnd-xStart), empiricalP) + 1 );

    return true;
  };
};

}

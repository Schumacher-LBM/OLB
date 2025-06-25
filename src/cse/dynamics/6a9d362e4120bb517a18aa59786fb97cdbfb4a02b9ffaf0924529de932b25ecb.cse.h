/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021-24 Adrian Kummerlaender, Shota Ito
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

/*  ========================================================
 *  ==  WARNING: This is an automatically generated file, ==
 *  ==                  do not modify.                    ==
 *  ========================================================
 */

#pragma once


namespace olb {

namespace dynamics {

template <typename T, typename... FIELDS>
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::FreeEnergyInletOutletDensity, momenta::FreeEnergyInletOutletMomentum<0, -1>, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately>, equilibria::FreeEnergy, collision::FreeEnergyInletOutlet<0, -1>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0]*x10;
auto x12 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x13 = cell.template getFieldComponent<descriptors::FORCE>(0)*(V{1.5}*x12 + V{-1}) + V{1};
auto x14 = cell[1]*x10;
auto x15 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::FORCE>(0)*x9;
auto x16 = cell.template getFieldComponent<descriptors::FORCE>(0)*(-V{3}*cell.template getFieldComponent<descriptors::FORCE>(1) + V{3}*x12 + V{1}) + V{-1};
auto x17 = x15 - V{0.0277777777777778}*x16*x9;
auto x18 = cell[2]*x10;
auto x19 = V{0.111111111111111}*x9;
auto x20 = cell.template getFieldComponent<descriptors::FORCE>(0)*x19;
auto x21 = cell[3]*x10;
auto x22 = cell[4]*x10;
auto x23 = x13*x19 + x20;
auto x24 = cell[8]*x10;
auto x25 = V{0.111111111111111}*x9;
auto x26 = V{0.166666666666667}*cell.template getFieldComponent<descriptors::FORCE>(0) + V{0.166666666666667}*x11 + x13*x25 + V{0.166666666666667}*x14 - x15 - V{0.0277777777777778}*x16*x9 + V{0.166666666666667}*x18 + V{0.166666666666667}*x21 + V{0.166666666666667}*x22 + V{0.166666666666667}*x24 + V{-0.166666666666667};
cell[0] = V{0.555555555555556}*cell.template getFieldComponent<descriptors::FORCE>(0)*x9 - x11 - V{0.444444444444444}*x13*x9;
cell[1] = -x14 - x17;
cell[2] = V{0.111111111111111}*x16*x9 - x18 - x20;
cell[3] = -x17 - x21;
cell[4] = -x22 - x23;
cell[5] = x26;
cell[6] = V{0.666666666666667}*cell.template getFieldComponent<descriptors::FORCE>(0) + V{0.666666666666667}*x11 + V{0.444444444444444}*x13*x9 + V{0.666666666666667}*x14 - x16*x25 + V{0.666666666666667}*x18 - x20 + V{0.666666666666667}*x21 + V{0.666666666666667}*x22 + V{0.666666666666667}*x24 + V{-0.666666666666667};
cell[7] = x26;
cell[8] = -x23 - x24;
return { cell.template getFieldComponent<descriptors::FORCE>(0), x12 };
}
};

}

}

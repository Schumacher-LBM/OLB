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
struct CSE<dynamics::Tuple<T, descriptors::D3Q7<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::BGK, AdvectionDiffusionExternalVelocityCollision>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x7 = parameters.template get<descriptors::OMEGA>();
auto x8 = x7 + V{-1};
auto x9 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6];
auto x10 = V{4}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x11 = x9 + V{1};
auto x12 = V{0.125}*x7;
auto x13 = V{4}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x14 = V{4}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x0 = -cell[0]*x8 + V{0.25}*x7*x9;
auto x1 = -cell[1]*x8 - x12*(x11*(x10 + V{-1}) + V{1});
auto x2 = -cell[2]*x8 - x12*(x11*(x13 + V{-1}) + V{1});
auto x3 = -cell[3]*x8 - x12*(x11*(x14 + V{-1}) + V{1});
auto x4 = -cell[4]*x8 + V{0.125}*x7*(x11*(x10 + V{1}) + V{-1});
auto x5 = -cell[5]*x8 + V{0.125}*x7*(x11*(x13 + V{1}) + V{-1});
auto x6 = -cell[6]*x8 + V{0.125}*x7*(x11*(x14 + V{1}) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
return { x9 + V{1}, cell.template getFieldComponent<descriptors::VELOCITY>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1) + cell.template getFieldComponent<descriptors::VELOCITY>(2)*cell.template getFieldComponent<descriptors::VELOCITY>(2) };
}
};

}

}

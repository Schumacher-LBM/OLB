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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x12 = x11 + V{1};
auto x13 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x14 = V{1.5}*x13;
auto x15 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x16 = V{1.5}*x15;
auto x17 = x16 + V{-1};
auto x18 = x14 + x17;
auto x19 = V{0.0277777777777778}*x9;
auto x20 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x22 = -x21;
auto x23 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x24 = x18 + x23;
auto x25 = V{0.111111111111111}*x9;
auto x26 = -x23;
auto x27 = V{1} - x14;
auto x28 = V{3}*x15 + x27;
auto x29 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x30 = V{4.5}*(x29*x29);
auto x31 = V{3}*x13;
auto x32 = -x16 + x20;
auto x0 = -cell[0]*x10 - V{0.444444444444444}*x9*(x12*x18 + V{1});
auto x1 = -(cell[1]*x10 + x19*(x12*(-x20 + x24 - V{4.5}*x22*x22) + V{1}));
auto x2 = -cell[2]*x10 + x25*(x12*(x26 + x28) + V{-1});
auto x3 = -cell[3]*x10 - x19*(x12*(x20 + x24 - x30) + V{1});
auto x4 = -cell[4]*x10 - x25*(x12*(x17 + x20 - x31) + V{1});
auto x5 = -(cell[5]*x10 + x19*(x12*(x18 + x20 + x26 - V{4.5}*x21*x21) + V{1}));
auto x6 = -cell[6]*x10 + x25*(x12*(x23 + x28) + V{-1});
auto x7 = -cell[7]*x10 + x19*(x12*(x23 + x27 + x30 + x32) + V{-1});
auto x8 = -cell[8]*x10 + x25*(x12*(x31 + x32 + V{1}) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x11 + V{1}, x13 + x15 };
}
};

}

}

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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = parameters.template get<statistics::AVERAGE_RHO>();
auto x11 = x9 + V{-1};
auto x12 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x13 = V{1.5}*x12;
auto x14 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x15 = V{1.5}*x14;
auto x16 = x15 + V{-1};
auto x17 = x13 + x16;
auto x18 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x19 = -(x10 + V{-1})/(x18 + V{1}) + V{1};
auto x20 = x18 + V{1};
auto x21 = x19*x20;
auto x22 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x23 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x24 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x25 = -x24;
auto x26 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x27 = x17 + x26;
auto x28 = -x23 + x27 - V{4.5}*x25*x25;
auto x29 = V{0.0277777777777778}*x21;
auto x30 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x31 = -x26;
auto x32 = V{1} - x13;
auto x33 = V{3}*x14 + x32;
auto x34 = x31 + x33;
auto x35 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x36 = V{4.5}*(x35*x35);
auto x37 = x23 + x27 - x36;
auto x38 = V{3}*x12;
auto x39 = x16 + x23 - x38;
auto x40 = x17 + x23 + x31 - V{4.5}*x24*x24;
auto x41 = x26 + x33;
auto x42 = -x15 + x23;
auto x43 = x26 + x32 + x36 + x42;
auto x44 = x38 + x42 + V{1};
auto x0 = -x11*(cell[0] + x17*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) + V{0.444444444444444}) - V{0.444444444444444}*x17*x21 + V{-0.444444444444444};
auto x1 = -x11*(cell[1] + x22*x28 + V{0.0277777777777778}) - x28*x29 + V{-0.0277777777777778};
auto x2 = -x11*(cell[2] - x30*x34 + V{0.111111111111111}) + V{0.111111111111111}*x19*x20*x34 + V{-0.111111111111111};
auto x3 = -x11*(cell[3] + x22*x37 + V{0.0277777777777778}) - x29*x37 + V{-0.0277777777777778};
auto x4 = -x11*(cell[4] + x30*x39 + V{0.111111111111111}) - V{0.111111111111111}*x21*x39 + V{-0.111111111111111};
auto x5 = -x11*(cell[5] + x22*x40 + V{0.0277777777777778}) - x29*x40 + V{-0.0277777777777778};
auto x6 = -x11*(cell[6] - x30*x41 + V{0.111111111111111}) + V{0.111111111111111}*x19*x20*x41 + V{-0.111111111111111};
auto x7 = -x11*(cell[7] - x22*x43 + V{0.0277777777777778}) + V{0.0277777777777778}*x19*x20*x43 + V{-0.0277777777777778};
auto x8 = -x11*(cell[8] - x30*x44 + V{0.111111111111111}) + V{0.111111111111111}*x19*x20*x44 + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -x10 + x18 + V{2}, x12 + x14 };
}
};

}

}

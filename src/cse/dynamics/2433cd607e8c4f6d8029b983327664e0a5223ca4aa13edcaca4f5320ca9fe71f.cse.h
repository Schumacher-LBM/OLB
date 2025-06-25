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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = parameters.template get<statistics::AVERAGE_RHO>();
auto x11 = x9 + V{-1};
auto x12 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1};
auto x13 = V{1} / (x12);
auto x14 = V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + V{1};
auto x15 = cell[0] + cell[4] + cell[8] + x14;
auto x16 = x13*x15;
auto x17 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x18 = V{1.5}*x17;
auto x19 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x20 = V{1.5}*x19;
auto x21 = x20 + V{-1};
auto x22 = x18 + x21;
auto x23 = x10 + V{-1};
auto x24 = x12*x23/x15 + V{1};
auto x25 = V{0.0277777777777778}*x16;
auto x26 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x27 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x28 = -x27;
auto x29 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x30 = x22 + x29;
auto x31 = -x26 + x30 - V{4.5}*x28*x28;
auto x32 = V{0.111111111111111}*x16;
auto x33 = -x29;
auto x34 = V{1} - x18;
auto x35 = V{3}*x19 + x34;
auto x36 = x33 + x35;
auto x37 = x24*x32;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x39 = V{4.5}*(x38*x38);
auto x40 = x26 + x30 - x39;
auto x41 = V{3}*x17;
auto x42 = x21 + x26 - x41;
auto x43 = x22 + x26 + x33 - V{4.5}*x27*x27;
auto x44 = x29 + x35;
auto x45 = -x20 + x26;
auto x46 = x29 + x34 + x39 + x45;
auto x47 = x41 + x45 + V{1};
auto x0 = -x11*(cell[0] - V{0.444444444444444}*x16*x22 + V{0.444444444444444}) + V{0.444444444444444}*x13*x15*x22*x24 + V{-0.444444444444444};
auto x1 = -x11*(cell[1] - x25*x31 + V{0.0277777777777778}) + V{0.0277777777777778}*x13*x15*x24*x31 + V{-0.0277777777777778};
auto x2 = -x11*(cell[2] + x32*x36 + V{0.111111111111111}) - x36*x37 + V{-0.111111111111111};
auto x3 = -x11*(cell[3] - x25*x40 + V{0.0277777777777778}) + V{0.0277777777777778}*x13*x15*x24*x40 + V{-0.0277777777777778};
auto x4 = -x11*(cell[4] - x32*x42 + V{0.111111111111111}) + V{0.111111111111111}*x13*x15*x24*x42 + V{-0.111111111111111};
auto x5 = -x11*(cell[5] - x25*x43 + V{0.0277777777777778}) + V{0.0277777777777778}*x13*x15*x24*x43 + V{-0.0277777777777778};
auto x6 = -x11*(cell[6] + x32*x44 + V{0.111111111111111}) - x37*x44 + V{-0.111111111111111};
auto x7 = -x11*(cell[7] + x25*x46 + V{0.0277777777777778}) - x24*x25*x46 + V{-0.0277777777777778};
auto x8 = -x11*(cell[8] + x32*x47 + V{0.111111111111111}) - x37*x47 + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -x13*(V{1}*cell[0] + V{1}*cell[4] + V{1}*cell[8] + x14) - x23, x17 + x19 };
}
};

}

}

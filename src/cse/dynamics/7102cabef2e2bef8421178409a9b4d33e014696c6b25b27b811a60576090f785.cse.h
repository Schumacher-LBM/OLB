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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = parameters.template get<statistics::AVERAGE_RHO>();
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x11 = x9 + V{-1};
auto x12 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1};
auto x13 = V{1} / (x12);
auto x14 = V{2}*cell[1] + V{2}*cell[7] + V{2}*cell[8] + V{1};
auto x15 = cell[0] + cell[2] + cell[6] + x14;
auto x16 = x13*x15;
auto x17 = V{0.444444444444444}*x16;
auto x18 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x19 = V{1.5}*x18;
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x21 = V{1.5}*x20;
auto x22 = x21 + V{-1};
auto x23 = x19 + x22;
auto x24 = -x12*(x10 + V{-1})/x15 + V{1};
auto x25 = V{0.0277777777777778}*x16;
auto x26 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x27 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x28 = -x27;
auto x29 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x30 = x23 + x29;
auto x31 = -x26 + x30 - V{4.5}*x28*x28;
auto x32 = x24*x25;
auto x33 = V{0.111111111111111}*x16;
auto x34 = -x29;
auto x35 = V{1} - x19;
auto x36 = V{3}*x20 + x35;
auto x37 = x34 + x36;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x39 = V{4.5}*(x38*x38);
auto x40 = x26 + x30 - x39;
auto x41 = V{3}*x18;
auto x42 = x22 + x26 - x41;
auto x43 = x23 + x26 + x34 - V{4.5}*x27*x27;
auto x44 = x29 + x36;
auto x45 = -x21 + x26;
auto x46 = x29 + x35 + x39 + x45;
auto x47 = x41 + x45 + V{1};
auto x0 = -x11*(cell[0] + x17*x23 + V{0.444444444444444}) - x17*x23*x24 + V{-0.444444444444444};
auto x1 = -x11*(cell[1] + x25*x31 + V{0.0277777777777778}) - x31*x32 + V{-0.0277777777777778};
auto x2 = -x11*(cell[2] - x33*x37 + V{0.111111111111111}) + V{0.111111111111111}*x13*x15*x24*x37 + V{-0.111111111111111};
auto x3 = -x11*(cell[3] + x25*x40 + V{0.0277777777777778}) - x32*x40 + V{-0.0277777777777778};
auto x4 = -x11*(cell[4] + x33*x42 + V{0.111111111111111}) - x24*x33*x42 + V{-0.111111111111111};
auto x5 = -x11*(cell[5] + x25*x43 + V{0.0277777777777778}) - x32*x43 + V{-0.0277777777777778};
auto x6 = -x11*(cell[6] - x33*x44 + V{0.111111111111111}) + V{0.111111111111111}*x13*x15*x24*x44 + V{-0.111111111111111};
auto x7 = -x11*(cell[7] - x25*x46 + V{0.0277777777777778}) + V{0.0277777777777778}*x13*x15*x24*x46 + V{-0.0277777777777778};
auto x8 = -x11*(cell[8] - x33*x47 + V{0.111111111111111}) + V{0.111111111111111}*x13*x15*x24*x47 + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -x10 + x13*(V{1}*cell[0] + V{1}*cell[2] + V{1}*cell[6] + x14) + V{1}, x18 + x20 };
}
};

}

}

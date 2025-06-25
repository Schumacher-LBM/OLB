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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x22 = x21 + V{1};
auto x23 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x24 = V{1.5}*x23;
auto x25 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x26 = V{1.5}*x25;
auto x27 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x28 = V{1.5}*x27;
auto x29 = x26 + x28 + V{-1};
auto x30 = x24 + x29;
auto x31 = V{0.0555555555555556}*x19;
auto x32 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x33 = V{3}*x23;
auto x34 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x35 = V{3}*x25;
auto x36 = x24 + V{-1};
auto x37 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x38 = V{3}*x27;
auto x39 = V{0.0277777777777778}*x19;
auto x40 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x41 = V{4.5}*(x40*x40);
auto x42 = x30 + x32;
auto x43 = -x34;
auto x44 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x45 = -x44;
auto x46 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x47 = V{4.5}*(x46*x46);
auto x48 = -x37;
auto x49 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x50 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x49;
auto x51 = -x50;
auto x52 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x53 = V{4.5}*(x52*x52);
auto x54 = x30 + x34;
auto x55 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x49;
auto x56 = -x55;
auto x57 = -x26;
auto x58 = V{1} - x28;
auto x59 = x57 + x58;
auto x60 = x32 + x59;
auto x61 = -x24;
auto x62 = x34 + x61;
auto x63 = x37 + x61;
auto x64 = -x32;
auto x65 = x30 + x37;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x22*x30 + V{1});
auto x1 = -cell[1]*x20 - x31*(x22*(x29 + x32 - x33) + V{1});
auto x2 = -cell[2]*x20 - x31*(x22*(x28 + x34 - x35 + x36) + V{1});
auto x3 = -cell[3]*x20 - x31*(x22*(x26 + x36 + x37 - x38) + V{1});
auto x4 = -cell[4]*x20 - x39*(x22*(x34 - x41 + x42) + V{1});
auto x5 = -(cell[5]*x20 + x39*(x22*(x42 + x43 - V{4.5}*x45*x45) + V{1}));
auto x6 = -cell[6]*x20 - x39*(x22*(x37 + x42 - x47) + V{1});
auto x7 = -(cell[7]*x20 + x39*(x22*(x42 + x48 - V{4.5}*x51*x51) + V{1}));
auto x8 = -cell[8]*x20 - x39*(x22*(x37 - x53 + x54) + V{1});
auto x9 = -(cell[9]*x20 + x39*(x22*(x48 + x54 - V{4.5}*x56*x56) + V{1}));
auto x10 = -cell[10]*x20 + x31*(x22*(x33 + x60) + V{-1});
auto x11 = -cell[11]*x20 + x31*(x22*(x35 + x58 + x62) + V{-1});
auto x12 = -cell[12]*x20 + x31*(x22*(x38 + x57 + x63 + V{1}) + V{-1});
auto x13 = -cell[13]*x20 + x39*(x22*(x41 + x60 + x62) + V{-1});
auto x14 = -(cell[14]*x20 + x39*(x22*(x54 + x64 - V{4.5}*x44*x44) + V{1}));
auto x15 = -cell[15]*x20 + x39*(x22*(x47 + x60 + x63) + V{-1});
auto x16 = -(cell[16]*x20 + x39*(x22*(x64 + x65 - V{4.5}*x50*x50) + V{1}));
auto x17 = -cell[17]*x20 + x39*(x22*(x37 + x53 + x59 + x62) + V{-1});
auto x18 = -(cell[18]*x20 + x39*(x22*(x43 + x65 - V{4.5}*x55*x55) + V{1}));
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { x21 + V{1}, x23 + x25 + x27 };
}
};

}

}

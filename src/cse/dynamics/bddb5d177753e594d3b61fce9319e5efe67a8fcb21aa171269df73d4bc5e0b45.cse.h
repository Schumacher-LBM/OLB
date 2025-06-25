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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x27 = parameters.template get<descriptors::OMEGA>();
auto x28 = x27 + V{-1};
auto x29 = (cell[0] + V{2}*cell[11] + V{2}*cell[13] + cell[14] + cell[15] + V{2}*cell[16] + cell[17] + cell[18] + V{2}*cell[19] + cell[1] + V{2}*cell[21] + V{2}*cell[23] + V{2}*cell[25] + cell[2] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + V{1})/(cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x30 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x31 = V{1.5}*x30;
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x33 = V{1.5}*x32;
auto x34 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x35 = V{1.5}*x34;
auto x36 = x33 + x35 + V{-1};
auto x37 = x31 + x36;
auto x38 = V{0.0740740740740741}*x27;
auto x39 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x40 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x41 = x31 + V{-1};
auto x42 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x43 = V{0.0185185185185185}*x27;
auto x44 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x45 = x44*x44;
auto x46 = x37 + x39;
auto x47 = x40 + x46;
auto x48 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x49 = -x48;
auto x50 = -x40;
auto x51 = x46 + x50;
auto x52 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x53 = x52*x52;
auto x54 = -x42;
auto x55 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x56 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x55;
auto x57 = -x56;
auto x58 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x59 = x58*x58;
auto x60 = x37 + x40;
auto x61 = x42 + x60;
auto x62 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x55;
auto x63 = -x62;
auto x64 = V{0.00462962962962963}*x27;
auto x65 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + x44;
auto x66 = x65*x65;
auto x67 = x44 + x55;
auto x68 = -x67;
auto x69 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + x48;
auto x70 = -x69;
auto x71 = -x39;
auto x72 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x58;
auto x73 = -x33;
auto x74 = V{1} - x35;
auto x75 = x73 + x74;
auto x76 = -x31;
auto x77 = x40 + x76;
auto x78 = x42 + x75 + x77;
auto x79 = x39 + x75;
auto x80 = x42 + x76;
auto x81 = x77 + x79;
auto x82 = x79 + x80;
auto x83 = x37 + x42;
auto x84 = -x72;
auto x0 = -cell[0]*x28 - V{0.296296296296296}*x27*(x29*x37 + V{1});
auto x1 = -cell[1]*x28 - x38*(-x29*(V{3}*x30 - x36 - x39) + V{1});
auto x2 = -cell[2]*x28 - x38*(-x29*(V{3}*x32 - x35 - x40 - x41) + V{1});
auto x3 = -cell[3]*x28 - x38*(-x29*(-x33 + V{3}*x34 - x41 - x42) + V{1});
auto x4 = -cell[4]*x28 - x43*(-x29*(V{4.5}*x45 - x47) + V{1});
auto x5 = -(cell[5]*x28 + x43*(x29*(x51 - V{4.5}*x49*x49) + V{1}));
auto x6 = -cell[6]*x28 - x43*(-x29*(-x42 - x46 + V{4.5}*x53) + V{1});
auto x7 = -(cell[7]*x28 + x43*(x29*(x46 + x54 - V{4.5}*x57*x57) + V{1}));
auto x8 = -cell[8]*x28 - x43*(-x29*(V{4.5}*x59 - x61) + V{1});
auto x9 = -(cell[9]*x28 + x43*(x29*(x54 + x60 - V{4.5}*x63*x63) + V{1}));
auto x10 = -cell[10]*x28 - x64*(-x29*(-x42 - x47 + V{4.5}*x66) + V{1});
auto x11 = -(cell[11]*x28 + x64*(x29*(x47 + x54 - V{4.5}*x68*x68) + V{1}));
auto x12 = -(cell[12]*x28 + x64*(x29*(x42 + x51 - V{4.5}*x70*x70) + V{1}));
auto x13 = -(cell[13]*x28 + x64*(-x29*(x71 + x78 + V{4.5}*(x72*x72)) + V{1}));
auto x14 = -cell[14]*x28 - x38*(-x29*(V{3}*x30 + x79) + V{1});
auto x15 = -cell[15]*x28 - x38*(-x29*(V{3}*x32 + x74 + x77) + V{1});
auto x16 = -cell[16]*x28 - x38*(-x29*(V{3}*x34 + x73 + x80 + V{1}) + V{1});
auto x17 = -cell[17]*x28 - x43*(-x29*(V{4.5}*x45 + x81) + V{1});
auto x18 = -(cell[18]*x28 + x43*(x29*(x60 + x71 - V{4.5}*x48*x48) + V{1}));
auto x19 = -cell[19]*x28 - x43*(-x29*(V{4.5}*x53 + x82) + V{1});
auto x20 = -(cell[20]*x28 + x43*(x29*(x71 + x83 - V{4.5}*x56*x56) + V{1}));
auto x21 = -cell[21]*x28 - x43*(-x29*(V{4.5}*x59 + x78) + V{1});
auto x22 = -(cell[22]*x28 + x43*(x29*(x50 + x83 - V{4.5}*x62*x62) + V{1}));
auto x23 = -cell[23]*x28 - x64*(-x29*(x42 + V{4.5}*x66 + x81) + V{1});
auto x24 = -(cell[24]*x28 + x64*(-x29*(x54 + x81 + V{4.5}*(x67*x67)) + V{1}));
auto x25 = -(cell[25]*x28 + x64*(-x29*(x50 + x82 + V{4.5}*(x69*x69)) + V{1}));
auto x26 = -(cell[26]*x28 + x64*(x29*(x61 + x71 - V{4.5}*x84*x84) + V{1}));
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
cell[19] = x19;
cell[20] = x20;
cell[21] = x21;
cell[22] = x22;
cell[23] = x23;
cell[24] = x24;
cell[25] = x25;
cell[26] = x26;
return { V{1}*x29, x30 + x32 + x34 };
}
};

}

}

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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{0.5}*cell[13];
auto x22 = V{0.5}*cell[14];
auto x23 = V{0.5}*cell[15];
auto x24 = V{0.5}*cell[16];
auto x25 = V{0.5}*cell[17];
auto x26 = V{0.5}*cell[18];
auto x27 = V{0.5}*cell[4];
auto x28 = V{0.5}*cell[5];
auto x29 = V{0.5}*cell[6];
auto x30 = V{0.5}*cell[7];
auto x31 = V{0.5}*cell[8];
auto x32 = V{0.5}*cell[9];
auto x33 = V{0.5}*cell[0];
auto x34 = V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x30 + x31 + x32 + x33 + V{0.5};
auto x35 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x36 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x37 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x38 = V{1.5}*x35;
auto x39 = V{1.5}*x36;
auto x40 = V{1.5}*x37;
auto x41 = x39 + x40 + V{-1};
auto x42 = x38 + x41;
auto x43 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x44 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x45 = V{3}*x35;
auto x46 = V{0.166666666666667}*cell[10];
auto x47 = V{0.166666666666667}*cell[1];
auto x48 = V{0.166666666666667}*cell[17];
auto x49 = V{0.166666666666667}*cell[18];
auto x50 = V{0.166666666666667}*cell[8];
auto x51 = V{0.166666666666667}*cell[9];
auto x52 = V{0.166666666666667}*cell[11];
auto x53 = V{0.166666666666667}*cell[12];
auto x54 = V{0.166666666666667}*cell[13];
auto x55 = V{0.166666666666667}*cell[14];
auto x56 = V{0.166666666666667}*cell[15];
auto x57 = V{0.166666666666667}*cell[16];
auto x58 = V{0.166666666666667}*cell[2];
auto x59 = V{0.166666666666667}*cell[3];
auto x60 = V{0.166666666666667}*cell[4];
auto x61 = V{0.166666666666667}*cell[5];
auto x62 = V{0.166666666666667}*cell[6];
auto x63 = V{0.166666666666667}*cell[7];
auto x64 = V{0.166666666666667}*cell[0] + x46 + x47 + x48 + x49 + x50 + x51 + x52 + x53 + x54 + x55 + x56 + x57 + x58 + x59 + x60 + x61 + x62 + x63 + V{0.166666666666667};
auto x65 = V{0.0833333333333333}*cell[13];
auto x66 = V{0.0833333333333333}*cell[14];
auto x67 = V{0.0833333333333333}*cell[4];
auto x68 = V{0.0833333333333333}*cell[5];
auto x69 = -V{6.93889390390723e-18}*cell[0];
auto x70 = V{0.0833333333333333}*cell[0] + V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[11] + V{0.0833333333333333}*cell[12] + V{0.0833333333333333}*cell[13] + V{0.0833333333333333}*cell[14] + V{0.0833333333333333}*cell[15] + V{0.0833333333333333}*cell[16] + V{0.0833333333333333}*cell[17] + V{0.0833333333333333}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0833333333333333}*cell[2] + V{0.0833333333333333}*cell[3] + V{0.0833333333333333}*cell[4] + V{0.0833333333333333}*cell[5] + V{0.0833333333333333}*cell[6] + V{0.0833333333333333}*cell[7] + V{0.0833333333333333}*cell[8] + V{0.0833333333333333}*cell[9] + V{0.0833333333333333};
auto x71 = V{0.0833333333333333}*cell[12];
auto x72 = V{0.0833333333333333}*cell[3];
auto x73 = -x71 - x72;
auto x74 = x37*x70 + x65 + x66 + x67 + x68 + x69 + x73;
auto x75 = V{0.0833333333333333}*cell[15];
auto x76 = V{0.0833333333333333}*cell[16];
auto x77 = V{0.0833333333333333}*cell[6];
auto x78 = V{0.0833333333333333}*cell[7];
auto x79 = V{0.0833333333333333}*cell[11];
auto x80 = V{0.0833333333333333}*cell[2];
auto x81 = -x79 - x80;
auto x82 = x36*x70 + x75 + x76 + x77 + x78 + x81;
auto x83 = x20*(-x35*x64 + x46 + x47 - x48 - x49 - x50 - x51 + x74 + x82) + V{0.0555555555555556};
auto x84 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x85 = V{3}*x36;
auto x86 = x38 + V{-1};
auto x87 = V{0.0833333333333333}*cell[17];
auto x88 = V{0.0833333333333333}*cell[18];
auto x89 = V{0.0833333333333333}*cell[8];
auto x90 = V{0.0833333333333333}*cell[9];
auto x91 = V{0.0833333333333333}*cell[10];
auto x92 = V{0.0833333333333333}*cell[1];
auto x93 = -x91 - x92;
auto x94 = x35*x70 + x87 + x88 + x89 + x90 + x93;
auto x95 = x20*(-x36*x64 + x52 - x56 - x57 + x58 - x62 - x63 + x74 + x94) + V{0.0555555555555556};
auto x96 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x97 = V{3}*x37;
auto x98 = x20*(-x37*x64 + x53 - x54 - x55 + x59 - x60 - x61 + x69 + x82 + x94) + V{0.0555555555555556};
auto x99 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x100 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x101 = V{4.5}*(x100*x100);
auto x102 = x42 + x44;
auto x103 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x104 = x103 + V{1};
auto x105 = V{0.25}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x104;
auto x106 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x105;
auto x107 = V{0.0416666666666667}*cell[0] + V{0.0416666666666667}*cell[10] + V{0.0416666666666667}*cell[11] + V{0.0416666666666667}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0416666666666667}*cell[1] + V{0.0416666666666667}*cell[2] + V{0.0416666666666667}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + V{0.0416666666666667};
auto x108 = -V{0.0416666666666667}*cell[0];
auto x109 = V{0.0833333333333333}*cell[0] + x65 + x66 + x67 + x68 + x71 + x72 + x75 + x76 + x77 + x78 + x79 + x80 + x87 + x88 + x89 + x90 + x91 + x92 + V{0.0833333333333333};
auto x110 = V{0.0416666666666667}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{0.0416666666666667}*cell[1] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] + x108 - x109*x35;
auto x111 = V{0.0416666666666667}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{0.0416666666666667}*cell[2] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] - x109*x36;
auto x112 = x107*x37 + x110 + x111 + x73;
auto x113 = x20*(V{0.375}*cell[13] - V{0.125}*cell[14] + V{0.375}*cell[4] - V{0.125}*cell[5] - x106 + x112) + V{0.0277777777777778};
auto x114 = -x84;
auto x115 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x116 = -x115;
auto x117 = x20*(-V{0.125}*cell[13] + V{0.375}*cell[14] - V{0.125}*cell[4] + V{0.375}*cell[5] + x106 + x112) + V{0.0277777777777778};
auto x118 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x119 = V{4.5}*(x118*x118);
auto x120 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x105;
auto x121 = V{0.0416666666666667}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{0.0416666666666667}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] - x109*x37;
auto x122 = x107*x36 + x110 + x121 + x81;
auto x123 = x20*(V{0.375}*cell[15] - V{0.125}*cell[16] + V{0.375}*cell[6] - V{0.125}*cell[7] - x120 + x122) + V{0.0277777777777778};
auto x124 = -x96;
auto x125 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x126 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x125;
auto x127 = -x126;
auto x128 = x20*(-V{0.125}*cell[15] + V{0.375}*cell[16] - V{0.125}*cell[6] + V{0.375}*cell[7] + x120 + x122) + V{0.0277777777777778};
auto x129 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x130 = V{4.5}*(x129*x129);
auto x131 = x42 + x84;
auto x132 = V{0.25}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x104;
auto x133 = x107*x35 + x108 + x111 + x121 + x93;
auto x134 = x20*(V{0.375}*cell[17] - V{0.125}*cell[18] + V{0.375}*cell[8] - V{0.125}*cell[9] - x132 + x133) + V{0.0277777777777778};
auto x135 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x125;
auto x136 = -x135;
auto x137 = x20*(-V{0.125}*cell[17] + V{0.375}*cell[18] - V{0.125}*cell[8] + V{0.375}*cell[9] + x132 + x133) + V{0.0277777777777778};
auto x138 = x103 + V{1};
auto x139 = -x39;
auto x140 = V{1} - x40;
auto x141 = x139 + x140;
auto x142 = x141 + x44;
auto x143 = -x38;
auto x144 = x143 + x84;
auto x145 = x143 + x96;
auto x146 = -x44;
auto x147 = x42 + x96;
auto x0 = x20*(V{4.16333634234434e-17}*cell[10] + V{4.16333634234434e-17}*cell[11] + V{4.16333634234434e-17}*cell[12] + V{4.16333634234434e-17}*cell[1] + V{4.16333634234434e-17}*cell[2] + V{4.16333634234434e-17}*cell[3] + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x30 + x31 + x32 - x33 - x34*x35 - x34*x36 - x34*x37) - x42*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + V{-0.333333333333333};
auto x1 = -x43*(x41 + x44 - x45) - x83;
auto x2 = -x43*(x40 + x84 - x85 + x86) - x95;
auto x3 = -x43*(x39 + x86 + x96 - x97) - x98;
auto x4 = -x113 - x99*(-x101 + x102 + x84);
auto x5 = -(x117 + x99*(x102 + x114 - V{4.5}*x116*x116));
auto x6 = -x123 - x99*(x102 - x119 + x96);
auto x7 = -(x128 + x99*(x102 + x124 - V{4.5}*x127*x127));
auto x8 = -x134 - x99*(-x130 + x131 + x96);
auto x9 = -(x137 + x99*(x124 + x131 - V{4.5}*x136*x136));
auto x10 = V{0.0555555555555556}*x138*(x142 + x45) - x83;
auto x11 = V{0.0555555555555556}*x138*(x140 + x144 + x85) - x95;
auto x12 = V{0.0555555555555556}*x138*(x139 + x145 + x97 + V{1}) - x98;
auto x13 = -x113 + V{0.0277777777777778}*x138*(x101 + x142 + x144);
auto x14 = -(x117 + x99*(x131 + x146 - V{4.5}*x115*x115));
auto x15 = -x123 + V{0.0277777777777778}*x138*(x119 + x142 + x145);
auto x16 = -(x128 + x99*(x146 + x147 - V{4.5}*x126*x126));
auto x17 = -x134 + V{0.0277777777777778}*x138*(x130 + x141 + x144 + x96);
auto x18 = -(x137 + x99*(x114 + x147 - V{4.5}*x135*x135));
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
return { x104, x35 + x36 + x37 };
}
};

}

}

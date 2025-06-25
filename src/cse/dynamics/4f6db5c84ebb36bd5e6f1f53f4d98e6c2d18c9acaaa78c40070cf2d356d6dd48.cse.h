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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x22 = x21 + V{1};
auto x23 = x21 + V{1};
auto x24 = V{1} / (x23);
auto x25 = V{1}*cell[14];
auto x26 = V{1}*cell[16];
auto x27 = V{1}*cell[5];
auto x28 = V{1}*cell[7];
auto x29 = V{1}*cell[13] - V{1}*cell[4];
auto x30 = V{1}*cell[15] - V{1}*cell[6];
auto x31 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0) + x24*(V{1}*cell[10] - V{1}*cell[1] + x25 + x26 - x27 - x28 + x29 + x30);
auto x32 = x31*x31;
auto x33 = V{1.5}*x32;
auto x34 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x35 = V{1}*cell[18];
auto x36 = V{1}*cell[9];
auto x37 = V{1}*cell[17] - V{1}*cell[8];
auto x38 = x24*(V{1}*cell[11] - V{1}*cell[2] - x25 + x27 + x29 + x35 - x36 + x37);
auto x39 = x34 + x38;
auto x40 = x39*x39;
auto x41 = V{1.5}*x40;
auto x42 = V{1}*cell[12] - V{1}*cell[3] - x26 + x28 + x30 - x35 + x36 + x37;
auto x43 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2) + x24*x42;
auto x44 = x43*x43;
auto x45 = V{1.5}*x44;
auto x46 = x33 + x41 + x45 + V{-1};
auto x47 = V{0.5}*x19 + V{-1};
auto x48 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x49 = V{3}*cell[14];
auto x50 = V{3}*cell[16];
auto x51 = V{3}*cell[5];
auto x52 = V{3}*cell[7];
auto x53 = V{3}*cell[13] - V{3}*cell[4];
auto x54 = V{3}*cell[15] - V{3}*cell[6];
auto x55 = x24*(V{3}*cell[10] - V{3}*cell[1] + x49 + x50 - x51 - x52 + x53 + x54);
auto x56 = x48 + x55;
auto x57 = cell.template getFieldComponent<descriptors::FORCE>(0)*x56;
auto x58 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x59 = V{3}*cell[18];
auto x60 = V{3}*cell[9];
auto x61 = V{3}*cell[17] - V{3}*cell[8];
auto x62 = x24*(V{3}*cell[11] - V{3}*cell[2] - x49 + x51 + x53 + x59 - x60 + x61);
auto x63 = x58 + x62;
auto x64 = cell.template getFieldComponent<descriptors::FORCE>(1)*x63;
auto x65 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x66 = x24*(V{3}*cell[12] - V{3}*cell[3] - x50 + x52 + x54 - x59 + x60 + x61);
auto x67 = x65 + x66;
auto x68 = cell.template getFieldComponent<descriptors::FORCE>(2)*x67;
auto x69 = x64 + x68;
auto x70 = V{0.0555555555555556}*x19;
auto x71 = V{4.5}*cell[14];
auto x72 = V{4.5}*cell[16];
auto x73 = V{4.5}*cell[5];
auto x74 = V{4.5}*cell[7];
auto x75 = V{4.5}*cell[13] - V{4.5}*cell[4];
auto x76 = V{4.5}*cell[15] - V{4.5}*cell[6];
auto x77 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0) + x24*(V{4.5}*cell[10] - V{4.5}*cell[1] + x71 + x72 - x73 - x74 + x75 + x76);
auto x78 = x31*x77;
auto x79 = x46 + x56;
auto x80 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x81 = V{6}*cell[14];
auto x82 = V{6}*cell[16];
auto x83 = V{6}*cell[5];
auto x84 = V{6}*cell[7];
auto x85 = V{6}*cell[13] - V{6}*cell[4];
auto x86 = V{6}*cell[15] - V{6}*cell[6];
auto x87 = x24*(V{6}*cell[10] - V{6}*cell[1] + x81 + x82 - x83 - x84 + x85 + x86);
auto x88 = x80 + x87;
auto x89 = x88 + V{-3};
auto x90 = x23*x47;
auto x91 = V{0.0555555555555556}*x90;
auto x92 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x93 = V{4.5}*cell[18];
auto x94 = V{4.5}*cell[9];
auto x95 = V{4.5}*cell[17] - V{4.5}*cell[8];
auto x96 = x24*(V{4.5}*cell[11] - V{4.5}*cell[2] - x71 + x73 + x75 + x93 - x94 + x95);
auto x97 = x92 + x96;
auto x98 = x39*x97;
auto x99 = x46 + x63;
auto x100 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x101 = V{6}*cell[18];
auto x102 = V{6}*cell[9];
auto x103 = V{6}*cell[17] - V{6}*cell[8];
auto x104 = x24*(V{6}*cell[11] - V{6}*cell[2] + x101 - x102 + x103 - x81 + x83 + x85);
auto x105 = x100 + x104;
auto x106 = x105 + V{-3};
auto x107 = x57 + x68;
auto x108 = V{4.5}*cell[12] - V{4.5}*cell[3] - x72 + x74 + x76 - x93 + x94 + x95;
auto x109 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(2) + x108*x24;
auto x110 = x109*x43;
auto x111 = x46 + x67;
auto x112 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x113 = x24*(V{6}*cell[12] - V{6}*cell[3] - x101 + x102 + x103 - x82 + x84 + x86);
auto x114 = x112 + x113;
auto x115 = x114 + V{-3};
auto x116 = x57 + x64;
auto x117 = V{0.0277777777777778}*x19;
auto x118 = (x31 + x39)*(x77 + x97);
auto x119 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x120 = V{9}*cell[18];
auto x121 = V{9}*cell[5];
auto x122 = V{9}*cell[14];
auto x123 = V{9}*cell[9];
auto x124 = V{9}*cell[13] - V{9}*cell[4];
auto x125 = V{9}*cell[17] - V{9}*cell[8];
auto x126 = x24*(V{9}*cell[11] - V{9}*cell[2] + x120 + x121 - x122 - x123 + x124 + x125);
auto x127 = x119 + x126;
auto x128 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x129 = V{9}*cell[16];
auto x130 = V{9}*cell[7];
auto x131 = V{9}*cell[15] - V{9}*cell[6];
auto x132 = x24*(V{9}*cell[10] - V{9}*cell[1] - x121 + x122 + x124 + x129 - x130 + x131);
auto x133 = x128 + x132;
auto x134 = -x68;
auto x135 = V{0.0277777777777778}*x90;
auto x136 = x31 - x34 - x38;
auto x137 = x77 - x92 - x96;
auto x138 = -x58 - x62;
auto x139 = x127 + V{3};
auto x140 = -x80 - x87;
auto x141 = x105 + V{3};
auto x142 = -x128 - x132;
auto x143 = (x109 + x77)*(x31 + x43);
auto x144 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x145 = x24*(V{9}*cell[12] - V{9}*cell[3] - x120 + x123 + x125 - x129 + x130 + x131);
auto x146 = x144 + x145;
auto x147 = -x64;
auto x148 = -V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2) - x24*x42;
auto x149 = x148 + x31;
auto x150 = -V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(2) - x108*x24;
auto x151 = x150 + x77;
auto x152 = -x65 - x66;
auto x153 = x146 + V{3};
auto x154 = x114 + V{3};
auto x155 = (x109 + x97)*(x39 + x43);
auto x156 = -x57;
auto x157 = x148 + x39;
auto x158 = x150 + x97;
auto x159 = -x100 - x104;
auto x160 = -x119 - x126;
auto x161 = -x33 - x41 - x45 + V{1};
auto x162 = x161 + x56;
auto x163 = x88 + V{3};
auto x164 = x161 + x63;
auto x165 = -x48 - x55;
auto x166 = x133 + V{3};
auto x167 = -x112 - x113;
auto x168 = -x144 - x145;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x22*x46 + V{1}) + V{0.333333333333333}*x23*x47*(x57 + x69);
auto x1 = -cell[1]*x20 - x70*(x22*(-x78 + x79) + V{1}) - x91*(cell.template getFieldComponent<descriptors::FORCE>(0)*x89 - x69);
auto x2 = -cell[2]*x20 - x70*(x22*(-x98 + x99) + V{1}) - x91*(cell.template getFieldComponent<descriptors::FORCE>(1)*x106 - x107);
auto x3 = -cell[3]*x20 - x70*(x22*(-x110 + x111) + V{1}) - x91*(cell.template getFieldComponent<descriptors::FORCE>(2)*x115 - x116);
auto x4 = -cell[4]*x20 - x117*(x22*(-x118 + x63 + x79) + V{1}) - x135*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x127 + x89) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x106 + x133) + x134);
auto x5 = -cell[5]*x20 - x117*(x22*(-x136*x137 + x138 + x79) + V{1}) - x135*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x139 + x140) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x141 + x142) - x68);
auto x6 = -cell[6]*x20 - x117*(x22*(-x143 + x67 + x79) + V{1}) - x135*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x146 + x89) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x115 + x133) + x147);
auto x7 = -cell[7]*x20 - x117*(x22*(-x149*x151 + x152 + x79) + V{1}) - x135*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x140 + x153) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x142 + x154) - x64);
auto x8 = -cell[8]*x20 - x117*(x22*(-x155 + x67 + x99) + V{1}) - x135*(cell.template getFieldComponent<descriptors::FORCE>(1)*(x106 + x146) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x115 + x127) + x156);
auto x9 = -cell[9]*x20 - x117*(x22*(x152 - x157*x158 + x99) + V{1}) - x135*(-cell.template getFieldComponent<descriptors::FORCE>(1)*(x153 + x159) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x154 + x160) - x57);
auto x10 = -cell[10]*x20 + V{0.0555555555555556}*x19*(x22*(x162 + x78) + V{-1}) - x91*(cell.template getFieldComponent<descriptors::FORCE>(0)*x163 - x69);
auto x11 = -cell[11]*x20 + V{0.0555555555555556}*x19*(x22*(x164 + x98) + V{-1}) - x91*(cell.template getFieldComponent<descriptors::FORCE>(1)*x141 - x107);
auto x12 = -cell[12]*x20 + V{0.0555555555555556}*x19*(x22*(x110 + x161 + x67) + V{-1}) - x91*(cell.template getFieldComponent<descriptors::FORCE>(2)*x154 - x116);
auto x13 = -cell[13]*x20 - x135*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x127 + x163) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x133 + x141) + x134) + V{0.0277777777777778}*x19*(x22*(x118 + x162 + x63) + V{-1});
auto x14 = -cell[14]*x20 - x117*(x22*(-x136*x137 + x165 + x99) + V{1}) - x135*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x160 + x163) - cell.template getFieldComponent<descriptors::FORCE>(1)*(x159 + x166) - x68);
auto x15 = -cell[15]*x20 - x135*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x146 + x163) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x133 + x154) + x147) + V{0.0277777777777778}*x19*(x22*(x143 + x162 + x67) + V{-1});
auto x16 = -cell[16]*x20 - x117*(x22*(x111 - x149*x151 + x165) + V{1}) - x135*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x163 + x168) - cell.template getFieldComponent<descriptors::FORCE>(2)*(x166 + x167) - x64);
auto x17 = -cell[17]*x20 - x135*(cell.template getFieldComponent<descriptors::FORCE>(1)*(x141 + x146) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x127 + x154) + x156) + V{0.0277777777777778}*x19*(x22*(x155 + x164 + x67) + V{-1});
auto x18 = -cell[18]*x20 - x117*(x22*(x111 + x138 - x157*x158) + V{1}) - x135*(cell.template getFieldComponent<descriptors::FORCE>(1)*(x141 + x168) - cell.template getFieldComponent<descriptors::FORCE>(2)*(x139 + x167) - x57);
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
return { x23, x32 + x40 + x44 };
}
};

}

}

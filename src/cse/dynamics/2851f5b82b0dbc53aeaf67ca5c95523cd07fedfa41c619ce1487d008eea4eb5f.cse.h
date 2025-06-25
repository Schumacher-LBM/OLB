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
struct CSE<dynamics::Tuple<T, descriptors::D3Q15<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, multiphase::OmegaFromCell<collision::BGK>, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x17 = parameters.template get<multiphase::OMEGA_LIQUID>();
auto x18 = parameters.template get<multiphase::RHO_VAPOR>();
auto x19 = parameters.template get<multiphase::RHO_LIQUID>();
auto x16 = parameters.template get<multiphase::OMEGA_VAPOR>();
auto x15 = -x19;
auto x20 = x15 + x18;
auto x21 = V{1} / (x20);
auto x22 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x23 = x22 + V{1};
auto x24 = x16*(x15 + x23);
auto x25 = x21*x24;
auto x26 = x17*(-x18 + x23);
auto x27 = -x26/x20;
auto x28 = x25 + x27 + V{-1};
auto x29 = x24 - x26;
auto x30 = x21*x29;
auto x31 = x22 + V{1};
auto x32 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x33 = V{1} / (x23);
auto x34 = V{1}*cell[14];
auto x35 = V{1}*cell[7];
auto x36 = V{1}*cell[11];
auto x37 = V{1}*cell[12];
auto x38 = -V{1}*cell[4];
auto x39 = V{1}*cell[5];
auto x40 = x36 + x37 + x38 - x39;
auto x41 = V{1}*cell[13];
auto x42 = V{1}*cell[6];
auto x43 = x41 - x42;
auto x44 = x33*(-V{1}*cell[1] + V{1}*cell[8] + x34 - x35 + x40 + x43);
auto x45 = x32 + x44;
auto x46 = x45*x45;
auto x47 = V{1.5}*x46;
auto x48 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x49 = -x34 + x35;
auto x50 = x33*(-V{1}*cell[2] + V{1}*cell[9] + x40 - x41 + x42 + x49);
auto x51 = x48 + x50;
auto x52 = x51*x51;
auto x53 = V{1.5}*x52;
auto x54 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x55 = x33*(V{1}*cell[10] - V{1}*cell[3] + x36 - x37 + x38 + x39 + x43 + x49);
auto x56 = x54 + x55;
auto x57 = x56*x56;
auto x58 = V{1.5}*x57;
auto x59 = x47 + x53 + x58 + V{-1};
auto x60 = V{0.5}*x25 + V{0.5}*x27 + V{-1};
auto x61 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x62 = V{3}*cell[14];
auto x63 = V{3}*cell[7];
auto x64 = V{3}*cell[11];
auto x65 = V{3}*cell[12];
auto x66 = -V{3}*cell[4];
auto x67 = V{3}*cell[5];
auto x68 = x64 + x65 + x66 - x67;
auto x69 = V{3}*cell[13];
auto x70 = V{3}*cell[6];
auto x71 = x69 - x70;
auto x72 = x33*(-V{3}*cell[1] + V{3}*cell[8] + x62 - x63 + x68 + x71);
auto x73 = x61 + x72;
auto x74 = cell.template getFieldComponent<descriptors::FORCE>(0)*x73;
auto x75 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x76 = -x62 + x63;
auto x77 = x33*(-V{3}*cell[2] + V{3}*cell[9] + x68 - x69 + x70 + x76);
auto x78 = x75 + x77;
auto x79 = cell.template getFieldComponent<descriptors::FORCE>(1)*x78;
auto x80 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x81 = x33*(V{3}*cell[10] - V{3}*cell[3] + x64 - x65 + x66 + x67 + x71 + x76);
auto x82 = x80 + x81;
auto x83 = cell.template getFieldComponent<descriptors::FORCE>(2)*x82;
auto x84 = x79 + x83;
auto x85 = V{0.111111111111111}*x30;
auto x86 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x87 = V{4.5}*cell[14];
auto x88 = V{4.5}*cell[7];
auto x89 = V{4.5}*cell[11];
auto x90 = V{4.5}*cell[12];
auto x91 = -V{4.5}*cell[4];
auto x92 = V{4.5}*cell[5];
auto x93 = x89 + x90 + x91 - x92;
auto x94 = V{4.5}*cell[13];
auto x95 = V{4.5}*cell[6];
auto x96 = x94 - x95;
auto x97 = x33*(-V{4.5}*cell[1] + V{4.5}*cell[8] + x87 - x88 + x93 + x96);
auto x98 = x86 + x97;
auto x99 = x45*x98;
auto x100 = x59 + x73;
auto x101 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x102 = V{6}*cell[14];
auto x103 = V{6}*cell[7];
auto x104 = V{6}*cell[11];
auto x105 = V{6}*cell[12];
auto x106 = -V{6}*cell[4];
auto x107 = V{6}*cell[5];
auto x108 = x104 + x105 + x106 - x107;
auto x109 = V{6}*cell[13];
auto x110 = V{6}*cell[6];
auto x111 = x109 - x110;
auto x112 = x33*(-V{6}*cell[1] + V{6}*cell[8] + x102 - x103 + x108 + x111);
auto x113 = x101 + x112;
auto x114 = x113 + V{-3};
auto x115 = x23*x60;
auto x116 = V{0.111111111111111}*x115;
auto x117 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x118 = -x87 + x88;
auto x119 = x33*(-V{4.5}*cell[2] + V{4.5}*cell[9] + x118 + x93 - x94 + x95);
auto x120 = x117 + x119;
auto x121 = x120*x51;
auto x122 = x59 + x78;
auto x123 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x124 = -x102 + x103;
auto x125 = x33*(-V{6}*cell[2] + V{6}*cell[9] + x108 - x109 + x110 + x124);
auto x126 = x123 + x125;
auto x127 = x126 + V{-3};
auto x128 = x74 + x83;
auto x129 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x130 = x33*(V{4.5}*cell[10] - V{4.5}*cell[3] + x118 + x89 - x90 + x91 + x92 + x96);
auto x131 = x129 + x130;
auto x132 = x131*x56;
auto x133 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x134 = x33*(V{6}*cell[10] - V{6}*cell[3] + x104 - x105 + x106 + x107 + x111 + x124);
auto x135 = x133 + x134;
auto x136 = x135 + V{-3};
auto x137 = x74 + x79;
auto x138 = V{0.0138888888888889}*x30;
auto x139 = x45 + x51;
auto x140 = x120 + x98;
auto x141 = (x131 + x140)*(x139 + x56);
auto x142 = x100 + x78;
auto x143 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x144 = V{9}*cell[5];
auto x145 = V{9}*cell[12];
auto x146 = V{9}*cell[11];
auto x147 = V{9}*cell[7];
auto x148 = V{9}*cell[14];
auto x149 = -V{9}*cell[4];
auto x150 = x146 + x147 - x148 + x149;
auto x151 = V{9}*cell[13];
auto x152 = V{9}*cell[6];
auto x153 = x151 - x152;
auto x154 = x33*(V{9}*cell[10] - V{9}*cell[3] + x144 - x145 + x150 + x153);
auto x155 = x143 + x154;
auto x156 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x157 = -x144 + x145;
auto x158 = x33*(-V{9}*cell[2] + V{9}*cell[9] + x150 - x151 + x152 + x157);
auto x159 = x156 + x158;
auto x160 = x114 + x159;
auto x161 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x162 = x33*(-V{9}*cell[1] + V{9}*cell[8] + x146 - x147 + x148 + x149 + x153 + x157);
auto x163 = x161 + x162;
auto x164 = x127 + x155;
auto x165 = x159 + x163;
auto x166 = V{0.0138888888888889}*x115;
auto x167 = x139 - x54 - x55;
auto x168 = -x129 - x130 + x140;
auto x169 = -x80 - x81;
auto x170 = -x143 - x154;
auto x171 = x163 + x170;
auto x172 = -x133 - x134 + x165;
auto x173 = x45 - x48 - x50 + x56;
auto x174 = -x117 - x119 + x131 + x98;
auto x175 = -x75 - x77 + x82;
auto x176 = -x156 - x158;
auto x177 = x155 + x176;
auto x178 = x163 + x176;
auto x179 = x155 + x163;
auto x180 = -x123 - x125 + x179;
auto x181 = x126 + V{3};
auto x182 = -x161 - x162;
auto x183 = x135 + V{3};
auto x184 = x159 + x182;
auto x185 = x155 + x159;
auto x186 = -x101 - x112 + x185;
auto x187 = -x32 - x44 + x51 + x56;
auto x188 = x120 + x131 - x86 - x97;
auto x189 = -x47 - x53 - x58 + V{1};
auto x190 = x189 + x78;
auto x191 = -x61 - x72 + x82;
auto x192 = x113 + V{3};
auto x193 = x189 + x73;
auto x194 = x193 + x78;
auto x0 = -cell[0]*x28 + V{0.222222222222222}*x23*x60*(x74 + x84) - V{0.222222222222222}*x30*(x31*x59 + V{1});
auto x1 = -cell[1]*x28 - x116*(cell.template getFieldComponent<descriptors::FORCE>(0)*x114 - x84) - x85*(x31*(x100 - x99) + V{1});
auto x2 = -cell[2]*x28 - x116*(cell.template getFieldComponent<descriptors::FORCE>(1)*x127 - x128) - x85*(x31*(-x121 + x122) + V{1});
auto x3 = -cell[3]*x28 - x116*(cell.template getFieldComponent<descriptors::FORCE>(2)*x136 - x137) - x85*(x31*(-x132 + x59 + x82) + V{1});
auto x4 = -cell[4]*x28 - x138*(x31*(-x141 + x142 + x82) + V{1}) - x166*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x155 + x160) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x163 + x164) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x136 + x165));
auto x5 = -cell[5]*x28 - x138*(x31*(x142 - x167*x168 + x169) + V{1}) - x166*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x160 + x170) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x127 + x171) - cell.template getFieldComponent<descriptors::FORCE>(2)*(x172 + V{-3}));
auto x6 = -cell[6]*x28 - x138*(x31*(x100 - x173*x174 + x175) + V{1}) - x166*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x114 + x177) - cell.template getFieldComponent<descriptors::FORCE>(1)*(x180 + V{-3}) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x136 + x178));
auto x7 = -cell[7]*x28 - x166*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x186 + V{3}) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x155 + x181 + x182) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x183 + x184)) + V{0.0138888888888889}*x21*x29*(x31*(x187*x188 + x190 + x191) + V{-1});
auto x8 = -cell[8]*x28 - x116*(cell.template getFieldComponent<descriptors::FORCE>(0)*x192 - x84) + V{0.111111111111111}*x21*x29*(x31*(x193 + x99) + V{-1});
auto x9 = -cell[9]*x28 - x116*(cell.template getFieldComponent<descriptors::FORCE>(1)*x181 - x128) + V{0.111111111111111}*x21*x29*(x31*(x121 + x190) + V{-1});
auto x10 = -cell[10]*x28 - x116*(cell.template getFieldComponent<descriptors::FORCE>(2)*x183 - x137) + V{0.111111111111111}*x21*x29*(x31*(x132 + x189 + x82) + V{-1});
auto x11 = -cell[11]*x28 - x166*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x185 + x192) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x179 + x181) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x165 + x183)) + V{0.0138888888888889}*x21*x29*(x31*(x141 + x194 + x82) + V{-1});
auto x12 = -cell[12]*x28 - x166*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x159 + x170 + x192) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x171 + x181) - cell.template getFieldComponent<descriptors::FORCE>(2)*(x172 + V{3})) + V{0.0138888888888889}*x21*x29*(x31*(x167*x168 + x169 + x194) + V{-1});
auto x13 = -cell[13]*x28 - x166*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x177 + x192) - cell.template getFieldComponent<descriptors::FORCE>(1)*(x180 + V{3}) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x178 + x183)) + V{0.0138888888888889}*x21*x29*(x31*(x173*x174 + x175 + x193) + V{-1});
auto x14 = -cell[14]*x28 - x138*(x31*(x122 - x187*x188 + x191) + V{1}) - x166*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x186 + V{-3}) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x164 + x182) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x136 + x184));
cell.template getFieldPointer<descriptors::OMEGA>()[0] = x30;
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
return { x23, x46 + x52 + x57 };
}
};

}

}

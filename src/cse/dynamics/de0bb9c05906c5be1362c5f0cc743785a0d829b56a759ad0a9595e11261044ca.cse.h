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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1};
auto x22 = -1/x21;
auto x23 = V{2}*cell[14] + V{2}*cell[4] + V{2}*cell[8] + V{2}*cell[9];
auto x24 = cell[0] + cell[10] + cell[12] + cell[15] + cell[16] + cell[1] + V{2}*cell[2] + cell[3] + cell[6] + cell[7] + x23 + V{1};
auto x25 = x22*x24;
auto x26 = V{0.0277777777777778}*x25;
auto x27 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x28 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x29 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x30 = V{1.5}*x29;
auto x31 = -x30;
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x33 = V{1.5}*x32;
auto x34 = V{1} - x33;
auto x35 = x31 + x34;
auto x36 = x28 + x35;
auto x37 = V{3}*x27 + x36;
auto x38 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x39 = V{1.5}*x27;
auto x40 = -x39;
auto x41 = x38 + x40;
auto x42 = V{3}*x29 + x34 + x41;
auto x43 = x30 + x33 + V{-1};
auto x44 = -V{3}*x27 + x28 + x43;
auto x45 = -x44;
auto x46 = x39 + V{-1};
auto x47 = -V{3}*x29 + x33 + x38 + x46;
auto x48 = -x47;
auto x49 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x50 = x49*x49;
auto x51 = x36 + x41 + V{4.5}*x50;
auto x52 = -x28;
auto x53 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x54 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x53;
auto x55 = -V{4.5}*x54*x54;
auto x56 = x39 + x43;
auto x57 = x38 + x56;
auto x58 = x52 + x55 + x57;
auto x59 = -x58;
auto x60 = -x54;
auto x61 = -x38 + x56;
auto x62 = x28 + x61;
auto x63 = x62 - V{4.5}*x60*x60;
auto x64 = -x63;
auto x65 = x28 - V{4.5}*x50 + x57;
auto x66 = -x65;
auto x67 = V{0.0555555555555556}*x25;
auto x68 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x69 = x30 - V{3}*x32 + x46 + x68;
auto x70 = -x69;
auto x71 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x72 = x56 + x68;
auto x73 = x52 + x72 - V{4.5}*x71*x71;
auto x74 = -x73;
auto x75 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x53;
auto x76 = -x75;
auto x77 = x61 + x68;
auto x78 = x77 - V{4.5}*x76*x76;
auto x79 = -x78;
auto x80 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x81 = x80*x80;
auto x82 = x28 + x72 - V{4.5}*x81;
auto x83 = -x82;
auto x84 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x85 = x84*x84;
auto x86 = x57 + x68 - V{4.5}*x85;
auto x87 = -x86;
auto x88 = -x56;
auto x89 = V{0.333333333333333}*x25*x88;
auto x90 = V{0.0555555555555555}*x25;
auto x91 = V{0.0555555555555555}*x25;
auto x92 = V{0.0277777777777778}*x25;
auto x93 = V{1.66533453693773e-16}*cell[10] + V{9.71445146547012e-17}*cell[12] + V{6.66133814775094e-16}*cell[14] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{1.66533453693773e-16}*cell[1] + V{4.44089209850063e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{2.22044604925031e-16}*cell[4] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{4.44089209850063e-16}*cell[9] + V{2.22044604925031e-16};
auto x94 = V{0.0555555555555555}*x25*x45 + V{0.0277777777777778}*x25*x64 + V{0.0277777777777778}*x25*x66 + V{0.111111111111111}*x25*x70 + x37*x67 + x42*x67 + x48*x90 + x51*x92 + x59*x92 + x74*x91 + x79*x90 + x83*x91 + x87*x90 + x89 + x93;
auto x95 = V{0.166666666666667}*cell[14];
auto x96 = V{0.166666666666667}*cell[4];
auto x97 = V{0.0833333333333333}*cell[12];
auto x98 = V{0.0833333333333333}*cell[3];
auto x99 = V{0.333333333333333}*cell[8];
auto x100 = V{0.333333333333333}*cell[9];
auto x101 = V{1} / (x21);
auto x102 = x101*x24;
auto x103 = V{0.00462962962962963}*x102;
auto x104 = V{0.00925925925925926}*x102;
auto x105 = V{0.0462962962962963}*x102;
auto x106 = V{0.00925925925925926}*x102;
auto x107 = x103*x42;
auto x108 = V{0.00462962962962963}*x102;
auto x109 = x108*x82;
auto x110 = x108*x73;
auto x111 = V{0.166666666666667}*cell[2];
auto x112 = V{0.00231481481481482}*x102;
auto x113 = x112*x65;
auto x114 = x112*x63;
auto x115 = x112*x58;
auto x116 = V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] + x104*x69 - x111 + x112*x51 - x113 - x114 - x115 + V{-0.0555555555555555};
auto x117 = V{0.0555555555555556}*x19;
auto x118 = -x95;
auto x119 = -x96;
auto x120 = -V{0.00462962962962963}*x101*x24*x47 + x107 + x109 + x110 + x118 + x119 + x97 + x98 + V{0.0555555555555555};
auto x121 = V{0.0833333333333333}*cell[10];
auto x122 = V{0.0833333333333333}*cell[1];
auto x123 = V{0.166666666666667}*cell[8];
auto x124 = -x123;
auto x125 = V{0.166666666666667}*cell[9];
auto x126 = -x125;
auto x127 = x103*x37;
auto x128 = x108*x86;
auto x129 = x108*x78;
auto x130 = -V{0.00462962962962963}*x101*x24*x44 + x121 + x122 + x124 + x126 + x127 + x128 + x129;
auto x131 = V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.333333333333333}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7];
auto x132 = V{0.333333333333333}*cell[14];
auto x133 = V{0.333333333333333}*cell[4];
auto x134 = V{0.00462962962962963}*x102;
auto x135 = V{0.00231481481481481}*x102;
auto x136 = V{0.00462962962962963}*x102;
auto x137 = V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] - x112*x78 - x112*x86 + x136*x37 - x136*x44 + V{0.0138888888888889};
auto x138 = V{0.00115740740740741}*x102;
auto x139 = V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.166666666666667}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x106*x69 + x138*x51 - x138*x58 - x138*x63 - x138*x65;
auto x140 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x135*x42 + x135*x47 + x137 + x139;
auto x141 = -V{0.166666666666667}*cell[14] + V{0.833333333333333}*cell[4] + x134*x73 + x140;
auto x142 = V{0.0277777777777778}*x19;
auto x143 = V{0.0277777777777778}*x102;
auto x144 = -x68;
auto x145 = -x71;
auto x146 = x144 + x28 + x56 - V{4.5}*x145*x145;
auto x147 = V{0.0231481481481481}*x102;
auto x148 = V{0.833333333333333}*cell[14] - V{0.166666666666667}*cell[4] + x134*x82 + x140;
auto x149 = V{0.0115740740740741}*x102;
auto x150 = V{0.0162037037037037}*x102;
auto x151 = V{0.00231481481481481}*x102;
auto x152 = V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] - x112*x73 - x112*x82 + x136*x42 - x136*x47;
auto x153 = -V{0.0833333333333333}*cell[2] + x103*x69 + x137 + x152;
auto x154 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x151*x58 + x151*x63 + x153;
auto x155 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x151*x51 + x151*x65 + x153;
auto x156 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x135*x37 + x135*x44 + x139 + x152 + V{0.0138888888888889};
auto x157 = V{0.833333333333333}*cell[8] - V{0.166666666666667}*cell[9] + x134*x78 + x156;
auto x158 = -V{0.166666666666667}*cell[8] + V{0.833333333333333}*cell[9] + x134*x86 + x156;
auto x159 = -V{0.0833333333333334}*cell[15] - V{0.0833333333333334}*cell[16] - V{0.0833333333333334}*cell[6] - V{0.0833333333333334}*cell[7] - V{0.00231481481481482}*x101*x24*x51 - V{0.00925925925925926}*x101*x24*x69 + x111 + x113 + x114 + x115;
auto x160 = V{0.00462962962962963}*x25;
auto x161 = V{0.00462962962962963}*x25;
auto x162 = x40 + x68;
auto x163 = x162 + x31 + V{3}*x32 + V{1};
auto x164 = x162 + x36 + V{4.5}*x81;
auto x165 = x35 + x41 + x68 + V{4.5}*x85;
auto x166 = -V{4.5}*x75*x75;
auto x167 = x144 + x166 + x57;
auto x168 = V{0.0555555555555555}*x102;
auto x169 = V{0.0555555555555555}*x102;
auto x170 = V{0.0277777777777778}*x102;
auto x171 = V{0.0555555555555556}*x102;
auto x0 = -V{0.333333333333333}*x19*(-x22*x88*x94 + V{1}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[12] + V{1}*cell[15] + V{1}*cell[16] + V{0.5}*cell[1] + V{1}*cell[2] + V{0.5}*cell[3] + V{1}*cell[6] + V{1}*cell[7] + x23 - x26*x37 - x26*x42 - x26*x45 - x26*x48 - x26*x51 - x26*x59 - x26*x64 - x26*x66 - x67*x70 - x67*x74 - x67*x79 - x67*x83 - x67*x87 - x89 + V{0.833333333333333});
auto x1 = -x117*(-x22*x45*x94 + V{1}) - x20*(V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] - x100 + x103*x47 + x104*x78 + x104*x86 + x105*x44 + x106*x37 - x107 - x109 - x110 + x116 + x95 + x96 - x97 - x98 - x99);
auto x2 = -x117*(-x22*x70*x94 + V{1}) - x20*(V{0.00462962962962963}*x101*x24*x58 + V{0.00462962962962963}*x101*x24*x63 + V{0.00462962962962963}*x101*x24*x65 + V{0.037037037037037}*x101*x24*x69 - x103*x51 - x120 - x130 - x131);
auto x3 = -x117*(-x22*x48*x94 + V{1}) - x20*(V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + x103*x44 + x104*x73 + x104*x82 + x105*x47 + x106*x42 + x116 - x121 - x122 + x123 + x125 - x127 - x128 - x129 - x132 - x133);
auto x4 = -x142*(-x22*x83*x94 + V{1}) - x20*(x103*x82 + x141);
auto x5 = -x142*(x146*x22*x94 + V{1}) - x20*(x143*x146 - x147*x73 + x148);
auto x6 = -x142*(-x22*x66*x94 + V{1}) - x20*(x149*x51 + x150*x65 + x154);
auto x7 = -x142*(-x22*x64*x94 + V{1}) - x20*(-x149*x58 + x150*x63 + x155);
auto x8 = -x142*(-x22*x87*x94 + V{1}) - x20*(x103*x86 + x157);
auto x9 = -x142*(-x22*x79*x94 + V{1}) - x20*(x103*x78 + x158);
auto x10 = -x117*(-x22*x37*x94 + V{1}) - x20*(V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] - x100 + V{0.00925925925925926}*x101*x24*x78 + V{0.00925925925925926}*x101*x24*x86 - x105*x37 - x106*x44 - x120 - x159 - x99);
auto x11 = -x117*(-x163*x22*x94 + V{1}) + x20*(x118 + x119 + x121 + x122 + x124 + x126 + x131 + x160*x74 + x160*x79 + x160*x83 + x160*x87 - x161*x37 - x161*x42 - x161*x45 - x161*x48 - x161*x51 - x161*x59 - x161*x64 - x161*x66 - x163*x67 + V{0.0185185185185185}*x25*x70 + x97 + x98 + V{0.0555555555555555});
auto x12 = -x117*(-x22*x42*x94 + V{1}) - x20*(V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + V{0.00925925925925926}*x101*x24*x73 + V{0.00925925925925926}*x101*x24*x82 - x105*x42 - x106*x47 - x130 - x132 - x133 - x159 + V{-0.0555555555555555});
auto x13 = -x142*(-x164*x22*x94 + V{1}) - x20*(x141 - x143*x164 - x147*x82);
auto x14 = -x142*(-x22*x74*x94 + V{1}) - x20*(x103*x73 + x148);
auto x15 = -x142*(-x22*x51*x94 + V{1}) - x20*(-x149*x65 - x150*x51 + x154);
auto x16 = -x142*(-x22*x59*x94 + V{1}) - x20*(-x149*x63 + x150*x58 + x155);
auto x17 = -x142*(-x165*x22*x94 + V{1}) - x20*(-x143*x165 - x147*x86 + x157);
auto x18 = -x142*(x167*x22*x94 + V{1}) - x20*(x143*x167 - x147*x78 + x158);
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
return { -V{1}*x101*(V{0.0555555555555555}*x102*x44 + V{0.333333333333333}*x102*x56 + V{0.0277777777777778}*x102*x65 + V{0.111111111111111}*x102*x69 + V{0.0277777777777778}*x102*(x55 + x62) + x168*x47 + x168*x86 + x168*(x166 + x77) + x169*x73 + x169*x82 - x170*x51 + x170*x58 - x171*x37 - x171*x42 + x93), x27 + x29 + x32 };
}
};

}

}

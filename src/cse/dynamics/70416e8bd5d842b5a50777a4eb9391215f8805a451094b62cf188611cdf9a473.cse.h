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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<2, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1};
auto x22 = -1/x21;
auto x23 = V{2}*cell[16] + V{2}*cell[18] + V{2}*cell[6] + V{2}*cell[8];
auto x24 = cell[0] + cell[10] + cell[11] + cell[13] + cell[14] + cell[1] + cell[2] + V{2}*cell[3] + cell[4] + cell[5] + x23 + V{1};
auto x25 = x22*x24;
auto x26 = V{0.0277777777777778}*x25;
auto x27 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x28 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x29 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x30 = V{1.5}*x29;
auto x31 = -x30;
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x33 = V{1.5}*x32;
auto x34 = V{1} - x33;
auto x35 = x31 + x34;
auto x36 = x28 + x35;
auto x37 = V{3}*x27 + x36;
auto x38 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
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
auto x49 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x50 = x49*x49;
auto x51 = x36 + x41 + V{4.5}*x50;
auto x52 = -x28;
auto x53 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x54 = -V{4.5}*x53*x53;
auto x55 = x39 + x43;
auto x56 = x38 + x55;
auto x57 = x52 + x54 + x56;
auto x58 = -x57;
auto x59 = -x53;
auto x60 = -x38 + x55;
auto x61 = x28 + x60;
auto x62 = x61 - V{4.5}*x59*x59;
auto x63 = -x62;
auto x64 = x28 - V{4.5}*x50 + x56;
auto x65 = -x64;
auto x66 = V{0.0555555555555556}*x25;
auto x67 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x68 = x30 - V{3}*x32 + x46 + x67;
auto x69 = -x68;
auto x70 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x71 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x70;
auto x72 = x55 + x67;
auto x73 = x52 + x72 - V{4.5}*x71*x71;
auto x74 = -x73;
auto x75 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x70;
auto x76 = x60 + x67 - V{4.5}*x75*x75;
auto x77 = -x76;
auto x78 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x79 = x78*x78;
auto x80 = x28 + x72 - V{4.5}*x79;
auto x81 = -x80;
auto x82 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x83 = x82*x82;
auto x84 = x56 + x67 - V{4.5}*x83;
auto x85 = -x84;
auto x86 = -x55;
auto x87 = V{0.333333333333333}*x25*x86;
auto x88 = V{0.0555555555555555}*x25;
auto x89 = V{0.0555555555555555}*x25;
auto x90 = V{0.0277777777777778}*x25;
auto x91 = V{1.66533453693773e-16}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{6.66133814775094e-16}*cell[16] + V{4.44089209850063e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{1.11022302462516e-16}*cell[2] + V{4.44089209850063e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{6.66133814775094e-16}*cell[6] + V{2.22044604925031e-16};
auto x92 = V{0.0555555555555555}*x25*x45 + V{0.0277777777777778}*x25*x63 + V{0.0277777777777778}*x25*x65 + V{0.111111111111111}*x25*x69 + x37*x66 + x42*x66 + x48*x88 + x51*x90 + x58*x90 + x74*x89 + x77*x88 + x81*x89 + x85*x88 + x87 + x91;
auto x93 = V{0.166666666666667}*cell[16];
auto x94 = V{0.166666666666667}*cell[6];
auto x95 = V{0.0833333333333333}*cell[11];
auto x96 = V{0.0833333333333333}*cell[2];
auto x97 = V{0.333333333333333}*cell[18];
auto x98 = V{0.333333333333333}*cell[8];
auto x99 = V{1} / (x21);
auto x100 = x24*x99;
auto x101 = V{0.00462962962962963}*x100;
auto x102 = V{0.00925925925925926}*x100;
auto x103 = V{0.0462962962962963}*x100;
auto x104 = V{0.00925925925925926}*x100;
auto x105 = x101*x42;
auto x106 = V{0.00462962962962963}*x100;
auto x107 = x106*x80;
auto x108 = x106*x73;
auto x109 = V{0.166666666666667}*cell[3];
auto x110 = V{0.00231481481481482}*x100;
auto x111 = x110*x64;
auto x112 = x110*x62;
auto x113 = x110*x57;
auto x114 = V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x102*x68 - x109 + x110*x51 - x111 - x112 - x113 + V{-0.0555555555555555};
auto x115 = V{0.0555555555555556}*x19;
auto x116 = V{0.166666666666667}*cell[18];
auto x117 = V{0.166666666666667}*cell[8];
auto x118 = V{0.0833333333333333}*cell[10];
auto x119 = V{0.0833333333333333}*cell[1];
auto x120 = V{0.333333333333333}*cell[16];
auto x121 = V{0.333333333333333}*cell[6];
auto x122 = x101*x37;
auto x123 = x106*x84;
auto x124 = x106*x76;
auto x125 = -x93;
auto x126 = -x94;
auto x127 = x105 + x107 + x108 + x125 + x126 - V{0.00462962962962963}*x24*x47*x99 + x95 + x96 + V{0.0555555555555555};
auto x128 = -x116;
auto x129 = -x117;
auto x130 = x118 + x119 + x122 + x123 + x124 + x128 + x129 - V{0.00462962962962963}*x24*x44*x99;
auto x131 = V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.333333333333333}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5];
auto x132 = V{0.0115740740740741}*x100;
auto x133 = V{0.0162037037037037}*x100;
auto x134 = V{0.00231481481481481}*x100;
auto x135 = V{0.00462962962962963}*x100;
auto x136 = V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] - x110*x76 - x110*x84 + x135*x37 - x135*x44 + V{0.0138888888888889};
auto x137 = V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] - x110*x73 - x110*x80 + x135*x42 - x135*x47;
auto x138 = -V{0.0833333333333333}*cell[3] + x101*x68 + x136 + x137;
auto x139 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x134*x57 + x134*x62 + x138;
auto x140 = V{0.0277777777777778}*x19;
auto x141 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x134*x51 + x134*x64 + x138;
auto x142 = V{0.00462962962962963}*x100;
auto x143 = V{0.00231481481481481}*x100;
auto x144 = V{0.00115740740740741}*x100;
auto x145 = V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.166666666666667}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x104*x68 + x144*x51 - x144*x57 - x144*x62 - x144*x64;
auto x146 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x136 - x143*x42 + x143*x47 + x145;
auto x147 = -V{0.166666666666667}*cell[16] + V{0.833333333333333}*cell[6] + x142*x73 + x146;
auto x148 = V{0.0277777777777778}*x100;
auto x149 = -x67;
auto x150 = -x71;
auto x151 = x149 + x28 + x55 - V{4.5}*x150*x150;
auto x152 = V{0.0231481481481481}*x100;
auto x153 = V{0.833333333333333}*cell[16] - V{0.166666666666667}*cell[6] + x142*x80 + x146;
auto x154 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x137 - x143*x37 + x143*x44 + x145 + V{0.0138888888888889};
auto x155 = -V{0.166666666666667}*cell[18] + V{0.833333333333333}*cell[8] + x142*x76 + x154;
auto x156 = -x75;
auto x157 = x149 + x56 - V{4.5}*x156*x156;
auto x158 = V{0.833333333333333}*cell[18] - V{0.166666666666667}*cell[8] + x142*x84 + x154;
auto x159 = -V{0.0833333333333334}*cell[13] - V{0.0833333333333334}*cell[14] - V{0.0833333333333334}*cell[4] - V{0.0833333333333334}*cell[5] + x109 + x111 + x112 + x113 - V{0.00231481481481482}*x24*x51*x99 - V{0.00925925925925926}*x24*x68*x99;
auto x160 = V{0.00462962962962963}*x25;
auto x161 = V{0.00462962962962963}*x25;
auto x162 = x40 + x67;
auto x163 = x162 + x31 + V{3}*x32 + V{1};
auto x164 = x162 + x36 + V{4.5}*x79;
auto x165 = x35 + x41 + x67 + V{4.5}*x83;
auto x166 = V{0.0555555555555555}*x100;
auto x167 = V{0.0555555555555555}*x100;
auto x168 = V{0.0277777777777778}*x100;
auto x169 = V{0.0555555555555556}*x100;
auto x0 = -V{0.333333333333333}*x19*(-x22*x86*x92 + V{1}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{1}*cell[13] + V{1}*cell[14] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{1}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + x23 - x26*x37 - x26*x42 - x26*x45 - x26*x48 - x26*x51 - x26*x58 - x26*x63 - x26*x65 - x66*x69 - x66*x74 - x66*x77 - x66*x81 - x66*x85 - x87 + V{0.833333333333333});
auto x1 = -x115*(-x22*x45*x92 + V{1}) - x20*(V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] + x101*x47 + x102*x76 + x102*x84 + x103*x44 + x104*x37 - x105 - x107 - x108 + x114 + x93 + x94 - x95 - x96 - x97 - x98);
auto x2 = -x115*(-x22*x48*x92 + V{1}) - x20*(V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + x101*x44 + x102*x73 + x102*x80 + x103*x47 + x104*x42 + x114 + x116 + x117 - x118 - x119 - x120 - x121 - x122 - x123 - x124);
auto x3 = -x115*(-x22*x69*x92 + V{1}) - x20*(-x101*x51 - x127 - x130 - x131 + V{0.00462962962962963}*x24*x57*x99 + V{0.00462962962962963}*x24*x62*x99 + V{0.00462962962962963}*x24*x64*x99 + V{0.037037037037037}*x24*x68*x99);
auto x4 = -x140*(-x22*x65*x92 + V{1}) - x20*(x132*x51 + x133*x64 + x139);
auto x5 = -x140*(-x22*x63*x92 + V{1}) - x20*(-x132*x57 + x133*x62 + x141);
auto x6 = -x140*(-x22*x81*x92 + V{1}) - x20*(x101*x80 + x147);
auto x7 = -x140*(x151*x22*x92 + V{1}) - x20*(x148*x151 - x152*x73 + x153);
auto x8 = -x140*(-x22*x85*x92 + V{1}) - x20*(x101*x84 + x155);
auto x9 = -x140*(x157*x22*x92 + V{1}) - x20*(x148*x157 - x152*x76 + x158);
auto x10 = -x115*(-x22*x37*x92 + V{1}) - x20*(V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] - x103*x37 - x104*x44 - x127 - x159 + V{0.00925925925925926}*x24*x76*x99 + V{0.00925925925925926}*x24*x84*x99 - x97 - x98);
auto x11 = -x115*(-x22*x42*x92 + V{1}) - x20*(V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] - x103*x42 - x104*x47 - x120 - x121 - x130 - x159 + V{0.00925925925925926}*x24*x73*x99 + V{0.00925925925925926}*x24*x80*x99 + V{-0.0555555555555555});
auto x12 = -x115*(-x163*x22*x92 + V{1}) + x20*(x118 + x119 + x125 + x126 + x128 + x129 + x131 + x160*x74 + x160*x77 + x160*x81 + x160*x85 - x161*x37 - x161*x42 - x161*x45 - x161*x48 - x161*x51 - x161*x58 - x161*x63 - x161*x65 - x163*x66 + V{0.0185185185185185}*x25*x69 + x95 + x96 + V{0.0555555555555555});
auto x13 = -x140*(-x22*x51*x92 + V{1}) - x20*(-x132*x64 - x133*x51 + x139);
auto x14 = -x140*(-x22*x58*x92 + V{1}) - x20*(-x132*x62 + x133*x57 + x141);
auto x15 = -x140*(-x164*x22*x92 + V{1}) - x20*(x147 - x148*x164 - x152*x80);
auto x16 = -x140*(-x22*x74*x92 + V{1}) - x20*(x101*x73 + x153);
auto x17 = -x140*(-x165*x22*x92 + V{1}) - x20*(-x148*x165 - x152*x84 + x155);
auto x18 = -x140*(-x22*x77*x92 + V{1}) - x20*(x101*x76 + x158);
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
return { -V{1}*x99*(V{0.0555555555555555}*x100*x44 + V{0.333333333333333}*x100*x55 + V{0.0277777777777778}*x100*x64 + V{0.111111111111111}*x100*x68 + V{0.0277777777777778}*x100*(x54 + x61) + x166*x47 + x166*x76 + x166*x84 + x167*x73 + x167*x80 - x168*x51 + x168*x57 - x169*x37 - x169*x42 + x91), x27 + x29 + x32 };
}
};

}

}

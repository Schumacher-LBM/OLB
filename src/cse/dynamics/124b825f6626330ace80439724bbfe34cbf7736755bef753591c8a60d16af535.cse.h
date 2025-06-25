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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1};
auto x22 = -1/x21;
auto x23 = V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7];
auto x24 = cell[0] + cell[11] + cell[12] + cell[17] + cell[18] + V{2}*cell[1] + cell[2] + cell[3] + cell[8] + cell[9] + x23 + V{1};
auto x25 = x22*x24;
auto x26 = V{0.0277777777777778}*x25;
auto x27 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x28 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x29 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x30 = V{1.5}*x29;
auto x31 = -x30;
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
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
auto x49 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x50 = x49*x49;
auto x51 = x36 + x41 + V{4.5}*x50;
auto x52 = -x28;
auto x53 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x54 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x53;
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
auto x68 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x69 = x30 - V{3}*x32 + x46 + x68;
auto x70 = -x69;
auto x71 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x72 = -x71;
auto x73 = x56 + x68;
auto x74 = x52 + x73;
auto x75 = x74 - V{4.5}*x72*x72;
auto x76 = -x75;
auto x77 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x53;
auto x78 = -x77;
auto x79 = x61 + x68;
auto x80 = x79 - V{4.5}*x78*x78;
auto x81 = -x80;
auto x82 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x83 = x82*x82;
auto x84 = x28 + x73 - V{4.5}*x83;
auto x85 = -x84;
auto x86 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x87 = x86*x86;
auto x88 = x57 + x68 - V{4.5}*x87;
auto x89 = -x88;
auto x90 = -x56;
auto x91 = V{0.333333333333333}*x25*x90;
auto x92 = V{0.0555555555555555}*x25;
auto x93 = V{0.0555555555555555}*x25;
auto x94 = V{0.0277777777777778}*x25;
auto x95 = V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{4.44089209850063e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{6.66133814775094e-16}*cell[4] + V{6.66133814775094e-16}*cell[5] + V{8.88178419700125e-16}*cell[6] + V{4.44089209850063e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + V{2.22044604925031e-16};
auto x96 = V{0.0555555555555555}*x25*x45 + V{0.0277777777777778}*x25*x64 + V{0.0277777777777778}*x25*x66 + V{0.111111111111111}*x25*x70 + x37*x67 + x42*x67 + x48*x92 + x51*x94 + x59*x94 + x76*x93 + x81*x92 + x85*x93 + x89*x92 + x91 + x95;
auto x97 = V{1} / (x21);
auto x98 = x24*x97;
auto x99 = V{0.00462962962962963}*x98;
auto x100 = V{0.0833333333333333}*cell[12];
auto x101 = V{0.0833333333333333}*cell[3];
auto x102 = V{0.166666666666667}*cell[4];
auto x103 = -x102;
auto x104 = V{0.166666666666667}*cell[5];
auto x105 = -x104;
auto x106 = x42*x99;
auto x107 = V{0.00462962962962963}*x98;
auto x108 = x107*x84;
auto x109 = x107*x75;
auto x110 = x100 + x101 + x103 + x105 + x106 + x108 + x109 - V{0.00462962962962963}*x24*x47*x97 + V{0.0555555555555555};
auto x111 = V{0.0833333333333333}*cell[11];
auto x112 = V{0.0833333333333333}*cell[2];
auto x113 = V{0.166666666666667}*cell[6];
auto x114 = -x113;
auto x115 = V{0.166666666666667}*cell[7];
auto x116 = -x115;
auto x117 = x37*x99;
auto x118 = x107*x88;
auto x119 = x107*x80;
auto x120 = x111 + x112 + x114 + x116 + x117 + x118 + x119 - V{0.00462962962962963}*x24*x44*x97;
auto x121 = V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.333333333333333}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9];
auto x122 = V{0.0555555555555556}*x19;
auto x123 = V{0.333333333333333}*cell[6];
auto x124 = V{0.333333333333333}*cell[7];
auto x125 = V{0.00925925925925926}*x98;
auto x126 = V{0.0462962962962963}*x98;
auto x127 = V{0.00925925925925926}*x98;
auto x128 = V{0.166666666666667}*cell[1];
auto x129 = V{0.00231481481481482}*x98;
auto x130 = x129*x65;
auto x131 = x129*x63;
auto x132 = x129*x58;
auto x133 = V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x125*x69 - x128 + x129*x51 - x130 - x131 - x132 + V{-0.0555555555555555};
auto x134 = V{0.333333333333333}*cell[4];
auto x135 = V{0.333333333333333}*cell[5];
auto x136 = V{0.00462962962962963}*x98;
auto x137 = V{0.00231481481481481}*x98;
auto x138 = V{0.00462962962962963}*x98;
auto x139 = V{0.0833333333333333}*cell[11] + V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] - x129*x80 - x129*x88 + x138*x37 - x138*x44 + V{0.0138888888888889};
auto x140 = V{0.00115740740740741}*x98;
auto x141 = V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.166666666666667}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] - x127*x69 + x140*x51 - x140*x58 - x140*x63 - x140*x65;
auto x142 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x137*x42 + x137*x47 + x139 + x141;
auto x143 = V{0.833333333333333}*cell[4] - V{0.166666666666667}*cell[5] + x136*x75 + x142;
auto x144 = V{0.0277777777777778}*x19;
auto x145 = -V{0.166666666666667}*cell[4] + V{0.833333333333333}*cell[5] + x136*x84 + x142;
auto x146 = V{0.0833333333333333}*cell[12] + V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] - x129*x75 - x129*x84 + x138*x42 - x138*x47;
auto x147 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x137*x37 + x137*x44 + x141 + x146 + V{0.0138888888888889};
auto x148 = V{0.833333333333333}*cell[6] - V{0.166666666666667}*cell[7] + x136*x80 + x147;
auto x149 = -V{0.166666666666667}*cell[6] + V{0.833333333333333}*cell[7] + x136*x88 + x147;
auto x150 = V{0.0115740740740741}*x98;
auto x151 = V{0.0162037037037037}*x98;
auto x152 = V{0.00231481481481481}*x98;
auto x153 = -V{0.0833333333333333}*cell[1] + x139 + x146 + x69*x99;
auto x154 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] + x152*x58 + x152*x63 + x153;
auto x155 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] - x152*x51 + x152*x65 + x153;
auto x156 = V{0.00462962962962963}*x25;
auto x157 = V{0.00462962962962963}*x25;
auto x158 = x40 + x68;
auto x159 = x158 + x31 + V{3}*x32 + V{1};
auto x160 = -V{0.0833333333333334}*cell[17] - V{0.0833333333333334}*cell[18] - V{0.0833333333333334}*cell[8] - V{0.0833333333333334}*cell[9] + x128 + x130 + x131 + x132 - V{0.00231481481481482}*x24*x51*x97 - V{0.00925925925925926}*x24*x69*x97;
auto x161 = V{0.0277777777777778}*x98;
auto x162 = x158 + x36 + V{4.5}*x83;
auto x163 = V{0.0231481481481481}*x98;
auto x164 = -x68;
auto x165 = -V{4.5}*x71*x71;
auto x166 = x164 + x165 + x28 + x56;
auto x167 = x35 + x41 + x68 + V{4.5}*x87;
auto x168 = -V{4.5}*x77*x77;
auto x169 = x164 + x168 + x57;
auto x170 = V{0.0555555555555555}*x98;
auto x171 = V{0.0555555555555555}*x98;
auto x172 = V{0.0277777777777778}*x98;
auto x173 = V{0.0555555555555556}*x98;
auto x0 = -V{0.333333333333333}*x19*(-x22*x90*x96 + V{1}) + x20*(V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[17] + V{1}*cell[18] + V{1}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[8] + V{1}*cell[9] + x23 - x26*x37 - x26*x42 - x26*x45 - x26*x48 - x26*x51 - x26*x59 - x26*x64 - x26*x66 - x67*x70 - x67*x76 - x67*x81 - x67*x85 - x67*x89 - x91 + V{0.833333333333333});
auto x1 = -x122*(-x22*x70*x96 + V{1}) - x20*(-x110 - x120 - x121 + V{0.00462962962962963}*x24*x58*x97 + V{0.00462962962962963}*x24*x63*x97 + V{0.00462962962962963}*x24*x65*x97 + V{0.037037037037037}*x24*x69*x97 - x51*x99);
auto x2 = -x122*(-x22*x45*x96 + V{1}) - x20*(V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] - x100 - x101 + x102 + x104 - x106 - x108 - x109 - x123 - x124 + x125*x80 + x125*x88 + x126*x44 + x127*x37 + x133 + x47*x99);
auto x3 = -x122*(-x22*x48*x96 + V{1}) - x20*(V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] - x111 - x112 + x113 + x115 - x117 - x118 - x119 + x125*x75 + x125*x84 + x126*x47 + x127*x42 + x133 - x134 - x135 + x44*x99);
auto x4 = -x144*(-x22*x85*x96 + V{1}) - x20*(x143 + x84*x99);
auto x5 = -x144*(-x22*x76*x96 + V{1}) - x20*(x145 + x75*x99);
auto x6 = -x144*(-x22*x89*x96 + V{1}) - x20*(x148 + x88*x99);
auto x7 = -x144*(-x22*x81*x96 + V{1}) - x20*(x149 + x80*x99);
auto x8 = -x144*(-x22*x66*x96 + V{1}) - x20*(x150*x51 + x151*x65 + x154);
auto x9 = -x144*(-x22*x64*x96 + V{1}) - x20*(-x150*x58 + x151*x63 + x155);
auto x10 = -x122*(-x159*x22*x96 + V{1}) + x20*(x100 + x101 + x103 + x105 + x111 + x112 + x114 + x116 + x121 + x156*x76 + x156*x81 + x156*x85 + x156*x89 - x157*x37 - x157*x42 - x157*x45 - x157*x48 - x157*x51 - x157*x59 - x157*x64 - x157*x66 - x159*x67 + V{0.0185185185185185}*x25*x70 + V{0.0555555555555555});
auto x11 = -x122*(-x22*x37*x96 + V{1}) - x20*(V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] - x110 - x123 - x124 - x126*x37 - x127*x44 - x160 + V{0.00925925925925926}*x24*x80*x97 + V{0.00925925925925926}*x24*x88*x97);
auto x12 = -x122*(-x22*x42*x96 + V{1}) - x20*(V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] - x120 - x126*x42 - x127*x47 - x134 - x135 - x160 + V{0.00925925925925926}*x24*x75*x97 + V{0.00925925925925926}*x24*x84*x97 + V{-0.0555555555555555});
auto x13 = -x144*(-x162*x22*x96 + V{1}) - x20*(x143 - x161*x162 - x163*x84);
auto x14 = -x144*(x166*x22*x96 + V{1}) - x20*(x145 + x161*x166 - x163*x75);
auto x15 = -x144*(-x167*x22*x96 + V{1}) - x20*(x148 - x161*x167 - x163*x88);
auto x16 = -x144*(x169*x22*x96 + V{1}) - x20*(x149 + x161*x169 - x163*x80);
auto x17 = -x144*(-x22*x51*x96 + V{1}) - x20*(-x150*x65 - x151*x51 + x154);
auto x18 = -x144*(-x22*x59*x96 + V{1}) - x20*(-x150*x63 + x151*x58 + x155);
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
return { -V{1}*x97*(x170*x47 + x170*x88 + x170*(x168 + x79) + x171*x84 + x171*(x165 + x74) - x172*x51 + x172*x58 - x173*x37 - x173*x42 + V{0.0555555555555555}*x44*x98 + V{0.333333333333333}*x56*x98 + V{0.0277777777777778}*x65*x98 + V{0.111111111111111}*x69*x98 + x95 + V{0.0277777777777778}*x98*(x55 + x62)), x27 + x29 + x32 };
}
};

}

}

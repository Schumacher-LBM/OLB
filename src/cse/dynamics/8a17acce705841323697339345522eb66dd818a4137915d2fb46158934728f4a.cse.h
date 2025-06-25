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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x22 = V{2}*cell[13] + V{2}*cell[17] + V{2}*cell[18] + V{2}*cell[5];
auto x23 = cell[0] + cell[10] + V{2}*cell[11] + cell[12] + cell[15] + cell[16] + cell[1] + cell[3] + cell[6] + cell[7] + x22 + V{1};
auto x24 = x21*x23;
auto x25 = V{0.0277777777777778}*x24;
auto x26 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x27 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x28 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x29 = V{1.5}*x28;
auto x30 = -x29;
auto x31 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x32 = V{1.5}*x31;
auto x33 = V{1} - x32;
auto x34 = x30 + x33;
auto x35 = x27 + x34;
auto x36 = V{3}*x26 + x35;
auto x37 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x38 = V{1.5}*x26;
auto x39 = -x38;
auto x40 = x37 + x39;
auto x41 = V{3}*x28 + x33 + x40;
auto x42 = x29 + x32 + V{-1};
auto x43 = -V{3}*x26 + x27 + x42;
auto x44 = -x43;
auto x45 = x38 + V{-1};
auto x46 = -V{3}*x28 + x32 + x37 + x45;
auto x47 = -x46;
auto x48 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x49 = x48*x48;
auto x50 = x35 + x40 + V{4.5}*x49;
auto x51 = -x27;
auto x52 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x53 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x52;
auto x54 = -V{4.5}*x53*x53;
auto x55 = x38 + x42;
auto x56 = x37 + x55;
auto x57 = x51 + x54 + x56;
auto x58 = -x57;
auto x59 = -x53;
auto x60 = -x37;
auto x61 = x27 + x55;
auto x62 = x60 + x61;
auto x63 = x62 - V{4.5}*x59*x59;
auto x64 = -x63;
auto x65 = x27 - V{4.5}*x49 + x56;
auto x66 = -x65;
auto x67 = V{0.0555555555555556}*x24;
auto x68 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x69 = x39 + x68;
auto x70 = x30 + V{3}*x31 + x69 + V{1};
auto x71 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x72 = V{4.5}*(x71*x71);
auto x73 = x35 + x69 + x72;
auto x74 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x75 = V{4.5}*(x74*x74);
auto x76 = x34 + x40 + x68 + x75;
auto x77 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x78 = -x77;
auto x79 = -x68;
auto x80 = x61 + x79;
auto x81 = x80 - V{4.5}*x78*x78;
auto x82 = -x81;
auto x83 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x52;
auto x84 = x56 + x79 - V{4.5}*x83*x83;
auto x85 = -x84;
auto x86 = V{0.333333333333333}*x24;
auto x87 = -x55;
auto x88 = x86*x87;
auto x89 = V{0.0277777777777778}*x24;
auto x90 = V{0.0555555555555555}*x24;
auto x91 = V{0.0555555555555555}*x24;
auto x92 = V{0.0277777777777778}*x24;
auto x93 = V{0.0555555555555555}*x24;
auto x94 = V{0.0277777777777778}*x24;
auto x95 = V{1.66533453693773e-16}*cell[10] + V{4.44089209850063e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{6.66133814775094e-16}*cell[13] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{8.88178419700125e-16}*cell[17] + V{4.44089209850063e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{1.11022302462516e-16}*cell[3] + V{6.66133814775094e-16}*cell[5] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{0.111111111111111}*x24*x70 + x36*x67 + x41*x67 + x50*x94 + x73*x93 + x76*x90 + V{2.22044604925031e-16};
auto x96 = x44*x91 + x47*x90 + x58*x94 + x64*x89 + x66*x92 + x82*x93 + x85*x90 + x88 + x95;
auto x97 = V{0.0462962962962963}*x24;
auto x98 = V{0.00925925925925926}*x24;
auto x99 = V{0.166666666666667}*cell[13];
auto x100 = V{0.166666666666667}*cell[5];
auto x101 = V{0.0833333333333333}*cell[12];
auto x102 = V{0.0833333333333333}*cell[3];
auto x103 = V{0.00462962962962963}*x24;
auto x104 = x103*x41;
auto x105 = V{0.00925925925925926}*x24;
auto x106 = V{0.00462962962962963}*x24;
auto x107 = x103*x46;
auto x108 = x106*x73;
auto x109 = V{0.00231481481481482}*x24;
auto x110 = -V{0.166666666666667}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] + x105*x70 - x109*x50 + x109*x57 + x109*x63 + x109*x65 + V{-0.0555555555555555};
auto x111 = V{0.166666666666667}*cell[10] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.166666666666667}*cell[1] + x100 - x101 - x102 + x104 + x105*x76 - x105*x84 + x106*x81 - x107 - x108 + x110 + x99;
auto x112 = V{0.0555555555555556}*x19;
auto x113 = -x29 + V{3}*x31 - x45 - x68;
auto x114 = V{0.0833333333333333}*cell[10];
auto x115 = V{0.0833333333333333}*cell[1];
auto x116 = V{0.166666666666667}*cell[17];
auto x117 = V{0.166666666666667}*cell[18];
auto x118 = x106*x76;
auto x119 = x103*x36;
auto x120 = -V{0.333333333333333}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] - x100 + x101 + x102 - x103*x50 - x104 + x108 + x114 + x115 - x116 - x117 + x118 - x119 - x99 + V{0.0555555555555555};
auto x121 = x103*x43;
auto x122 = V{0.166666666666667}*cell[12] - V{0.333333333333333}*cell[13] + V{0.166666666666667}*cell[3] - V{0.333333333333333}*cell[5] + x105*x73 - x105*x81 + x106*x84 + x110 - x114 - x115 + x116 + x117 - x118 + x119 - x121;
auto x123 = x61 + x68 - x72;
auto x124 = V{0.0231481481481481}*x24;
auto x125 = V{0.00462962962962963}*x24;
auto x126 = V{0.00231481481481481}*x24;
auto x127 = V{0.00462962962962963}*x24;
auto x128 = V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333333}*cell[1] - x109*x76 + x109*x84 - x127*x36 + x127*x43 + V{0.0138888888888889};
auto x129 = V{0.00115740740740741}*x24;
auto x130 = V{0.166666666666667}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x129*x50 + x129*x57 + x129*x63 + x129*x65 - x70*x98;
auto x131 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x126*x41 - x126*x46 + x128 + x130;
auto x132 = V{0.833333333333333}*cell[13] - V{0.166666666666667}*cell[5] - x125*x81 + x131;
auto x133 = V{0.0277777777777778}*x19;
auto x134 = -V{0.166666666666667}*cell[13] + V{0.833333333333333}*cell[5] + x125*x73 + x131;
auto x135 = V{0.0115740740740741}*x24;
auto x136 = V{0.0162037037037037}*x24;
auto x137 = V{0.00231481481481481}*x24;
auto x138 = V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[5] - x109*x73 + x109*x81 - x127*x41 + x127*x46;
auto x139 = -V{0.0833333333333333}*cell[11] + x103*x70 + x128 + x138;
auto x140 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x137*x57 - x137*x63 + x139;
auto x141 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x137*x50 - x137*x65 + x139;
auto x142 = x56 + x68 - x75;
auto x143 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x126*x36 - x126*x43 + x130 + x138 + V{0.0138888888888889};
auto x144 = V{0.833333333333333}*cell[17] - V{0.166666666666667}*cell[18] - x125*x84 + x143;
auto x145 = -x83;
auto x146 = x55 + x68;
auto x147 = x146 + x60 - V{4.5}*x145*x145;
auto x148 = -V{0.166666666666667}*cell[17] + V{0.833333333333333}*cell[18] + x125*x76 + x143;
auto x149 = -V{4.5}*x77*x77;
auto x150 = x146 + x149 + x51;
auto x0 = -V{0.333333333333333}*x19*(-x21*x87*x96 + V{1}) + x20*(V{0.5}*cell[10] + V{1}*cell[11] + V{0.5}*cell[12] + V{1}*cell[15] + V{1}*cell[16] + V{0.5}*cell[1] + V{0.5}*cell[3] + V{1}*cell[6] + V{1}*cell[7] + x22 - x25*x36 - x25*x41 - x25*x44 - x25*x47 - x25*x50 - x25*x58 - x25*x64 - x25*x66 - x67*x70 - x67*x73 - x67*x76 - x67*x82 - x67*x85 - x88 + V{0.833333333333333});
auto x1 = -x112*(-x21*x44*x96 + V{1}) - x20*(x111 - x36*x98 - x43*x97);
auto x2 = -x112*(-x113*x21*x96 + V{1}) + x20*(-x103*x44 - x103*x47 - x103*x58 - x103*x64 - x103*x66 + x106*x82 + x106*x85 - x113*x67 + x120 + V{0.0185185185185185}*x24*x70);
auto x3 = -x112*(-x21*x47*x96 + V{1}) - x20*(x122 - x41*x98 - x46*x97);
auto x4 = -x133*(x123*x21*x96 + V{1}) - x20*(-x123*x25 - x124*x73 + x132);
auto x5 = -x133*(-x21*x82*x96 + V{1}) - x20*(-x103*x81 + x134);
auto x6 = -x133*(-x21*x66*x96 + V{1}) - x20*(-x135*x50 - x136*x65 + x140);
auto x7 = -x133*(-x21*x64*x96 + V{1}) - x20*(x135*x57 - x136*x63 + x141);
auto x8 = -x133*(x142*x21*x96 + V{1}) - x20*(-x124*x76 - x142*x25 + x144);
auto x9 = -x133*(x147*x21*x96 + V{1}) - x20*(x124*x84 - x147*x25 + x148);
auto x10 = -x112*(-x21*x36*x96 + V{1}) - x20*(x111 + x36*x97 + x43*x98);
auto x11 = -x112*(-x21*x70*x96 + V{1}) - x20*(-x103*x57 - x103*x63 - x103*x65 - x107 - x120 - x121 + V{0.037037037037037}*x21*x23*x70 + V{0.00462962962962963}*x21*x23*x81 + V{0.00462962962962963}*x21*x23*x84);
auto x12 = -x112*(-x21*x41*x96 + V{1}) - x20*(x122 + x41*x97 + x46*x98);
auto x13 = -x133*(-x21*x73*x96 + V{1}) - x20*(x103*x73 + x132);
auto x14 = -x133*(x150*x21*x96 + V{1}) - x20*(x124*x81 + x134 - x150*x25);
auto x15 = -x133*(-x21*x50*x96 + V{1}) - x20*(x135*x65 + x136*x50 + x140);
auto x16 = -x133*(-x21*x58*x96 + V{1}) - x20*(x135*x63 - x136*x57 + x141);
auto x17 = -x133*(-x21*x76*x96 + V{1}) - x20*(x103*x76 + x144);
auto x18 = -x133*(-x21*x85*x96 + V{1}) - x20*(-x103*x84 + x148);
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
return { V{1}*x21*(-x43*x91 - x46*x90 - x55*x86 - x57*x94 - x65*x92 - x84*x90 - x89*(x54 + x62) - x93*(x149 + x80) + x95), x26 + x28 + x31 };
}
};

}

}

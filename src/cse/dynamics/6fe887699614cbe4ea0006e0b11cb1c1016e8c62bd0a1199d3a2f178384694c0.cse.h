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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<2, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x22 = V{2}*cell[15] + V{2}*cell[17] + V{2}*cell[7] + V{2}*cell[9];
auto x23 = cell[0] + cell[10] + cell[11] + V{2}*cell[12] + cell[13] + cell[14] + cell[1] + cell[2] + cell[4] + cell[5] + x22 + V{1};
auto x24 = x21*x23;
auto x25 = V{0.0277777777777778}*x24;
auto x26 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x27 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x28 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x29 = V{1.5}*x28;
auto x30 = -x29;
auto x31 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x32 = V{1.5}*x31;
auto x33 = V{1} - x32;
auto x34 = x30 + x33;
auto x35 = x27 + x34;
auto x36 = V{3}*x26 + x35;
auto x37 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
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
auto x48 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x49 = x48*x48;
auto x50 = x35 + x40 + V{4.5}*x49;
auto x51 = -x27;
auto x52 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x53 = -V{4.5}*x52*x52;
auto x54 = x38 + x42;
auto x55 = x37 + x54;
auto x56 = x51 + x53 + x55;
auto x57 = -x56;
auto x58 = -x52;
auto x59 = -x37;
auto x60 = x27 + x54;
auto x61 = x59 + x60;
auto x62 = x61 - V{4.5}*x58*x58;
auto x63 = -x62;
auto x64 = x27 - V{4.5}*x49 + x55;
auto x65 = -x64;
auto x66 = V{0.0555555555555556}*x24;
auto x67 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x68 = x39 + x67;
auto x69 = x30 + V{3}*x31 + x68 + V{1};
auto x70 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x71 = V{4.5}*(x70*x70);
auto x72 = x35 + x68 + x71;
auto x73 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x74 = V{4.5}*(x73*x73);
auto x75 = x34 + x40 + x67 + x74;
auto x76 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x77 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x76;
auto x78 = -x77;
auto x79 = -x67;
auto x80 = x60 + x79;
auto x81 = x80 - V{4.5}*x78*x78;
auto x82 = -x81;
auto x83 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x76;
auto x84 = -x83;
auto x85 = x55 + x79;
auto x86 = x85 - V{4.5}*x84*x84;
auto x87 = -x86;
auto x88 = V{0.333333333333333}*x24;
auto x89 = -x54;
auto x90 = x88*x89;
auto x91 = V{0.0277777777777778}*x24;
auto x92 = V{0.0555555555555555}*x24;
auto x93 = V{0.0555555555555555}*x24;
auto x94 = V{0.0277777777777778}*x24;
auto x95 = V{0.0555555555555555}*x24;
auto x96 = V{0.0277777777777778}*x24;
auto x97 = V{1.66533453693773e-16}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{4.44089209850063e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{6.66133814775094e-16}*cell[15] + V{4.44089209850063e-16}*cell[17] + V{1.66533453693773e-16}*cell[1] + V{1.11022302462516e-16}*cell[2] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{6.66133814775094e-16}*cell[7] + V{0.111111111111111}*x24*x69 + x36*x66 + x41*x66 + x50*x96 + x72*x95 + x75*x92 + V{2.22044604925031e-16};
auto x98 = x44*x93 + x47*x92 + x57*x96 + x63*x91 + x65*x94 + x82*x95 + x87*x92 + x90 + x97;
auto x99 = V{0.0462962962962963}*x24;
auto x100 = V{0.00925925925925926}*x24;
auto x101 = V{0.166666666666667}*cell[15];
auto x102 = V{0.166666666666667}*cell[7];
auto x103 = V{0.0833333333333333}*cell[11];
auto x104 = V{0.0833333333333333}*cell[2];
auto x105 = V{0.00462962962962963}*x24;
auto x106 = x105*x41;
auto x107 = V{0.00925925925925926}*x24;
auto x108 = V{0.00462962962962963}*x24;
auto x109 = x105*x46;
auto x110 = x108*x72;
auto x111 = V{0.00231481481481482}*x24;
auto x112 = -V{0.166666666666667}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x107*x69 - x111*x50 + x111*x56 + x111*x62 + x111*x64 + V{-0.0555555555555555};
auto x113 = V{0.166666666666667}*cell[10] - V{0.333333333333333}*cell[17] + V{0.166666666666667}*cell[1] - V{0.333333333333333}*cell[9] + x101 + x102 - x103 - x104 + x106 + x107*x75 - x107*x86 + x108*x81 - x109 - x110 + x112;
auto x114 = V{0.0555555555555556}*x19;
auto x115 = V{0.166666666666667}*cell[17];
auto x116 = V{0.166666666666667}*cell[9];
auto x117 = V{0.0833333333333333}*cell[10];
auto x118 = V{0.0833333333333333}*cell[1];
auto x119 = x105*x36;
auto x120 = x105*x43;
auto x121 = x108*x75;
auto x122 = V{0.166666666666667}*cell[11] - V{0.333333333333333}*cell[15] + V{0.166666666666667}*cell[2] - V{0.333333333333333}*cell[7] + x107*x72 - x107*x81 + x108*x86 + x112 + x115 + x116 - x117 - x118 + x119 - x120 - x121;
auto x123 = -x29 + V{3}*x31 - x45 - x67;
auto x124 = -V{0.333333333333333}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] - x101 - x102 + x103 + x104 - x105*x50 - x106 + x110 - x115 - x116 + x117 + x118 - x119 + x121 + V{0.0555555555555555};
auto x125 = V{0.0115740740740741}*x24;
auto x126 = V{0.0162037037037037}*x24;
auto x127 = V{0.00231481481481481}*x24;
auto x128 = V{0.00462962962962963}*x24;
auto x129 = V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[9] - x111*x75 + x111*x86 - x128*x36 + x128*x43 + V{0.0138888888888889};
auto x130 = V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[7] - x111*x72 + x111*x81 - x128*x41 + x128*x46;
auto x131 = -V{0.0833333333333333}*cell[12] + x105*x69 + x129 + x130;
auto x132 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x127*x56 - x127*x62 + x131;
auto x133 = V{0.0277777777777778}*x19;
auto x134 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x127*x50 - x127*x64 + x131;
auto x135 = x60 + x67 - x71;
auto x136 = V{0.0231481481481481}*x24;
auto x137 = V{0.00462962962962963}*x24;
auto x138 = V{0.00231481481481481}*x24;
auto x139 = V{0.00115740740740741}*x24;
auto x140 = V{0.166666666666667}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x100*x69 - x139*x50 + x139*x56 + x139*x62 + x139*x64;
auto x141 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x129 + x138*x41 - x138*x46 + x140;
auto x142 = V{0.833333333333333}*cell[15] - V{0.166666666666667}*cell[7] - x137*x81 + x141;
auto x143 = -V{0.166666666666667}*cell[15] + V{0.833333333333333}*cell[7] + x137*x72 + x141;
auto x144 = x55 + x67 - x74;
auto x145 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x130 + x138*x36 - x138*x43 + x140 + V{0.0138888888888889};
auto x146 = V{0.833333333333333}*cell[17] - V{0.166666666666667}*cell[9] - x137*x86 + x145;
auto x147 = -V{0.166666666666667}*cell[17] + V{0.833333333333333}*cell[9] + x137*x75 + x145;
auto x148 = -V{4.5}*x77*x77;
auto x149 = x54 + x67;
auto x150 = x148 + x149 + x51;
auto x151 = -V{4.5}*x83*x83;
auto x152 = x149 + x151 + x59;
auto x0 = -V{0.333333333333333}*x19*(-x21*x89*x98 + V{1}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{1}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{1}*cell[4] + V{1}*cell[5] + x22 - x25*x36 - x25*x41 - x25*x44 - x25*x47 - x25*x50 - x25*x57 - x25*x63 - x25*x65 - x66*x69 - x66*x72 - x66*x75 - x66*x82 - x66*x87 - x90 + V{0.833333333333333});
auto x1 = -x114*(-x21*x44*x98 + V{1}) - x20*(-x100*x36 + x113 - x43*x99);
auto x2 = -x114*(-x21*x47*x98 + V{1}) - x20*(-x100*x41 + x122 - x46*x99);
auto x3 = -x114*(-x123*x21*x98 + V{1}) + x20*(-x105*x44 - x105*x47 - x105*x57 - x105*x63 - x105*x65 + x108*x82 + x108*x87 - x123*x66 + x124 + V{0.0185185185185185}*x24*x69);
auto x4 = -x133*(-x21*x65*x98 + V{1}) - x20*(-x125*x50 - x126*x64 + x132);
auto x5 = -x133*(-x21*x63*x98 + V{1}) - x20*(x125*x56 - x126*x62 + x134);
auto x6 = -x133*(x135*x21*x98 + V{1}) - x20*(-x135*x25 - x136*x72 + x142);
auto x7 = -x133*(-x21*x82*x98 + V{1}) - x20*(-x105*x81 + x143);
auto x8 = -x133*(x144*x21*x98 + V{1}) - x20*(-x136*x75 - x144*x25 + x146);
auto x9 = -x133*(-x21*x87*x98 + V{1}) - x20*(-x105*x86 + x147);
auto x10 = -x114*(-x21*x36*x98 + V{1}) - x20*(x100*x43 + x113 + x36*x99);
auto x11 = -x114*(-x21*x41*x98 + V{1}) - x20*(x100*x46 + x122 + x41*x99);
auto x12 = -x114*(-x21*x69*x98 + V{1}) - x20*(-x105*x56 - x105*x62 - x105*x64 - x109 - x120 - x124 + V{0.037037037037037}*x21*x23*x69 + V{0.00462962962962963}*x21*x23*x81 + V{0.00462962962962963}*x21*x23*x86);
auto x13 = -x133*(-x21*x50*x98 + V{1}) - x20*(x125*x64 + x126*x50 + x132);
auto x14 = -x133*(-x21*x57*x98 + V{1}) - x20*(x125*x62 - x126*x56 + x134);
auto x15 = -x133*(-x21*x72*x98 + V{1}) - x20*(x105*x72 + x142);
auto x16 = -x133*(x150*x21*x98 + V{1}) - x20*(x136*x81 + x143 - x150*x25);
auto x17 = -x133*(-x21*x75*x98 + V{1}) - x20*(x105*x75 + x146);
auto x18 = -x133*(x152*x21*x98 + V{1}) - x20*(x136*x86 + x147 - x152*x25);
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
return { V{1}*x21*(-x43*x93 - x46*x92 - x54*x88 - x56*x96 - x64*x94 - x91*(x53 + x61) - x92*(x151 + x85) - x95*(x148 + x80) + x97), x26 + x28 + x31 };
}
};

}

}

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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x22 = V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16];
auto x23 = cell[0] + V{2}*cell[10] + cell[11] + cell[12] + cell[17] + cell[18] + cell[2] + cell[3] + cell[8] + cell[9] + x22 + V{1};
auto x24 = x21*x23;
auto x25 = V{0.0277777777777778}*x24;
auto x26 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x27 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x28 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x29 = V{1.5}*x28;
auto x30 = -x29;
auto x31 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
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
auto x48 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x49 = x48*x48;
auto x50 = x35 + x40 + V{4.5}*x49;
auto x51 = -x27;
auto x52 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x53 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x52;
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
auto x68 = V{3}*x31;
auto x69 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x70 = x39 + x69;
auto x71 = x30 + x68 + x70 + V{1};
auto x72 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x73 = V{4.5}*(x72*x72);
auto x74 = x35 + x70 + x73;
auto x75 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x76 = V{4.5}*(x75*x75);
auto x77 = x34 + x40 + x69 + x76;
auto x78 = -x69;
auto x79 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x80 = x61 + x78 - V{4.5}*x79*x79;
auto x81 = -x80;
auto x82 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x52;
auto x83 = x56 + x78 - V{4.5}*x82*x82;
auto x84 = -x83;
auto x85 = V{0.333333333333333}*x24;
auto x86 = -x55;
auto x87 = x85*x86;
auto x88 = V{0.0277777777777778}*x24;
auto x89 = V{0.0555555555555555}*x24;
auto x90 = V{0.0555555555555555}*x24;
auto x91 = V{0.0277777777777778}*x24;
auto x92 = V{0.0555555555555555}*x24;
auto x93 = V{0.0277777777777778}*x24;
auto x94 = V{4.44089209850063e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{6.66133814775094e-16}*cell[13] + V{6.66133814775094e-16}*cell[14] + V{8.88178419700125e-16}*cell[15] + V{4.44089209850063e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + V{0.111111111111111}*x24*x71 + x36*x67 + x41*x67 + x50*x93 + x74*x92 + x77*x89 + V{2.22044604925031e-16};
auto x95 = x44*x90 + x47*x89 + x58*x93 + x64*x88 + x66*x91 + x81*x92 + x84*x89 + x87 + x94;
auto x96 = x29 + x45 - x68 + x69;
auto x97 = V{0.0833333333333333}*cell[11];
auto x98 = V{0.0833333333333333}*cell[12];
auto x99 = V{0.0833333333333333}*cell[2];
auto x100 = V{0.0833333333333333}*cell[3];
auto x101 = V{0.00462962962962963}*x24;
auto x102 = x101*x46;
auto x103 = x101*x43;
auto x104 = V{0.00462962962962963}*x24;
auto x105 = x104*x74;
auto x106 = x104*x77;
auto x107 = -V{0.333333333333333}*cell[10] - V{0.166666666666667}*cell[13] - V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[15] - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x100 + x101*x57 + x101*x63 + x101*x65 + x102 + x103 + x105 + x106 - V{0.00462962962962963}*x21*x23*x36 - V{0.00462962962962963}*x21*x23*x41 - V{0.00462962962962963}*x21*x23*x50 - V{0.00462962962962963}*x21*x23*x80 - V{0.00462962962962963}*x21*x23*x83 + x97 + x98 + x99 + V{0.0555555555555555};
auto x108 = V{0.0555555555555556}*x19;
auto x109 = V{0.0462962962962963}*x24;
auto x110 = V{0.00925925925925926}*x24;
auto x111 = V{0.00925925925925926}*x24;
auto x112 = V{0.00231481481481482}*x24;
auto x113 = -V{0.166666666666667}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x111*x71 - x112*x50 + x112*x57 + x112*x63 + x112*x65 + V{-0.0555555555555555};
auto x114 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.333333333333333}*cell[15] - V{0.333333333333333}*cell[16] + V{0.166666666666667}*cell[2] - x100 + x101*x41 - x102 + x104*x80 - x105 + x111*x77 - x111*x83 + x113 - x98;
auto x115 = V{0.166666666666667}*cell[12] - V{0.333333333333333}*cell[13] - V{0.333333333333333}*cell[14] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[3] + x101*x36 - x103 + x104*x83 - x106 + x111*x74 - x111*x80 + x113 - x97 - x99;
auto x116 = x61 + x69 - x73;
auto x117 = V{0.0231481481481481}*x24;
auto x118 = V{0.00462962962962963}*x24;
auto x119 = V{0.00231481481481481}*x24;
auto x120 = V{0.00462962962962963}*x24;
auto x121 = V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] + V{0.0833333333333333}*cell[2] - x112*x77 + x112*x83 - x120*x36 + x120*x43 + V{0.0138888888888889};
auto x122 = V{0.00115740740740741}*x24;
auto x123 = V{0.166666666666667}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] - x110*x71 - x122*x50 + x122*x57 + x122*x63 + x122*x65;
auto x124 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x119*x41 - x119*x46 + x121 + x123;
auto x125 = V{0.833333333333333}*cell[13] - V{0.166666666666667}*cell[14] - x118*x80 + x124;
auto x126 = V{0.0277777777777778}*x19;
auto x127 = -x79;
auto x128 = x55 + x69;
auto x129 = x128 + x51 - V{4.5}*x127*x127;
auto x130 = -V{0.166666666666667}*cell[13] + V{0.833333333333333}*cell[14] + x118*x74 + x124;
auto x131 = x56 + x69 - x76;
auto x132 = V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] + V{0.0833333333333333}*cell[3] - x112*x74 + x112*x80 - x120*x41 + x120*x46;
auto x133 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x119*x36 - x119*x43 + x123 + x132 + V{0.0138888888888889};
auto x134 = V{0.833333333333333}*cell[15] - V{0.166666666666667}*cell[16] - x118*x83 + x133;
auto x135 = -x82;
auto x136 = x128 + x60 - V{4.5}*x135*x135;
auto x137 = -V{0.166666666666667}*cell[15] + V{0.833333333333333}*cell[16] + x118*x77 + x133;
auto x138 = V{0.0115740740740741}*x24;
auto x139 = V{0.0162037037037037}*x24;
auto x140 = V{0.00231481481481481}*x24;
auto x141 = -V{0.0833333333333333}*cell[10] + x101*x71 + x121 + x132;
auto x142 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x140*x57 - x140*x63 + x141;
auto x143 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x140*x50 - x140*x65 + x141;
auto x0 = -V{0.333333333333333}*x19*(-x21*x86*x95 + V{1}) + x20*(V{1}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[8] + V{1}*cell[9] + x22 - x25*x36 - x25*x41 - x25*x44 - x25*x47 - x25*x50 - x25*x58 - x25*x64 - x25*x66 - x67*x71 - x67*x74 - x67*x77 - x67*x81 - x67*x84 - x87 + V{0.833333333333333});
auto x1 = -x108*(x21*x95*x96 + V{1}) - x20*(-x107 - V{0.0185185185185185}*x24*x71 - x67*x96);
auto x2 = -x108*(-x21*x44*x95 + V{1}) - x20*(-x109*x43 - x110*x36 + x114);
auto x3 = -x108*(-x21*x47*x95 + V{1}) - x20*(-x109*x46 - x110*x41 + x115);
auto x4 = -x126*(x116*x21*x95 + V{1}) - x20*(-x116*x25 - x117*x74 + x125);
auto x5 = -x126*(x129*x21*x95 + V{1}) - x20*(x117*x80 - x129*x25 + x130);
auto x6 = -x126*(x131*x21*x95 + V{1}) - x20*(-x117*x77 - x131*x25 + x134);
auto x7 = -x126*(x136*x21*x95 + V{1}) - x20*(x117*x83 - x136*x25 + x137);
auto x8 = -x126*(-x21*x66*x95 + V{1}) - x20*(-x138*x50 - x139*x65 + x142);
auto x9 = -x126*(-x21*x64*x95 + V{1}) - x20*(x138*x57 - x139*x63 + x143);
auto x10 = -x108*(-x21*x71*x95 + V{1}) - x20*(-x107 + V{0.037037037037037}*x21*x23*x71);
auto x11 = -x108*(-x21*x36*x95 + V{1}) - x20*(x109*x36 + x110*x43 + x114);
auto x12 = -x108*(-x21*x41*x95 + V{1}) - x20*(x109*x41 + x110*x46 + x115);
auto x13 = -x126*(-x21*x74*x95 + V{1}) - x20*(x101*x74 + x125);
auto x14 = -x126*(-x21*x81*x95 + V{1}) - x20*(-x101*x80 + x130);
auto x15 = -x126*(-x21*x77*x95 + V{1}) - x20*(x101*x77 + x134);
auto x16 = -x126*(-x21*x84*x95 + V{1}) - x20*(-x101*x83 + x137);
auto x17 = -x126*(-x21*x50*x95 + V{1}) - x20*(x138*x65 + x139*x50 + x142);
auto x18 = -x126*(-x21*x58*x95 + V{1}) - x20*(x138*x63 - x139*x57 + x143);
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
return { V{1}*x21*(-x43*x90 - x46*x89 - x55*x85 - x57*x93 - x65*x91 - x80*x92 - x83*x89 - x88*(x54 + x62) + x94), x26 + x28 + x31 };
}
};

}

}

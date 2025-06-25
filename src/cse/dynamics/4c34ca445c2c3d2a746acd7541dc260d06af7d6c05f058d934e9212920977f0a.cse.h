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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<1, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<1, 1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{0.5}*cell[10];
auto x22 = V{0.5}*cell[12];
auto x23 = V{2}*cell[15];
auto x24 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x25 = V{0.25}*cell[0] + V{0.25}*cell[11] + V{0.5}*cell[15] + V{0.25}*cell[2] + V{0.25};
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x27 = cell[0] + cell[11] + cell[2] + x23 + V{1};
auto x28 = V{2}*cell[10] + cell[12] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[16] + cell[17] + cell[18] + cell[3] + cell[8] + cell[9] + x27;
auto x29 = cell[10] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[17] + cell[1] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + x27;
auto x30 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x31 = V{1.5}*x30;
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x33 = V{1.5}*x32;
auto x34 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x35 = V{1.5}*x34;
auto x36 = x33 + x35 + V{-1};
auto x37 = x31 + x36;
auto x38 = x37*(V{0.166666666666667}*x24*x28 + V{0.166666666666667}*x26*x29);
auto x39 = V{0.25}*x24*x28 + V{0.25}*x26*x29;
auto x40 = V{0.0138888888888889}*x24*x28 + V{0.0138888888888889}*x26*x29;
auto x41 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x42 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x43 = V{4.5}*(x42*x42);
auto x44 = -x31;
auto x45 = V{1} - x33;
auto x46 = x44 + x45;
auto x47 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x48 = -x35;
auto x49 = x47 + x48;
auto x50 = x41 + x43 + x46 + x49;
auto x51 = x37 + x47;
auto x52 = x41 - x43 + x51;
auto x53 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x54 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x53;
auto x55 = -x54;
auto x56 = -x47;
auto x57 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x58 = x37 + x57;
auto x59 = x56 + x58;
auto x60 = x59 - V{4.5}*x55*x55;
auto x61 = V{1.11022302462516e-16}*cell[12];
auto x62 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x63 = V{4.5}*(x62*x62);
auto x64 = x51 + x57 - x63;
auto x65 = -x40*x64;
auto x66 = x46 + x57;
auto x67 = x49 + x63 + x66;
auto x68 = x40*x67;
auto x69 = -x57;
auto x70 = -V{4.5}*x54*x54;
auto x71 = x51 + x69 + x70;
auto x72 = -x40*x71;
auto x73 = V{0.0277777777777778}*x24*x28 + V{0.0277777777777778}*x26*x29;
auto x74 = V{3}*x30;
auto x75 = x45 + x49 + x74;
auto x76 = x73*x75;
auto x77 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x78 = V{4.5}*(x77*x77);
auto x79 = x41 + x48;
auto x80 = x66 + x78 + x79;
auto x81 = V{0.0555555555555556}*x24*x28 + V{0.0555555555555556}*x26*x29;
auto x82 = V{3}*x32;
auto x83 = x44 + x79 + x82 + V{1};
auto x84 = V{5.55111512312578e-17}*x24*x28 + V{5.55111512312578e-17}*x26*x29;
auto x85 = V{1.11022302462516e-16}*x24*x28 + V{1.11022302462516e-16}*x26*x29;
auto x86 = V{8.32667268468867e-17}*x24*x28 + V{8.32667268468867e-17}*x26*x29;
auto x87 = x36 + x47 - x74;
auto x88 = x73*x87;
auto x89 = -x41;
auto x90 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x91 = -V{4.5}*x90*x90;
auto x92 = x58 + x89 + x91;
auto x93 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x53;
auto x94 = -V{4.5}*x93*x93;
auto x95 = x51 + x89 + x94;
auto x96 = V{3}*x34;
auto x97 = x31 + V{-1};
auto x98 = x33 + x57 - x96 + x97;
auto x99 = -x73*x98;
auto x100 = x66 + x96;
auto x101 = x100*x73;
auto x102 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[11] + V{1.11022302462516e-16}*cell[15] + V{5.55111512312578e-17}*cell[2] + V{5.55111512312578e-17};
auto x103 = x101 - x24*(V{1.11022302462516e-16}*cell[10] + V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{1.11022302462516e-16}*cell[14] + V{1.11022302462516e-16}*cell[16] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{5.55111512312578e-17}*cell[3] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x102) - x26*(V{5.55111512312578e-17}*cell[10] + V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[17] + V{5.55111512312578e-17}*cell[1] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[7] + V{1.11022302462516e-16}*cell[9] + x102 + x61) - x38 + x50*(V{0.0277777777777778}*x24*x28 + V{0.0277777777777778}*x26*x29) + x99 + V{-2.22044604925031e-16};
auto x104 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{8.88178419700125e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x103 - x30*x84 - x32*x85 - x34*x86 - x52*(V{6.16790569236198e-18}*x24*x28 + V{6.16790569236198e-18}*x26*x29) + x61 + x65 + x68 + x72 + x73*x80 - x73*x92 - x73*x95 + x76 + x81*x83 - x88;
auto x105 = -x90;
auto x106 = x37 + x41;
auto x107 = x106 + x69;
auto x108 = x107 - V{4.5}*x105*x105;
auto x109 = -x93;
auto x110 = x106 + x56;
auto x111 = x110 - V{4.5}*x109*x109;
auto x112 = x41 + x58 - x78;
auto x113 = -x112*x40;
auto x114 = x40*x80;
auto x115 = -x40*x92;
auto x116 = x73*x83;
auto x117 = x35 + x41 - x82 + x97;
auto x118 = x117*x73;
auto x119 = V{1.66533453693773e-16}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{6.66133814775094e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[17] + V{1.66533453693773e-16}*cell[1] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + x103 + x113 + x114 + x115 + x116 - x118 - x30*x85 - x32*x86 - x34*x84 - x52*(V{4.62592926927149e-18}*x24*x28 + V{4.62592926927149e-18}*x26*x29) + x67*x73 + x75*x81;
auto x120 = x24*(x104 - x40*x60) + x26*(-x108*x40 - x111*x73 + x119 - x60*x73);
auto x121 = V{0.0833333333333333}*cell[12];
auto x122 = V{0.0833333333333333}*cell[3];
auto x123 = V{0.0833333333333334}*cell[13];
auto x124 = V{0.0833333333333334}*cell[14];
auto x125 = V{0.0833333333333334}*cell[4];
auto x126 = V{0.0833333333333334}*cell[5];
auto x127 = V{0.0833333333333333}*x24*x28 + V{0.0833333333333333}*x26*x29;
auto x128 = V{0.0416666666666667}*x24*x28 + V{0.0416666666666667}*x26*x29;
auto x129 = x128*x30;
auto x130 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[2] + V{3.46944695195361e-18};
auto x131 = x24*(V{6.93889390390723e-18}*cell[10] + V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[3] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x130);
auto x132 = x26*(V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[1] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[7] + V{6.93889390390723e-18}*cell[9] + x130);
auto x133 = V{0.00115740740740741}*x24*x28 + V{0.00115740740740741}*x26*x29;
auto x134 = V{0.0833333333333333}*cell[11] - V{0.166666666666667}*cell[15] - V{0.0833333333333334}*cell[16] + V{0.0833333333333333}*cell[2] - V{0.0833333333333334}*cell[7] - x128*x34 + x131 + x132 + x133*x50 + x133*x52 + V{0.0555555555555555};
auto x135 = -V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x121 + x122 - x123 - x124 - x125 - x126 + x127*x32 - x129 + x134;
auto x136 = V{0.0833333333333334}*cell[17];
auto x137 = V{0.0833333333333334}*cell[18];
auto x138 = V{0.0833333333333334}*cell[8];
auto x139 = V{0.0833333333333334}*cell[9];
auto x140 = V{0.0833333333333333}*cell[10];
auto x141 = V{0.0833333333333333}*cell[1];
auto x142 = V{0.00231481481481481}*x24*x28 + V{0.00231481481481481}*x26*x29;
auto x143 = x128*x32;
auto x144 = V{0.166666666666667}*cell[11] - V{0.333333333333333}*cell[15] - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[2] - V{0.166666666666667}*cell[7] - x121 - x122 + x123 + x124 + x125 + x126 - x127*x34 + x129 - x131 - x132 + x136 + x137 + x138 + x139 - x140 - x141 + x142*x50 + x142*x52 + x143 + V{-0.0555555555555555};
auto x145 = -V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x127*x30 + x134 - x136 - x137 - x138 - x139 + x140 + x141 - x143;
auto x146 = x24*x28 + x26*x29;
auto x147 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x146;
auto x148 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x147;
auto x149 = V{0.0208333333333333}*x24*x28 + V{0.0208333333333333}*x26*x29;
auto x150 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0208333333333333}*cell[2] + V{0.0208333333333333};
auto x151 = -x24*(V{0.0416666666666667}*cell[10] + V{0.0208333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0416666666666667}*cell[16] + V{0.0208333333333333}*cell[17] + V{0.0208333333333333}*cell[18] + V{0.0208333333333333}*cell[3] + V{0.0208333333333333}*cell[8] + V{0.0208333333333333}*cell[9] + x150);
auto x152 = -x26*(V{0.0208333333333333}*cell[10] + V{0.0416666666666667}*cell[12] + V{0.0208333333333333}*cell[13] + V{0.0208333333333333}*cell[14] + V{0.0416666666666667}*cell[17] + V{0.0208333333333333}*cell[1] + V{0.0208333333333333}*cell[4] + V{0.0208333333333333}*cell[5] + V{0.0416666666666667}*cell[7] + V{0.0416666666666667}*cell[9] + x150);
auto x153 = V{0.0416666666666667}*x24*x28 + V{0.0416666666666667}*x26*x29;
auto x154 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x151 + x152 - x153*x32 + V{0.0138888888888889};
auto x155 = V{0.000578703703703704}*x24*x28 + V{0.000578703703703704}*x26*x29;
auto x156 = V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[7] - x153*x34 - x155*x50 - x155*x52;
auto x157 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x149*x30 + x154 + x156;
auto x158 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x148 + x157;
auto x159 = -x40*(x107 + x91);
auto x160 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x148 + x157;
auto x161 = V{0.00578703703703704}*x24*x28 + V{0.00578703703703704}*x26*x29;
auto x162 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x147;
auto x163 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x153*x30;
auto x164 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x149*x34 + x154 + x163;
auto x165 = V{0.833333333333333}*cell[15] - V{0.0833333333333333}*cell[16] - V{0.0833333333333333}*cell[7] - x162 + x164;
auto x166 = x110 + x94;
auto x167 = V{0.00115740740740741}*x24*x28 + V{0.00115740740740741}*x26*x29;
auto x168 = -V{0.166666666666667}*cell[15] + V{0.416666666666667}*cell[16] + V{0.416666666666667}*cell[7] + x162 + x164 + x167*x50 + x167*x52;
auto x169 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x146;
auto x170 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x149*x32 + x151 + x152 + x156 + x163 + V{0.0138888888888889};
auto x171 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x169 + x170;
auto x172 = x59 + x70;
auto x173 = -x172*x40;
auto x174 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x169 + x170;
auto x175 = x24*(x104 + x173) + x26*(x119 + x159 - x166*x73 - x172*x73);
auto x0 = -x19*(V{0.166666666666667}*x120*x37 + V{0.333333333333333}) + x20*(V{0.5}*cell[11] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21 + x22 + x23 - x24*(V{0.25}*cell[12] + V{0.5}*cell[13] + V{0.5}*cell[14] + V{0.5}*cell[16] + V{0.25}*cell[17] + V{0.25}*cell[18] + V{0.25}*cell[3] + V{0.25}*cell[8] + V{0.25}*cell[9] + x21 + x25) - x26*(V{0.25}*cell[10] + V{0.25}*cell[13] + V{0.25}*cell[14] + V{0.5}*cell[17] + V{0.25}*cell[1] + V{0.25}*cell[4] + V{0.25}*cell[5] + V{0.5}*cell[7] + V{0.5}*cell[9] + x22 + x25) - x30*x39 - x32*x39 - x34*x39 + x38 - x40*x50 - x40*x52 + V{0.833333333333333});
auto x1 = -x19*(V{0.0277777777777778}*x117*x120 + V{0.0555555555555556}) + x20*(x118 + x135);
auto x2 = -x19*(V{0.0277777777777778}*x120*x98 + V{0.0555555555555556}) - x20*(x144 + x99);
auto x3 = -x19*(V{0.0277777777777778}*x120*x87 + V{0.0555555555555556}) + x20*(x145 + x88);
auto x4 = -x19*(V{0.0138888888888889}*x112*x120 + V{0.0277777777777778}) - x20*(x113 + x158);
auto x5 = -x19*(V{0.0138888888888889}*x108*x120 + V{0.0277777777777778}) - x20*(x159 + x160);
auto x6 = -x19*(V{0.0138888888888889}*x120*x52 + V{0.0277777777777778}) - x20*(-x161*x50 + x165 - x52*(V{0.0196759259259259}*x24*x28 + V{0.0196759259259259}*x26*x29));
auto x7 = -x19*(V{0.0138888888888889}*x111*x120 + V{0.0277777777777778}) - x20*(-x166*x40 + x168);
auto x8 = -x19*(V{0.0138888888888889}*x120*x64 + V{0.0277777777777778}) - x20*(x171 + x65);
auto x9 = -x19*(V{0.0138888888888889}*x120*x60 + V{0.0277777777777778}) - x20*(x173 + x174);
auto x10 = x19*(V{0.0277777777777778}*x175*x83 + V{-0.0555555555555556}) + x20*(-x116 + x135);
auto x11 = x19*(V{0.0277777777777778}*x100*x175 + V{-0.0555555555555556}) - x20*(x101 + x144);
auto x12 = x19*(V{0.0277777777777778}*x175*x75 + V{-0.0555555555555556}) + x20*(x145 - x76);
auto x13 = x19*(V{0.0138888888888889}*x175*x80 + V{-0.0277777777777778}) - x20*(x114 + x158);
auto x14 = -x19*(V{0.0138888888888889}*x120*x92 + V{0.0277777777777778}) - x20*(x115 + x160);
auto x15 = x19*(V{0.0138888888888889}*x175*x50 + V{-0.0277777777777778}) - x20*(-x161*x52 + x165 + x50*(V{0.00810185185185185}*x24*x28 + V{0.00810185185185185}*x26*x29));
auto x16 = -x19*(V{0.0138888888888889}*x120*x95 + V{0.0277777777777778}) - x20*(x168 - x40*x95);
auto x17 = x19*(V{0.0138888888888889}*x175*x67 + V{-0.0277777777777778}) - x20*(x171 + x68);
auto x18 = -x19*(V{0.0138888888888889}*x120*x71 + V{0.0277777777777778}) - x20*(x174 + x72);
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
return { V{0.5}*x175, x30 + x32 + x34 };
}
};

}

}

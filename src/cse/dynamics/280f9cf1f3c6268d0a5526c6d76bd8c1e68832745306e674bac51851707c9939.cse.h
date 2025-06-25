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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<1, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<1, 1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x22 = V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9];
auto x23 = V{0.166666666666667}*cell[0];
auto x24 = V{0.166666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.166666666666667}*cell[3] + x23 + V{0.166666666666667};
auto x25 = V{0.166666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.166666666666667}*cell[2];
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x27 = V{0.166666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.166666666666667}*cell[1];
auto x28 = V{0.166666666666667}*cell[15];
auto x29 = V{0.166666666666667}*cell[16];
auto x30 = V{0.166666666666667}*cell[6];
auto x31 = V{0.166666666666667}*cell[7];
auto x32 = x28 + x29 + x30 + x31;
auto x33 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x34 = V{0.166666666666667}*cell[13];
auto x35 = V{0.166666666666667}*cell[14];
auto x36 = V{0.166666666666667}*cell[4];
auto x37 = V{0.166666666666667}*cell[5];
auto x38 = x34 + x35 + x36 + x37;
auto x39 = cell[0] + cell[12] + V{2}*cell[13] + cell[3] + V{1};
auto x40 = cell[11] + V{2}*cell[15] + cell[2];
auto x41 = V{2}*cell[10] + V{2}*cell[14] + V{2}*cell[16] + cell[17] + cell[18] + cell[8] + cell[9] + x39 + x40;
auto x42 = cell[10] + V{2}*cell[17] + cell[1];
auto x43 = V{2}*cell[11] + cell[15] + cell[16] + V{2}*cell[18] + V{2}*cell[5] + cell[6] + cell[7] + x39 + x42;
auto x44 = cell[0] + V{2}*cell[12] + cell[13] + cell[14] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + x40 + x42 + V{1};
auto x45 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x46 = V{1.5}*x45;
auto x47 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x48 = V{1.5}*x47;
auto x49 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x50 = V{1.5}*x49;
auto x51 = x48 + x50 + V{-1};
auto x52 = x46 + x51;
auto x53 = x52*(V{0.111111111111111}*x21*x41 + V{0.111111111111111}*x26*x43 + V{0.111111111111111}*x33*x44);
auto x54 = V{0.166666666666667}*x21*x41 + V{0.166666666666667}*x26*x43 + V{0.166666666666667}*x33*x44;
auto x55 = V{0.00925925925925926}*x21*x41 + V{0.00925925925925926}*x26*x43 + V{0.00925925925925926}*x33*x44;
auto x56 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x57 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x56;
auto x58 = -x57;
auto x59 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x60 = -x59;
auto x61 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x62 = x52 + x61;
auto x63 = x60 + x62;
auto x64 = x63 - V{4.5}*x58*x58;
auto x65 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x66 = V{4.5}*(x65*x65);
auto x67 = x59 + x62 - x66;
auto x68 = -x55*x67;
auto x69 = -x46;
auto x70 = V{1} - x48;
auto x71 = x69 + x70;
auto x72 = x61 + x71;
auto x73 = -x50;
auto x74 = x59 + x73;
auto x75 = x66 + x72 + x74;
auto x76 = x55*x75;
auto x77 = -x61;
auto x78 = -V{4.5}*x57*x57;
auto x79 = x52 + x59;
auto x80 = x77 + x78 + x79;
auto x81 = -x55*x80;
auto x82 = V{0.037037037037037}*x21*x41 + V{0.037037037037037}*x26*x43 + V{0.037037037037037}*x33*x44;
auto x83 = V{3}*x47;
auto x84 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x85 = x73 + x84;
auto x86 = x69 + x83 + x85 + V{1};
auto x87 = V{5.55111512312578e-17}*x21*x41 + V{5.55111512312578e-17}*x26*x43 + V{5.55111512312578e-17}*x33*x44;
auto x88 = V{0.0185185185185185}*x21*x41 + V{0.0185185185185185}*x26*x43 + V{0.0185185185185185}*x33*x44;
auto x89 = -x84;
auto x90 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x91 = -V{4.5}*x90*x90;
auto x92 = x62 + x89 + x91;
auto x93 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x56;
auto x94 = -V{4.5}*x93*x93;
auto x95 = x79 + x89 + x94;
auto x96 = V{7.40148683083438e-17}*x21*x41 + V{7.40148683083438e-17}*x26*x43 + V{7.40148683083438e-17}*x33*x44;
auto x97 = V{3}*x49;
auto x98 = x72 + x97;
auto x99 = x88*x98;
auto x100 = V{5.55111512312578e-17}*cell[0];
auto x101 = V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{5.55111512312578e-17}*cell[3] + x100 + V{5.55111512312578e-17};
auto x102 = V{5.55111512312578e-17}*cell[11] + V{1.11022302462516e-16}*cell[15] + V{5.55111512312578e-17}*cell[2];
auto x103 = -x21*(V{1.11022302462516e-16}*cell[10] + V{1.11022302462516e-16}*cell[14] + V{1.11022302462516e-16}*cell[16] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x101 + x102);
auto x104 = V{1.11022302462516e-16}*cell[11];
auto x105 = V{5.55111512312578e-17}*cell[10] + V{1.11022302462516e-16}*cell[17] + V{5.55111512312578e-17}*cell[1];
auto x106 = -x26*(V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{1.11022302462516e-16}*cell[18] + V{1.11022302462516e-16}*cell[5] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + x101 + x104 + x105);
auto x107 = V{1.11022302462516e-16}*cell[12];
auto x108 = -x33*(V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[7] + V{1.11022302462516e-16}*cell[9] + x100 + x102 + x105 + x107 + V{5.55111512312578e-17});
auto x109 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x110 = V{4.5}*(x109*x109);
auto x111 = x110 + x71 + x74 + x84;
auto x112 = x46 + V{-1};
auto x113 = x112 + x48 + x61 - x97;
auto x114 = x113*x88;
auto x115 = -x53;
auto x116 = x103 + x106 + x108 + x111*x88 - x114 + x115 + x99 + V{-2.22044604925031e-16};
auto x117 = V{3}*x45;
auto x118 = x117 + x70 + x74;
auto x119 = x118*x88;
auto x120 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x121 = V{4.5}*(x120*x120);
auto x122 = x121 + x72 + x85;
auto x123 = -x117 + x51 + x59;
auto x124 = x123*x88;
auto x125 = V{3.70074341541719e-17}*x21*x41 + V{3.70074341541719e-17}*x26*x43 + V{3.70074341541719e-17}*x33*x44;
auto x126 = V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + x107 + x119 + x122*x88 - x124 - x125*x45;
auto x127 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x116 + x126 - x47*x96 - x49*x87 + x68 + x76 + x81 + x82*x86 - x88*x92 - x88*x95;
auto x128 = -x93;
auto x129 = x52 + x84;
auto x130 = x129 + x60;
auto x131 = x130 - V{4.5}*x128*x128;
auto x132 = -x90;
auto x133 = x129 + x77;
auto x134 = x133 - V{4.5}*x132*x132;
auto x135 = -x110 + x79 + x84;
auto x136 = -x135*x55;
auto x137 = x111*x55;
auto x138 = -x55*x95;
auto x139 = x86*x88;
auto x140 = x112 + x50 - x83 + x84;
auto x141 = x140*x88;
auto x142 = V{1.66533453693773e-16}*cell[10] + V{1.66533453693773e-16}*cell[1] + x139 - x141 - x47*x87 + x75*x88;
auto x143 = V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] + x103 + x106 + x108 + x115 + x126 + x136 + x137 + x138 + x142 - x49*x96 - x80*x88 + x82*x98 + V{-2.22044604925031e-16};
auto x144 = -x121 + x62 + x84;
auto x145 = -x144*x55;
auto x146 = x122*x55;
auto x147 = -x55*x92;
auto x148 = V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[17] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + x104 + x116 + x118*x82 - x125*x49 + x142 + x145 + x146 + x147 - x45*x96;
auto x149 = x21*(x127 - x55*x64) + x26*(-x131*x55 - x134*x88 + x143) + x33*(-x131*x88 - x134*x55 + x148 - x64*x88);
auto x150 = V{0.0555555555555556}*x21*x41 + V{0.0555555555555556}*x26*x43 + V{0.0555555555555556}*x33*x44;
auto x151 = V{0.0833333333333333}*cell[12];
auto x152 = V{0.0833333333333333}*cell[3];
auto x153 = V{0.0833333333333334}*cell[13];
auto x154 = V{0.0833333333333334}*cell[14];
auto x155 = V{0.0833333333333334}*cell[4];
auto x156 = V{0.0833333333333334}*cell[5];
auto x157 = V{3.46944695195361e-18}*cell[0];
auto x158 = V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[3] + x157 + V{3.46944695195361e-18};
auto x159 = V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[2];
auto x160 = x21*(V{6.93889390390723e-18}*cell[10] + V{6.93889390390723e-18}*cell[14] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x158 + x159);
auto x161 = V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[1];
auto x162 = x26*(V{6.93889390390723e-18}*cell[11] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[18] + V{6.93889390390723e-18}*cell[5] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x158 + x161);
auto x163 = x33*(V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[7] + V{6.93889390390723e-18}*cell[9] + x157 + x159 + x161 + V{3.46944695195361e-18});
auto x164 = V{0.0277777777777778}*x21*x41 + V{0.0277777777777778}*x26*x43 + V{0.0277777777777778}*x33*x44;
auto x165 = x164*x45;
auto x166 = x151 + x152 - x153 - x154 - x155 - x156 + x160 + x162 + x163 - x165 + V{0.0555555555555555};
auto x167 = V{0.0833333333333333}*cell[11];
auto x168 = V{0.0833333333333333}*cell[2];
auto x169 = V{0.0833333333333334}*cell[15];
auto x170 = V{0.0833333333333334}*cell[16];
auto x171 = V{0.0833333333333334}*cell[6];
auto x172 = V{0.0833333333333334}*cell[7];
auto x173 = x164*x49;
auto x174 = x167 + x168 - x169 - x170 - x171 - x172 - x173;
auto x175 = -V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[1] + x150*x47 + x166 + x174 + x22;
auto x176 = V{0.166666666666667}*cell[11];
auto x177 = V{0.166666666666667}*cell[2];
auto x178 = x150*x49;
auto x179 = V{0.0833333333333333}*cell[10];
auto x180 = V{0.0833333333333333}*cell[1];
auto x181 = V{0.0833333333333334}*cell[17];
auto x182 = V{0.0833333333333334}*cell[18];
auto x183 = V{0.0833333333333334}*cell[8];
auto x184 = V{0.0833333333333334}*cell[9];
auto x185 = x164*x47;
auto x186 = x179 + x180 - x181 - x182 - x183 - x184 - x185;
auto x187 = V{0.166666666666667}*cell[12];
auto x188 = V{0.166666666666667}*cell[3];
auto x189 = x150*x45;
auto x190 = x21*x41 + x26*x43 + x33*x44;
auto x191 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x190;
auto x192 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x191;
auto x193 = V{0.0138888888888889}*x21*x41 + V{0.0138888888888889}*x26*x43 + V{0.0138888888888889}*x33*x44;
auto x194 = V{0.0138888888888889}*cell[0];
auto x195 = V{0.0138888888888889}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0138888888888889}*cell[3] + x194 + V{0.0138888888888889};
auto x196 = V{0.0138888888888889}*cell[11] + V{0.0277777777777778}*cell[15] + V{0.0138888888888889}*cell[2];
auto x197 = -x21*(V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[16] + V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + x195 + x196);
auto x198 = V{0.0138888888888889}*cell[10] + V{0.0277777777777778}*cell[17] + V{0.0138888888888889}*cell[1];
auto x199 = -x26*(V{0.0277777777777778}*cell[11] + V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[5] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + x195 + x198);
auto x200 = -x33*(V{0.0277777777777778}*cell[12] + V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[9] + x194 + x196 + x198 + V{0.0138888888888889});
auto x201 = V{0.0277777777777778}*x21*x41 + V{0.0277777777777778}*x26*x43 + V{0.0277777777777778}*x33*x44;
auto x202 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x197 + x199 + x200 - x201*x47 + V{0.0138888888888889};
auto x203 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x201*x49;
auto x204 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x193*x45 + x202 + x203;
auto x205 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x192 + x204;
auto x206 = x133 + x91;
auto x207 = -x206*x55;
auto x208 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x192 + x204;
auto x209 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x191;
auto x210 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x201*x45;
auto x211 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x193*x49 + x202 + x210;
auto x212 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x209 + x211;
auto x213 = x130 + x94;
auto x214 = -x213*x55;
auto x215 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x209 + x211;
auto x216 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x190;
auto x217 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x193*x47 + x197 + x199 + x200 + x203 + x210 + V{0.0138888888888889};
auto x218 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x216 + x217;
auto x219 = x63 + x78;
auto x220 = -x219*x55;
auto x221 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x216 + x217;
auto x222 = x21*(x127 + x220) + x26*(x143 - x206*x88 + x214) + x33*(x148 + x207 - x213*x88 - x219*x88);
auto x223 = -x160 - x162 - x163 - x179 - x180 + x181 + x182 + x183 + x184 + x185 + V{-0.0555555555555555};
auto x0 = -x19*(V{0.111111111111111}*x149*x52 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] - x21*(V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[16] + x22 + x24 + x25) - x26*(V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[5] + x24 + x27 + x32) - x33*(V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[9] + x23 + x25 + x27 + x38 + V{0.166666666666667}) - x45*x54 - x47*x54 - x49*x54 + x53 + V{0.833333333333333});
auto x1 = -x19*(V{0.0185185185185185}*x140*x149 + V{0.0555555555555556}) + x20*(x141 + x175);
auto x2 = -x19*(V{0.0185185185185185}*x113*x149 + V{0.0555555555555556}) + x20*(x114 + x166 - x176 - x177 + x178 + x186 + x32);
auto x3 = -x19*(V{0.0185185185185185}*x123*x149 + V{0.0555555555555556}) + x20*(x124 + x160 + x162 + x163 + x174 + x186 - x187 - x188 + x189 + x38 + V{0.0555555555555555});
auto x4 = -x19*(V{0.00925925925925926}*x144*x149 + V{0.0277777777777778}) - x20*(x145 + x205);
auto x5 = -x19*(V{0.00925925925925926}*x134*x149 + V{0.0277777777777778}) - x20*(x207 + x208);
auto x6 = -x19*(V{0.00925925925925926}*x135*x149 + V{0.0277777777777778}) - x20*(x136 + x212);
auto x7 = -x19*(V{0.00925925925925926}*x131*x149 + V{0.0277777777777778}) - x20*(x214 + x215);
auto x8 = -x19*(V{0.00925925925925926}*x149*x67 + V{0.0277777777777778}) - x20*(x218 + x68);
auto x9 = -x19*(V{0.00925925925925926}*x149*x64 + V{0.0277777777777778}) - x20*(x220 + x221);
auto x10 = x19*(V{0.0185185185185185}*x222*x86 + V{-0.0555555555555556}) + x20*(-x139 + x175);
auto x11 = x19*(V{0.0185185185185185}*x222*x98 + V{-0.0555555555555556}) - x20*(-x151 - x152 + x153 + x154 + x155 + x156 + x165 + x176 + x177 - x178 + x223 - x28 - x29 - x30 - x31 + x99);
auto x12 = x19*(V{0.0185185185185185}*x118*x222 + V{-0.0555555555555556}) - x20*(x119 - x167 - x168 + x169 + x170 + x171 + x172 + x173 + x187 + x188 - x189 + x223 - x34 - x35 - x36 - x37);
auto x13 = x19*(V{0.00925925925925926}*x122*x222 + V{-0.0277777777777778}) - x20*(x146 + x205);
auto x14 = -x19*(V{0.00925925925925926}*x149*x92 + V{0.0277777777777778}) - x20*(x147 + x208);
auto x15 = x19*(V{0.00925925925925926}*x111*x222 + V{-0.0277777777777778}) - x20*(x137 + x212);
auto x16 = -x19*(V{0.00925925925925926}*x149*x95 + V{0.0277777777777778}) - x20*(x138 + x215);
auto x17 = x19*(V{0.00925925925925926}*x222*x75 + V{-0.0277777777777778}) - x20*(x218 + x76);
auto x18 = -x19*(V{0.00925925925925926}*x149*x80 + V{0.0277777777777778}) - x20*(x221 + x81);
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
return { V{0.333333333333333}*x222, x45 + x47 + x49 };
}
};

}

}

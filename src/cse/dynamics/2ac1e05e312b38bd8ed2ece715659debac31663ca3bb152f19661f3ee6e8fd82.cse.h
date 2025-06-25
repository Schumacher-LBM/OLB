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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<1, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<1, 1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1});
auto x22 = V{0.166666666666667}*cell[13];
auto x23 = V{0.166666666666667}*cell[14];
auto x24 = V{0.166666666666667}*cell[4];
auto x25 = V{0.166666666666667}*cell[5];
auto x26 = V{0.166666666666667}*cell[0];
auto x27 = V{0.166666666666667}*cell[11] + V{0.333333333333333}*cell[16] + V{0.166666666666667}*cell[2] + x26 + V{0.166666666666667};
auto x28 = V{0.166666666666667}*cell[10] + V{0.333333333333333}*cell[18] + V{0.166666666666667}*cell[1];
auto x29 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x30 = V{0.166666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.166666666666667}*cell[3];
auto x31 = V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9];
auto x32 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x33 = V{0.166666666666667}*cell[15];
auto x34 = V{0.166666666666667}*cell[16];
auto x35 = V{0.166666666666667}*cell[6];
auto x36 = V{0.166666666666667}*cell[7];
auto x37 = cell[0] + cell[12] + V{2}*cell[13] + cell[3] + V{1};
auto x38 = cell[11] + V{2}*cell[16] + cell[2];
auto x39 = V{2}*cell[10] + V{2}*cell[14] + V{2}*cell[15] + cell[17] + cell[18] + cell[8] + cell[9] + x37 + x38;
auto x40 = cell[10] + V{2}*cell[18] + cell[1];
auto x41 = V{2}*cell[11] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[5] + cell[6] + cell[7] + x37 + x40;
auto x42 = cell[0] + cell[13] + cell[14] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + V{2}*cell[8] + x38 + x40 + V{1};
auto x43 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x44 = V{1.5}*x43;
auto x45 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x46 = V{1.5}*x45;
auto x47 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x48 = V{1.5}*x47;
auto x49 = x46 + x48 + V{-1};
auto x50 = x44 + x49;
auto x51 = x50*(-V{0.111111111111111}*x21*x42 + V{0.111111111111111}*x29*x39 + V{0.111111111111111}*x32*x41);
auto x52 = -V{0.166666666666667}*x21*x42 + V{0.166666666666667}*x29*x39 + V{0.166666666666667}*x32*x41;
auto x53 = -V{0.00925925925925926}*x21*x42 + V{0.00925925925925926}*x29*x39 + V{0.00925925925925926}*x32*x41;
auto x54 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x55 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x54;
auto x56 = -x55;
auto x57 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x58 = -x57;
auto x59 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x60 = x50 + x59;
auto x61 = x58 + x60;
auto x62 = x61 - V{4.5}*x56*x56;
auto x63 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x64 = V{4.5}*(x63*x63);
auto x65 = x57 + x60 - x64;
auto x66 = -x53*x65;
auto x67 = -x44;
auto x68 = V{1} - x46;
auto x69 = x67 + x68;
auto x70 = x59 + x69;
auto x71 = -x48;
auto x72 = x57 + x71;
auto x73 = x64 + x70 + x72;
auto x74 = x53*x73;
auto x75 = -x59;
auto x76 = -V{4.5}*x55*x55;
auto x77 = x50 + x57;
auto x78 = x75 + x76 + x77;
auto x79 = -x53*x78;
auto x80 = -V{0.0185185185185185}*x21*x42 + V{0.0185185185185185}*x29*x39 + V{0.0185185185185185}*x32*x41;
auto x81 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x82 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x83 = V{4.5}*(x82*x82);
auto x84 = x69 + x72 + x81 + x83;
auto x85 = -V{0.037037037037037}*x21*x42 + V{0.037037037037037}*x29*x39 + V{0.037037037037037}*x32*x41;
auto x86 = V{3}*x45;
auto x87 = x71 + x81;
auto x88 = x67 + x86 + x87 + V{1};
auto x89 = -V{5.55111512312578e-17}*x21*x42 + V{5.55111512312578e-17}*x29*x39 + V{5.55111512312578e-17}*x32*x41;
auto x90 = -x81;
auto x91 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x92 = -V{4.5}*x91*x91;
auto x93 = x60 + x90 + x92;
auto x94 = -V{7.40148683083438e-17}*x21*x42 + V{7.40148683083438e-17}*x29*x39 + V{7.40148683083438e-17}*x32*x41;
auto x95 = V{3}*x47;
auto x96 = x44 + V{-1};
auto x97 = x46 + x59 - x95 + x96;
auto x98 = -x80*x97;
auto x99 = x70 + x95;
auto x100 = x80*x99;
auto x101 = V{1.11022302462516e-16}*cell[3];
auto x102 = V{5.55111512312578e-17}*cell[0];
auto x103 = V{5.55111512312578e-17}*cell[11] + V{1.11022302462516e-16}*cell[16] + V{5.55111512312578e-17}*cell[2] + x102 + V{5.55111512312578e-17};
auto x104 = V{5.55111512312578e-17}*cell[10] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[1];
auto x105 = x21*(V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[6] + V{1.11022302462516e-16}*cell[8] + x101 + x103 + x104);
auto x106 = V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{5.55111512312578e-17}*cell[3];
auto x107 = -x29*(V{1.11022302462516e-16}*cell[10] + V{1.11022302462516e-16}*cell[14] + V{1.11022302462516e-16}*cell[15] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x103 + x106);
auto x108 = V{1.11022302462516e-16}*cell[11];
auto x109 = -x32*(V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{1.11022302462516e-16}*cell[17] + V{1.11022302462516e-16}*cell[5] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + x102 + x104 + x106 + x108 + V{5.55111512312578e-17});
auto x110 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x54;
auto x111 = -V{4.5}*x110*x110;
auto x112 = x111 + x77 + x90;
auto x113 = -x51;
auto x114 = x100 + x105 + x107 + x109 - x112*x80 + x113 + x98 + V{-2.22044604925031e-16};
auto x115 = V{3}*x43;
auto x116 = -x115 + x49 + x57;
auto x117 = -x116*x80;
auto x118 = x115 + x68 + x72;
auto x119 = x118*x80;
auto x120 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x121 = V{4.5}*(x120*x120);
auto x122 = x121 + x70 + x87;
auto x123 = -V{3.70074341541719e-17}*x21*x42 + V{3.70074341541719e-17}*x29*x39 + V{3.70074341541719e-17}*x32*x41;
auto x124 = V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + x101 + x117 + x119 + x122*x80 - x123*x43;
auto x125 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x114 + x124 - x45*x94 - x47*x89 + x66 + x74 + x79 + x80*x84 - x80*x93 + x85*x88;
auto x126 = -x110;
auto x127 = x50 + x81;
auto x128 = x127 + x58;
auto x129 = x128 - V{4.5}*x126*x126;
auto x130 = -x91;
auto x131 = x127 + x75;
auto x132 = x131 - V{4.5}*x130*x130;
auto x133 = x77 + x81 - x83;
auto x134 = -x133*x53;
auto x135 = x53*x84;
auto x136 = -x112*x53;
auto x137 = x80*x88;
auto x138 = x48 + x81 - x86 + x96;
auto x139 = x138*x80;
auto x140 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[9] + x137 - x139 - x45*x89 - x78*x80;
auto x141 = V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + x105 + x107 + x109 + x113 + x124 + x134 + x135 + x136 + x140 - x47*x94 + x73*x80 + x85*x99 + V{-2.22044604925031e-16};
auto x142 = -x121 + x60 + x81;
auto x143 = -x142*x53;
auto x144 = x122*x53;
auto x145 = -x53*x93;
auto x146 = V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x108 + x114 - x116*x85 - x123*x47 - x133*x80 + x140 + x143 + x144 + x145 - x43*x94 - x65*x80;
auto x147 = -x21*(-x132*x53 + x146) + x29*(x125 - x53*x62) + x32*(-x129*x53 - x132*x80 + x141);
auto x148 = V{0.0833333333333333}*cell[11];
auto x149 = V{0.0833333333333333}*cell[12];
auto x150 = V{0.0833333333333333}*cell[2];
auto x151 = V{0.0833333333333333}*cell[3];
auto x152 = V{0.166666666666667}*cell[10];
auto x153 = V{0.166666666666667}*cell[1];
auto x154 = V{0.0833333333333334}*cell[13];
auto x155 = V{0.0833333333333334}*cell[14];
auto x156 = V{0.0833333333333334}*cell[15];
auto x157 = V{0.0833333333333334}*cell[16];
auto x158 = V{0.0833333333333334}*cell[4];
auto x159 = V{0.0833333333333334}*cell[5];
auto x160 = V{0.0833333333333334}*cell[6];
auto x161 = V{0.0833333333333334}*cell[7];
auto x162 = V{3.46944695195361e-18}*cell[0];
auto x163 = V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[3] + x162 + V{3.46944695195361e-18};
auto x164 = V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[2];
auto x165 = x29*(V{6.93889390390723e-18}*cell[10] + V{6.93889390390723e-18}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x163 + x164);
auto x166 = V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[1];
auto x167 = x32*(V{6.93889390390723e-18}*cell[11] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[5] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x163 + x166);
auto x168 = x21*(V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[3] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[8] + x162 + x164 + x166 + V{3.46944695195361e-18});
auto x169 = -V{0.0555555555555556}*x21*x42 + V{0.0555555555555556}*x29*x39 + V{0.0555555555555556}*x32*x41;
auto x170 = -V{0.0277777777777778}*x21*x42 + V{0.0277777777777778}*x29*x39 + V{0.0277777777777778}*x32*x41;
auto x171 = x170*x47;
auto x172 = x170*x43;
auto x173 = -x165;
auto x174 = -x167;
auto x175 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x168 + x170*x45 + x173 + x174 + V{-0.0555555555555555};
auto x176 = -x149 - x151 + x154 + x155 + x158 + x159 + x172;
auto x177 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] - x169*x47 + x175 + x176 - x33 - x34 - x35 - x36;
auto x178 = -x148 - x150 + x156 + x157 + x160 + x161 + x171;
auto x179 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] - x169*x43 + x175 + x178 - x22 - x23 - x24 - x25;
auto x180 = -x21*x42 + x29*x39 + x32*x41;
auto x181 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x180;
auto x182 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x181;
auto x183 = -V{0.0138888888888889}*x21*x42 + V{0.0138888888888889}*x29*x39 + V{0.0138888888888889}*x32*x41;
auto x184 = V{0.0138888888888889}*cell[0];
auto x185 = V{0.0138888888888889}*cell[11] + V{0.0277777777777778}*cell[16] + V{0.0138888888888889}*cell[2] + x184 + V{0.0138888888888889};
auto x186 = V{0.0138888888888889}*cell[10] + V{0.0277777777777778}*cell[18] + V{0.0138888888888889}*cell[1];
auto x187 = x21*(V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0277777777777778}*cell[3] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[8] + x185 + x186);
auto x188 = V{0.0138888888888889}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0138888888888889}*cell[3];
auto x189 = -x29*(V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + x185 + x188);
auto x190 = -x32*(V{0.0277777777777778}*cell[11] + V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[5] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + x184 + x186 + x188 + V{0.0138888888888889});
auto x191 = -V{0.0277777777777778}*x21*x42 + V{0.0277777777777778}*x29*x39 + V{0.0277777777777778}*x32*x41;
auto x192 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x187 + x189 + x190 - x191*x45 + V{0.0138888888888889};
auto x193 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x191*x47;
auto x194 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x183*x43 + x192 + x193;
auto x195 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x182 + x194;
auto x196 = x131 + x92;
auto x197 = -x196*x53;
auto x198 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x182 + x194;
auto x199 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x181;
auto x200 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x191*x43;
auto x201 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x183*x47 + x192 + x200;
auto x202 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x199 + x201;
auto x203 = -x53*(x111 + x128);
auto x204 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x199 + x201;
auto x205 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x180;
auto x206 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x183*x45 + x187 + x189 + x190 + x193 + x200 + V{0.0138888888888889};
auto x207 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x205 + x206;
auto x208 = -x53*(x61 + x76);
auto x209 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x205 + x206;
auto x210 = -x21*(x146 + x197) + x29*(x125 + x208) + x32*(x141 - x196*x80 + x203);
auto x0 = -x19*(V{0.111111111111111}*x147*x50 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21*(V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[8] + x22 + x23 + x24 + x25 + x27 + x28) - x29*(V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + x27 + x30 + x31) - x32*(V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[5] + x26 + x28 + x30 + x33 + x34 + x35 + x36 + V{0.166666666666667}) - x43*x52 - x45*x52 - x47*x52 + x51 + V{0.833333333333333});
auto x1 = -x19*(V{0.0185185185185185}*x138*x147 + V{0.0555555555555556}) + x20*(x139 + x148 + x149 + x150 + x151 - x152 - x153 - x154 - x155 - x156 - x157 - x158 - x159 - x160 - x161 + x165 + x167 - x168 + x169*x45 - x171 - x172 + x31 + V{0.0555555555555555});
auto x2 = -x19*(V{0.0185185185185185}*x147*x97 + V{0.0555555555555556}) - x20*(x177 + x98);
auto x3 = -x19*(V{0.0185185185185185}*x116*x147 + V{0.0555555555555556}) - x20*(x117 + x179);
auto x4 = -x19*(V{0.00925925925925926}*x142*x147 + V{0.0277777777777778}) - x20*(x143 + x195);
auto x5 = -x19*(V{0.00925925925925926}*x132*x147 + V{0.0277777777777778}) - x20*(x197 + x198);
auto x6 = -x19*(V{0.00925925925925926}*x133*x147 + V{0.0277777777777778}) - x20*(x134 + x202);
auto x7 = -x19*(V{0.00925925925925926}*x129*x147 + V{0.0277777777777778}) - x20*(x203 + x204);
auto x8 = -x19*(V{0.00925925925925926}*x147*x65 + V{0.0277777777777778}) - x20*(x207 + x66);
auto x9 = -x19*(V{0.00925925925925926}*x147*x62 + V{0.0277777777777778}) - x20*(x208 + x209);
auto x10 = x19*(V{0.0185185185185185}*x210*x88 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x137 - x152 - x153 - x168 - x173 - x174 - x176 - x178 + V{0.0555555555555556}*x180*x45 + V{0.0555555555555555});
auto x11 = x19*(V{0.0185185185185185}*x210*x99 + V{-0.0555555555555556}) - x20*(x100 + x177);
auto x12 = x19*(V{0.0185185185185185}*x118*x210 + V{-0.0555555555555556}) - x20*(x119 + x179);
auto x13 = x19*(V{0.00925925925925926}*x122*x210 + V{-0.0277777777777778}) - x20*(x144 + x195);
auto x14 = -x19*(V{0.00925925925925926}*x147*x93 + V{0.0277777777777778}) - x20*(x145 + x198);
auto x15 = x19*(V{0.00925925925925926}*x210*x84 + V{-0.0277777777777778}) - x20*(x135 + x202);
auto x16 = -x19*(V{0.00925925925925926}*x112*x147 + V{0.0277777777777778}) - x20*(x136 + x204);
auto x17 = x19*(V{0.00925925925925926}*x210*x73 + V{-0.0277777777777778}) - x20*(x207 + x74);
auto x18 = -x19*(V{0.00925925925925926}*x147*x78 + V{0.0277777777777778}) - x20*(x209 + x79);
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
return { V{0.333333333333333}*x210, x43 + x45 + x47 };
}
};

}

}

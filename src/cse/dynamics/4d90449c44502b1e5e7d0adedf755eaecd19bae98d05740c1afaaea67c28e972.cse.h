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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<1, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<1, -1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1});
auto x22 = V{0.166666666666667}*cell[15];
auto x23 = V{0.166666666666667}*cell[16];
auto x24 = V{0.166666666666667}*cell[6];
auto x25 = V{0.166666666666667}*cell[7];
auto x26 = V{0.166666666666667}*cell[0];
auto x27 = V{0.166666666666667}*cell[12] + V{0.333333333333333}*cell[14] + V{0.166666666666667}*cell[3] + x26 + V{0.166666666666667};
auto x28 = V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] + V{0.333333333333333}*cell[9];
auto x29 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x30 = V{0.166666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.166666666666667}*cell[2];
auto x31 = V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9];
auto x32 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x33 = V{0.166666666666667}*cell[13];
auto x34 = V{0.166666666666667}*cell[14];
auto x35 = V{0.166666666666667}*cell[4];
auto x36 = V{0.166666666666667}*cell[5];
auto x37 = cell[0] + cell[11] + V{2}*cell[15] + cell[2] + V{1};
auto x38 = cell[12] + V{2}*cell[14] + cell[3];
auto x39 = V{2}*cell[10] + V{2}*cell[13] + V{2}*cell[16] + cell[17] + cell[18] + cell[8] + cell[9] + x37 + x38;
auto x40 = cell[10] + cell[1] + V{2}*cell[9];
auto x41 = V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[17] + cell[4] + cell[5] + V{2}*cell[7] + x37 + x40;
auto x42 = cell[0] + cell[15] + cell[16] + V{2}*cell[2] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[8] + x38 + x40 + V{1};
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
auto x63 = V{1.11022302462516e-16}*cell[12];
auto x64 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x65 = V{4.5}*(x64*x64);
auto x66 = x57 + x60 - x65;
auto x67 = -x53*x66;
auto x68 = -x44;
auto x69 = V{1} - x46;
auto x70 = x68 + x69;
auto x71 = x59 + x70;
auto x72 = -x48;
auto x73 = x57 + x72;
auto x74 = x65 + x71 + x73;
auto x75 = x53*x74;
auto x76 = -x59;
auto x77 = -V{4.5}*x55*x55;
auto x78 = x50 + x57;
auto x79 = x76 + x77 + x78;
auto x80 = -x53*x79;
auto x81 = -V{0.0185185185185185}*x21*x42 + V{0.0185185185185185}*x29*x39 + V{0.0185185185185185}*x32*x41;
auto x82 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x83 = V{4.5}*(x82*x82);
auto x84 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x85 = x72 + x84;
auto x86 = x71 + x83 + x85;
auto x87 = -V{0.037037037037037}*x21*x42 + V{0.037037037037037}*x29*x39 + V{0.037037037037037}*x32*x41;
auto x88 = V{3}*x45;
auto x89 = x68 + x85 + x88 + V{1};
auto x90 = -V{5.55111512312578e-17}*x21*x42 + V{5.55111512312578e-17}*x29*x39 + V{5.55111512312578e-17}*x32*x41;
auto x91 = -x84;
auto x92 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x54;
auto x93 = -V{4.5}*x92*x92;
auto x94 = x78 + x91 + x93;
auto x95 = -V{3.70074341541719e-17}*x21*x42 + V{3.70074341541719e-17}*x29*x39 + V{3.70074341541719e-17}*x32*x41;
auto x96 = -V{7.40148683083438e-17}*x21*x42 + V{7.40148683083438e-17}*x29*x39 + V{7.40148683083438e-17}*x32*x41;
auto x97 = V{3}*x47;
auto x98 = x44 + V{-1};
auto x99 = x46 + x59 - x97 + x98;
auto x100 = -x81*x99;
auto x101 = x71 + x97;
auto x102 = x101*x81;
auto x103 = V{1.11022302462516e-16}*cell[2];
auto x104 = V{1.11022302462516e-16}*cell[4];
auto x105 = V{5.55111512312578e-17}*cell[0];
auto x106 = V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[14] + V{5.55111512312578e-17}*cell[3] + x105 + V{5.55111512312578e-17};
auto x107 = V{5.55111512312578e-17}*cell[10] + V{5.55111512312578e-17}*cell[1] + V{1.11022302462516e-16}*cell[9];
auto x108 = x21*(V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + V{1.11022302462516e-16}*cell[8] + x103 + x104 + x106 + x107);
auto x109 = V{1.11022302462516e-16}*cell[13];
auto x110 = V{5.55111512312578e-17}*cell[11] + V{1.11022302462516e-16}*cell[15] + V{5.55111512312578e-17}*cell[2];
auto x111 = -x29*(V{1.11022302462516e-16}*cell[10] + V{1.11022302462516e-16}*cell[16] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x106 + x109 + x110);
auto x112 = -x32*(V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[17] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[7] + x105 + x107 + x110 + x63 + V{5.55111512312578e-17});
auto x113 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x114 = V{4.5}*(x113*x113);
auto x115 = x114 + x70 + x73 + x84;
auto x116 = -x51;
auto x117 = x100 + x102 + x108 + x111 + x112 + x115*x81 + x116 + V{-2.22044604925031e-16};
auto x118 = V{3}*x43;
auto x119 = -x118 + x49 + x57;
auto x120 = -x119*x81;
auto x121 = x118 + x69 + x73;
auto x122 = x121*x81;
auto x123 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x124 = -V{4.5}*x123*x123;
auto x125 = x124 + x60 + x91;
auto x126 = V{3.33066907387547e-16}*cell[14] + V{3.33066907387547e-16}*cell[5] + x120 + x122 - x125*x81;
auto x127 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{3.33066907387547e-16}*cell[13] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x117 + x126 - x43*x95 - x45*x96 - x47*x90 + x63 + x67 + x75 + x80 + x81*x86 - x81*x94 + x87*x89;
auto x128 = -x123;
auto x129 = x50 + x84;
auto x130 = x129 + x76;
auto x131 = x130 - V{4.5}*x128*x128;
auto x132 = -x92;
auto x133 = x129 + x58;
auto x134 = x133 - V{4.5}*x132*x132;
auto x135 = x81*x89;
auto x136 = x48 + x84 - x88 + x98;
auto x137 = x136*x81;
auto x138 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[17] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[8] + x135 - x137 - x45*x90;
auto x139 = x138 - x62*x81;
auto x140 = x60 - x83 + x84;
auto x141 = -x140*x53;
auto x142 = x53*x86;
auto x143 = -x125*x53;
auto x144 = V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x103 + x117 + x121*x87 + x141 + x142 + x143 - x43*x96 - x47*x95 + x74*x81;
auto x145 = -x114 + x78 + x84;
auto x146 = -x145*x53;
auto x147 = x115*x53;
auto x148 = -x53*x94;
auto x149 = V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[18] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{2.22044604925031e-16}*cell[9] + x104 + x108 + x109 + x111 + x112 + x116 + x126 - x140*x81 + x146 + x147 + x148 - x43*(-V{3.23815048849004e-17}*x21*x42 + V{3.23815048849004e-17}*x29*x39 + V{3.23815048849004e-17}*x32*x41) - x47*x96 - x66*x81 - x87*x99 + V{-2.22044604925031e-16};
auto x150 = -x21*(-x134*x53 + x139 + x149) + x29*(x127 - x53*x62) + x32*(-x131*x53 - x134*x81 + x139 + x144);
auto x151 = V{0.0833333333333333}*cell[11];
auto x152 = V{0.0833333333333333}*cell[12];
auto x153 = V{0.0833333333333333}*cell[2];
auto x154 = V{0.0833333333333333}*cell[3];
auto x155 = V{0.166666666666667}*cell[10];
auto x156 = V{0.166666666666667}*cell[1];
auto x157 = V{0.0833333333333334}*cell[13];
auto x158 = V{0.0833333333333334}*cell[14];
auto x159 = V{0.0833333333333334}*cell[15];
auto x160 = V{0.0833333333333334}*cell[16];
auto x161 = V{0.0833333333333334}*cell[4];
auto x162 = V{0.0833333333333334}*cell[5];
auto x163 = V{0.0833333333333334}*cell[6];
auto x164 = V{0.0833333333333334}*cell[7];
auto x165 = V{3.46944695195361e-18}*cell[0];
auto x166 = V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[2] + x165 + V{3.46944695195361e-18};
auto x167 = V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[14] + V{3.46944695195361e-18}*cell[3];
auto x168 = x29*(V{6.93889390390723e-18}*cell[10] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x166 + x167);
auto x169 = V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[9];
auto x170 = x32*(V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[7] + x166 + x169);
auto x171 = x21*(V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{6.93889390390723e-18}*cell[4] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + V{6.93889390390723e-18}*cell[8] + x165 + x167 + x169 + V{3.46944695195361e-18});
auto x172 = -V{0.0555555555555556}*x21*x42 + V{0.0555555555555556}*x29*x39 + V{0.0555555555555556}*x32*x41;
auto x173 = -V{0.0277777777777778}*x21*x42 + V{0.0277777777777778}*x29*x39 + V{0.0277777777777778}*x32*x41;
auto x174 = x173*x47;
auto x175 = x173*x43;
auto x176 = -x168;
auto x177 = -x170;
auto x178 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x171 + x173*x45 + x176 + x177 + V{-0.0555555555555555};
auto x179 = -x152 - x154 + x157 + x158 + x161 + x162 + x175;
auto x180 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] - x172*x47 + x178 + x179 - x22 - x23 - x24 - x25;
auto x181 = -x151 - x153 + x159 + x160 + x163 + x164 + x174;
auto x182 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] - x172*x43 + x178 + x181 - x33 - x34 - x35 - x36;
auto x183 = -x21*x42 + x29*x39 + x32*x41;
auto x184 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x183;
auto x185 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x184;
auto x186 = -V{0.0138888888888889}*x21*x42 + V{0.0138888888888889}*x29*x39 + V{0.0138888888888889}*x32*x41;
auto x187 = V{0.0138888888888889}*cell[0];
auto x188 = V{0.0138888888888889}*cell[12] + V{0.0277777777777778}*cell[14] + V{0.0138888888888889}*cell[3] + x187 + V{0.0138888888888889};
auto x189 = V{0.0138888888888889}*cell[10] + V{0.0138888888888889}*cell[1] + V{0.0277777777777778}*cell[9];
auto x190 = x21*(V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[4] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + V{0.0277777777777778}*cell[8] + x188 + x189);
auto x191 = V{0.0138888888888889}*cell[11] + V{0.0277777777777778}*cell[15] + V{0.0138888888888889}*cell[2];
auto x192 = -x29*(V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[16] + V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + x188 + x191);
auto x193 = -x32*(V{0.0277777777777778}*cell[12] + V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0277777777777778}*cell[17] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + V{0.0277777777777778}*cell[7] + x187 + x189 + x191 + V{0.0138888888888889});
auto x194 = -V{0.0277777777777778}*x21*x42 + V{0.0277777777777778}*x29*x39 + V{0.0277777777777778}*x32*x41;
auto x195 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x190 + x192 + x193 - x194*x45 + V{0.0138888888888889};
auto x196 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x194*x47;
auto x197 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x186*x43 + x195 + x196;
auto x198 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x185 + x197;
auto x199 = -x53*(x124 + x130);
auto x200 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x185 + x197;
auto x201 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x184;
auto x202 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x194*x43;
auto x203 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x186*x47 + x195 + x202;
auto x204 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x201 + x203;
auto x205 = x133 + x93;
auto x206 = -x205*x53;
auto x207 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x201 + x203;
auto x208 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x183;
auto x209 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x186*x45 + x190 + x192 + x193 + x196 + x202 + V{0.0138888888888889};
auto x210 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x208 + x209;
auto x211 = x61 + x77;
auto x212 = -x211*x53;
auto x213 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x208 + x209;
auto x214 = x138 - x211*x81;
auto x215 = -x21*(x149 + x206 + x214) + x29*(x127 + x212) + x32*(x144 + x199 - x205*x81 + x214);
auto x0 = -x19*(V{0.111111111111111}*x150*x50 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21*(V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[8] + x22 + x23 + x24 + x25 + x27 + x28) - x29*(V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[16] + x27 + x30 + x31) - x32*(V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[7] + x26 + x28 + x30 + x33 + x34 + x35 + x36 + V{0.166666666666667}) - x43*x52 - x45*x52 - x47*x52 + x51 + V{0.833333333333333});
auto x1 = -x19*(V{0.0185185185185185}*x136*x150 + V{0.0555555555555556}) + x20*(x137 + x151 + x152 + x153 + x154 - x155 - x156 - x157 - x158 - x159 - x160 - x161 - x162 - x163 - x164 + x168 + x170 - x171 + x172*x45 - x174 - x175 + x31 + V{0.0555555555555555});
auto x2 = -x19*(V{0.0185185185185185}*x150*x99 + V{0.0555555555555556}) - x20*(x100 + x180);
auto x3 = -x19*(V{0.0185185185185185}*x119*x150 + V{0.0555555555555556}) - x20*(x120 + x182);
auto x4 = -x19*(V{0.00925925925925926}*x140*x150 + V{0.0277777777777778}) - x20*(x141 + x198);
auto x5 = -x19*(V{0.00925925925925926}*x131*x150 + V{0.0277777777777778}) - x20*(x199 + x200);
auto x6 = -x19*(V{0.00925925925925926}*x145*x150 + V{0.0277777777777778}) - x20*(x146 + x204);
auto x7 = -x19*(V{0.00925925925925926}*x134*x150 + V{0.0277777777777778}) - x20*(x206 + x207);
auto x8 = -x19*(V{0.00925925925925926}*x150*x66 + V{0.0277777777777778}) - x20*(x210 + x67);
auto x9 = -x19*(V{0.00925925925925926}*x150*x62 + V{0.0277777777777778}) - x20*(x212 + x213);
auto x10 = x19*(V{0.0185185185185185}*x215*x89 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x135 - x155 - x156 - x171 - x176 - x177 - x179 - x181 + V{0.0555555555555556}*x183*x45 + V{0.0555555555555555});
auto x11 = x19*(V{0.0185185185185185}*x101*x215 + V{-0.0555555555555556}) - x20*(x102 + x180);
auto x12 = x19*(V{0.0185185185185185}*x121*x215 + V{-0.0555555555555556}) - x20*(x122 + x182);
auto x13 = x19*(V{0.00925925925925926}*x215*x86 + V{-0.0277777777777778}) - x20*(x142 + x198);
auto x14 = -x19*(V{0.00925925925925926}*x125*x150 + V{0.0277777777777778}) - x20*(x143 + x200);
auto x15 = x19*(V{0.00925925925925926}*x115*x215 + V{-0.0277777777777778}) - x20*(x147 + x204);
auto x16 = -x19*(V{0.00925925925925926}*x150*x94 + V{0.0277777777777778}) - x20*(x148 + x207);
auto x17 = x19*(V{0.00925925925925926}*x215*x74 + V{-0.0277777777777778}) - x20*(x210 + x75);
auto x18 = -x19*(V{0.00925925925925926}*x150*x79 + V{0.0277777777777778}) - x20*(x213 + x80);
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
return { V{0.333333333333333}*x215, x43 + x45 + x47 };
}
};

}

}

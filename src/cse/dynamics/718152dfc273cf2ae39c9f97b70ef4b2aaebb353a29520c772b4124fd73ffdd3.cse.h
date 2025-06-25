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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<-1, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<-1, 1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1});
auto x22 = V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9];
auto x23 = V{0.166666666666667}*cell[0];
auto x24 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + V{0.333333333333333}*cell[5] + x23 + V{0.166666666666667};
auto x25 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + V{0.333333333333333}*cell[7];
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x27 = V{0.166666666666667}*cell[15];
auto x28 = V{0.166666666666667}*cell[16];
auto x29 = V{0.166666666666667}*cell[6];
auto x30 = V{0.166666666666667}*cell[7];
auto x31 = V{0.166666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.166666666666667}*cell[1];
auto x32 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x33 = V{0.166666666666667}*cell[13];
auto x34 = V{0.166666666666667}*cell[14];
auto x35 = V{0.166666666666667}*cell[4];
auto x36 = V{0.166666666666667}*cell[5];
auto x37 = cell[0] + cell[10] + V{2}*cell[17] + cell[1] + V{1};
auto x38 = cell[12] + cell[3] + V{2}*cell[5];
auto x39 = V{2}*cell[11] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[18] + cell[6] + cell[7] + x37 + x38;
auto x40 = cell[11] + cell[2] + V{2}*cell[7];
auto x41 = V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + cell[4] + cell[5] + V{2}*cell[9] + x37 + x40;
auto x42 = cell[0] + cell[17] + cell[18] + V{2}*cell[1] + V{2}*cell[4] + V{2}*cell[6] + cell[8] + cell[9] + x38 + x40 + V{1};
auto x43 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x44 = V{1.5}*x43;
auto x45 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x46 = V{1.5}*x45;
auto x47 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x48 = V{1.5}*x47;
auto x49 = x46 + x48 + V{-1};
auto x50 = x44 + x49;
auto x51 = x50*(-V{0.111111111111111}*x21*x42 + V{0.111111111111111}*x26*x39 + V{0.111111111111111}*x32*x41);
auto x52 = -V{0.166666666666667}*x21*x42 + V{0.166666666666667}*x26*x39 + V{0.166666666666667}*x32*x41;
auto x53 = -V{0.00925925925925926}*x21*x42 + V{0.00925925925925926}*x26*x39 + V{0.00925925925925926}*x32*x41;
auto x54 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x55 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x54;
auto x56 = -x55;
auto x57 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x58 = -x57;
auto x59 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x60 = x50 + x59;
auto x61 = x58 + x60;
auto x62 = x61 - V{4.5}*x56*x56;
auto x63 = -V{0.0185185185185185}*x21*x42 + V{0.0185185185185185}*x26*x39 + V{0.0185185185185185}*x32*x41;
auto x64 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x65 = -x64;
auto x66 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x67 = -x66;
auto x68 = x60 + x67;
auto x69 = x68 - V{4.5}*x65*x65;
auto x70 = V{1.11022302462516e-16}*cell[12];
auto x71 = V{3}*x43;
auto x72 = x49 + x57 - x71;
auto x73 = -x63*x72;
auto x74 = V{1} - x48;
auto x75 = -x46;
auto x76 = x57 + x75;
auto x77 = x71 + x74 + x76;
auto x78 = x63*x77;
auto x79 = -V{3.70074341541719e-17}*x21*x42 + V{3.70074341541719e-17}*x26*x39 + V{3.70074341541719e-17}*x32*x41;
auto x80 = V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] - x43*x79 + x70 + x73 + x78;
auto x81 = -x63*x69 + x80;
auto x82 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x83 = V{4.5}*(x82*x82);
auto x84 = x57 + x60 - x83;
auto x85 = -x53*x84;
auto x86 = -x44;
auto x87 = x74 + x86;
auto x88 = x59 + x87;
auto x89 = x76 + x83 + x88;
auto x90 = x53*x89;
auto x91 = -x59;
auto x92 = -V{4.5}*x55*x55;
auto x93 = x50 + x57;
auto x94 = x91 + x92 + x93;
auto x95 = -x53*x94;
auto x96 = V{5.55111512312578e-17}*cell[0];
auto x97 = V{5.55111512312578e-17}*cell[12] + V{5.55111512312578e-17}*cell[3] + V{1.11022302462516e-16}*cell[5] + x96 + V{5.55111512312578e-17};
auto x98 = V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[2] + V{1.11022302462516e-16}*cell[7];
auto x99 = x21*(V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{1.11022302462516e-16}*cell[1] + V{1.11022302462516e-16}*cell[4] + V{1.11022302462516e-16}*cell[6] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x97 + x98);
auto x100 = V{1.11022302462516e-16}*cell[11];
auto x101 = V{5.55111512312578e-17}*cell[10] + V{1.11022302462516e-16}*cell[17] + V{5.55111512312578e-17}*cell[1];
auto x102 = -x26*(V{1.11022302462516e-16}*cell[13] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + x100 + x101 + x97);
auto x103 = -x32*(V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[15] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[9] + x101 + x70 + x96 + x98 + V{5.55111512312578e-17});
auto x104 = -x51;
auto x105 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x106 = V{4.5}*(x105*x105);
auto x107 = x66 + x75;
auto x108 = x106 + x107 + x88;
auto x109 = -V{0.037037037037037}*x21*x42 + V{0.037037037037037}*x26*x39 + V{0.037037037037037}*x32*x41;
auto x110 = V{3}*x47;
auto x111 = x107 + x110 + x86 + V{1};
auto x112 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x54;
auto x113 = -V{4.5}*x112*x112;
auto x114 = x113 + x67 + x93;
auto x115 = -V{7.40148683083438e-17}*x21*x42 + V{7.40148683083438e-17}*x26*x39 + V{7.40148683083438e-17}*x32*x41;
auto x116 = V{3}*x45;
auto x117 = x116 + x88;
auto x118 = x117*x63;
auto x119 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x120 = V{4.5}*(x119*x119);
auto x121 = x120 + x66 + x76 + x87;
auto x122 = -V{5.55111512312578e-17}*x21*x42 + V{5.55111512312578e-17}*x26*x39 + V{5.55111512312578e-17}*x32*x41;
auto x123 = x44 + V{-1};
auto x124 = -x116 + x123 + x48 + x59;
auto x125 = x124*x63;
auto x126 = V{1.66533453693773e-16}*cell[10] + V{1.66533453693773e-16}*cell[1] + x118 + x121*x63 - x122*x45 - x125;
auto x127 = V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] + x102 + x103 + x104 + x108*x63 + x109*x111 - x114*x63 - x115*x47 + x126 + x85 + x90 + x95 + x99 + V{-2.22044604925031e-16};
auto x128 = -x112;
auto x129 = x50 + x66;
auto x130 = x129 + x58;
auto x131 = x130 - V{4.5}*x128*x128;
auto x132 = -x110 + x123 + x46 + x66;
auto x133 = -x132*x63;
auto x134 = x111*x63;
auto x135 = x102 + x103 + x104 + x133 + x134 + x99 + V{-2.22044604925031e-16};
auto x136 = x135 - x62*x63;
auto x137 = -x106 + x60 + x66;
auto x138 = -x137*x53;
auto x139 = x108*x53;
auto x140 = -V{4.5}*x64*x64;
auto x141 = x129 + x140 + x91;
auto x142 = -x141*x53;
auto x143 = V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[17] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + x100 + x109*x77 - x115*x43 + x126 + x138 + x139 + x142 - x47*x79 + x63*x89;
auto x144 = -x120 + x66 + x93;
auto x145 = -x144*x53;
auto x146 = x121*x53;
auto x147 = -x114*x53;
auto x148 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] - x109*x124 - x115*x45 - x122*x47 - x137*x63 + x145 + x146 + x147 - x63*x84;
auto x149 = -x21*(-x131*x53 + x136 + x148 + x81) + x26*(x127 - x53*x62 + x81) + x32*(-x131*x63 + x136 + x143 - x53*x69);
auto x150 = V{0.0833333333333333}*cell[11];
auto x151 = V{0.0833333333333333}*cell[12];
auto x152 = V{0.0833333333333333}*cell[2];
auto x153 = V{0.0833333333333333}*cell[3];
auto x154 = V{0.166666666666667}*cell[10];
auto x155 = V{0.166666666666667}*cell[1];
auto x156 = V{0.0833333333333334}*cell[13];
auto x157 = V{0.0833333333333334}*cell[14];
auto x158 = V{0.0833333333333334}*cell[15];
auto x159 = V{0.0833333333333334}*cell[16];
auto x160 = V{0.0833333333333334}*cell[4];
auto x161 = V{0.0833333333333334}*cell[5];
auto x162 = V{0.0833333333333334}*cell[6];
auto x163 = V{0.0833333333333334}*cell[7];
auto x164 = V{3.46944695195361e-18}*cell[0];
auto x165 = V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[1] + x164 + V{3.46944695195361e-18};
auto x166 = V{3.46944695195361e-18}*cell[12] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[5];
auto x167 = x26*(V{6.93889390390723e-18}*cell[11] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x165 + x166);
auto x168 = V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[7];
auto x169 = x32*(V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[9] + x165 + x168);
auto x170 = x21*(V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[6] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x164 + x166 + x168 + V{3.46944695195361e-18});
auto x171 = -V{0.0555555555555556}*x21*x42 + V{0.0555555555555556}*x26*x39 + V{0.0555555555555556}*x32*x41;
auto x172 = -V{0.0277777777777778}*x21*x42 + V{0.0277777777777778}*x26*x39 + V{0.0277777777777778}*x32*x41;
auto x173 = x172*x47;
auto x174 = x172*x43;
auto x175 = -x167;
auto x176 = -x169;
auto x177 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x170 + x172*x45 + x175 + x176 + V{-0.0555555555555555};
auto x178 = -x151 - x153 + x156 + x157 + x160 + x161 + x174;
auto x179 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] - x171*x47 + x177 + x178 - x27 - x28 - x29 - x30;
auto x180 = -x150 - x152 + x158 + x159 + x162 + x163 + x173;
auto x181 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] - x171*x43 + x177 + x180 - x33 - x34 - x35 - x36;
auto x182 = -x21*x42 + x26*x39 + x32*x41;
auto x183 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x182;
auto x184 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x183;
auto x185 = -V{0.0138888888888889}*x21*x42 + V{0.0138888888888889}*x26*x39 + V{0.0138888888888889}*x32*x41;
auto x186 = V{0.0138888888888889}*cell[0];
auto x187 = V{0.0138888888888889}*cell[12] + V{0.0138888888888889}*cell[3] + V{0.0277777777777778}*cell[5] + x186 + V{0.0138888888888889};
auto x188 = V{0.0138888888888889}*cell[11] + V{0.0138888888888889}*cell[2] + V{0.0277777777777778}*cell[7];
auto x189 = x21*(V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[6] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + x187 + x188);
auto x190 = V{0.0138888888888889}*cell[10] + V{0.0277777777777778}*cell[17] + V{0.0138888888888889}*cell[1];
auto x191 = -x26*(V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[13] + V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0277777777777778}*cell[18] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + x187 + x190);
auto x192 = -x32*(V{0.0277777777777778}*cell[12] + V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + V{0.0277777777777778}*cell[9] + x186 + x188 + x190 + V{0.0138888888888889});
auto x193 = -V{0.0277777777777778}*x21*x42 + V{0.0277777777777778}*x26*x39 + V{0.0277777777777778}*x32*x41;
auto x194 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x189 + x191 + x192 - x193*x45 + V{0.0138888888888889};
auto x195 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x193*x47;
auto x196 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x185*x43 + x194 + x195;
auto x197 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x184 + x196;
auto x198 = x140 + x68;
auto x199 = -x198*x53;
auto x200 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x184 + x196;
auto x201 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x183;
auto x202 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x193*x43;
auto x203 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x185*x47 + x194 + x202;
auto x204 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x201 + x203;
auto x205 = x61 + x92;
auto x206 = -x205*x53;
auto x207 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x201 + x203;
auto x208 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x182;
auto x209 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x185*x45 + x189 + x191 + x192 + x195 + x202 + V{0.0138888888888889};
auto x210 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x208 + x209;
auto x211 = x113 + x130;
auto x212 = -x211*x53;
auto x213 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x208 + x209;
auto x214 = -x198*x63 + x80;
auto x215 = x135 - x205*x63;
auto x216 = -x21*(x148 + x212 + x214 + x215) + x26*(x127 + x206 + x214) + x32*(x143 + x199 - x211*x63 + x215);
auto x0 = -x19*(V{0.111111111111111}*x149*x50 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21*(V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[6] + x22 + x24 + x25) - x26*(V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[18] + x24 + x27 + x28 + x29 + x30 + x31) - x32*(V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[9] + x23 + x25 + x31 + x33 + x34 + x35 + x36 + V{0.166666666666667}) - x43*x52 - x45*x52 - x47*x52 + x51 + V{0.833333333333333});
auto x1 = -x19*(V{0.0185185185185185}*x124*x149 + V{0.0555555555555556}) + x20*(x125 + x150 + x151 + x152 + x153 - x154 - x155 - x156 - x157 - x158 - x159 - x160 - x161 - x162 - x163 + x167 + x169 - x170 + x171*x45 - x173 - x174 + x22 + V{0.0555555555555555});
auto x2 = -x19*(V{0.0185185185185185}*x132*x149 + V{0.0555555555555556}) - x20*(x133 + x179);
auto x3 = -x19*(V{0.0185185185185185}*x149*x72 + V{0.0555555555555556}) - x20*(x181 + x73);
auto x4 = -x19*(V{0.00925925925925926}*x137*x149 + V{0.0277777777777778}) - x20*(x138 + x197);
auto x5 = -x19*(V{0.00925925925925926}*x149*x69 + V{0.0277777777777778}) - x20*(x199 + x200);
auto x6 = -x19*(V{0.00925925925925926}*x149*x84 + V{0.0277777777777778}) - x20*(x204 + x85);
auto x7 = -x19*(V{0.00925925925925926}*x149*x62 + V{0.0277777777777778}) - x20*(x206 + x207);
auto x8 = -x19*(V{0.00925925925925926}*x144*x149 + V{0.0277777777777778}) - x20*(x145 + x210);
auto x9 = -x19*(V{0.00925925925925926}*x131*x149 + V{0.0277777777777778}) - x20*(x212 + x213);
auto x10 = x19*(V{0.0185185185185185}*x117*x216 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x118 - x154 - x155 - x170 - x175 - x176 - x178 - x180 + V{0.0555555555555556}*x182*x45 + V{0.0555555555555555});
auto x11 = x19*(V{0.0185185185185185}*x111*x216 + V{-0.0555555555555556}) - x20*(x134 + x179);
auto x12 = x19*(V{0.0185185185185185}*x216*x77 + V{-0.0555555555555556}) - x20*(x181 + x78);
auto x13 = x19*(V{0.00925925925925926}*x108*x216 + V{-0.0277777777777778}) - x20*(x139 + x197);
auto x14 = -x19*(V{0.00925925925925926}*x141*x149 + V{0.0277777777777778}) - x20*(x142 + x200);
auto x15 = x19*(V{0.00925925925925926}*x216*x89 + V{-0.0277777777777778}) - x20*(x204 + x90);
auto x16 = -x19*(V{0.00925925925925926}*x149*x94 + V{0.0277777777777778}) - x20*(x207 + x95);
auto x17 = x19*(V{0.00925925925925926}*x121*x216 + V{-0.0277777777777778}) - x20*(x146 + x210);
auto x18 = -x19*(V{0.00925925925925926}*x114*x149 + V{0.0277777777777778}) - x20*(x147 + x213);
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
return { V{0.333333333333333}*x216, x43 + x45 + x47 };
}
};

}

}

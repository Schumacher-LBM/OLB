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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<1, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<1, 1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{0.5}*cell[12];
auto x22 = V{0.5}*cell[1];
auto x23 = V{2}*cell[7];
auto x24 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1});
auto x25 = V{0.25}*cell[0] + V{0.25}*cell[11] + V{0.25}*cell[2] + V{0.5}*cell[7] + V{0.25};
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x27 = cell[0] + cell[11] + cell[2] + x23 + V{1};
auto x28 = cell[10] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + V{2}*cell[17] + cell[1] + cell[4] + cell[5] + V{2}*cell[9] + x27;
auto x29 = cell[12] + cell[17] + cell[18] + V{2}*cell[1] + cell[3] + V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + cell[8] + cell[9] + x27;
auto x30 = -V{0.0138888888888889}*x24*x29 + V{0.0138888888888889}*x26*x28;
auto x31 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x31;
auto x33 = -x32;
auto x34 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x35 = -x34;
auto x36 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x37 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x38 = V{1.5}*x37;
auto x39 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x40 = V{1.5}*x39;
auto x41 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x42 = V{1.5}*x41;
auto x43 = x40 + x42 + V{-1};
auto x44 = x38 + x43;
auto x45 = x36 + x44;
auto x46 = x35 + x45;
auto x47 = x46 - V{4.5}*x33*x33;
auto x48 = x44*(-V{0.166666666666667}*x24*x29 + V{0.166666666666667}*x26*x28);
auto x49 = -V{0.25}*x24*x29 + V{0.25}*x26*x28;
auto x50 = -V{4.5}*x32*x32;
auto x51 = -x36 + x44;
auto x52 = x34 + x50 + x51;
auto x53 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x54 = -x53;
auto x55 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x56 = -x55;
auto x57 = x45 + x56;
auto x58 = x57 - V{4.5}*x54*x54;
auto x59 = -V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28;
auto x60 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x31;
auto x61 = -x60;
auto x62 = x44 + x55;
auto x63 = x35 + x62;
auto x64 = x63 - V{4.5}*x61*x61;
auto x65 = -V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28;
auto x66 = V{3}*x39;
auto x67 = -x38;
auto x68 = V{1} - x42;
auto x69 = x67 + x68;
auto x70 = x36 + x69;
auto x71 = x66 + x70;
auto x72 = x59*x71;
auto x73 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x74 = V{4.5}*(x73*x73);
auto x75 = x45 + x55 - x74;
auto x76 = -x30*x75;
auto x77 = -x40;
auto x78 = x55 + x77;
auto x79 = x70 + x74 + x78;
auto x80 = x30*x79;
auto x81 = -V{4.5}*x53*x53;
auto x82 = x51 + x55 + x81;
auto x83 = -x30*x82;
auto x84 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x85 = V{4.5}*(x84*x84);
auto x86 = x34 + x77;
auto x87 = x70 + x85 + x86;
auto x88 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x89 = V{4.5}*(x88*x88);
auto x90 = x34 + x69 + x78 + x89;
auto x91 = -V{0.0555555555555556}*x24*x29 + V{0.0555555555555556}*x26*x28;
auto x92 = V{3}*x41;
auto x93 = x67 + x86 + x92 + V{1};
auto x94 = -V{5.55111512312578e-17}*x24*x29 + V{5.55111512312578e-17}*x26*x28;
auto x95 = -V{1.11022302462516e-16}*x24*x29 + V{1.11022302462516e-16}*x26*x28;
auto x96 = -V{8.32667268468867e-17}*x24*x29 + V{8.32667268468867e-17}*x26*x28;
auto x97 = x38 + V{-1};
auto x98 = x36 + x42 - x66 + x97;
auto x99 = x59*x98;
auto x100 = V{3}*x37;
auto x101 = -x100 + x43 + x55;
auto x102 = -x101*x59;
auto x103 = x100 + x68 + x78;
auto x104 = x103*x59;
auto x105 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[2] + V{1.11022302462516e-16}*cell[7] + V{5.55111512312578e-17};
auto x106 = V{1.11022302462516e-16}*cell[12];
auto x107 = x102 + x104 + x24*(V{5.55111512312578e-17}*cell[12] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{1.11022302462516e-16}*cell[1] + V{5.55111512312578e-17}*cell[3] + V{1.11022302462516e-16}*cell[4] + V{1.11022302462516e-16}*cell[5] + V{1.11022302462516e-16}*cell[6] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x105) - x26*(V{5.55111512312578e-17}*cell[10] + V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[15] + V{1.11022302462516e-16}*cell[17] + V{5.55111512312578e-17}*cell[1] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[9] + x105 + x106) - x48 + V{-2.22044604925031e-16};
auto x108 = V{1.66533453693773e-16}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{2.22044604925031e-16}*cell[17] + V{1.66533453693773e-16}*cell[1] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{6.66133814775094e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + x107 - x37*x94 - x39*x96 - x41*x95 - x52*(-V{4.62592926927149e-18}*x24*x29 + V{4.62592926927149e-18}*x26*x28) + x59*x87 + x59*x90 + x72 + x76 + x80 + x83 + x91*x93 - x99;
auto x109 = -V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28;
auto x110 = x34 + x40 - x92 + x97;
auto x111 = -x110*x59;
auto x112 = x59*x93;
auto x113 = x34 + x62 - x89;
auto x114 = -x113*x30;
auto x115 = x30*x90;
auto x116 = -V{4.5}*x60*x60;
auto x117 = x116 + x34 + x44 + x56;
auto x118 = -x117*x30;
auto x119 = x34 + x45 - x85;
auto x120 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.44089209850063e-16}*cell[15] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{4.44089209850063e-16}*cell[6] + V{4.44089209850063e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x106 + x107 + x111 + x112 + x114 + x115 + x118 - x119*x59 - x37*x96 - x39*x95 - x41*x94 - x52*(-V{3.08395284618099e-18}*x24*x29 + V{3.08395284618099e-18}*x26*x28) - x59*x75 - x91*x98;
auto x121 = -x24*(-x109*x47 + x120 - x30*x64 - x58*x59) + x26*(x108 - x30*x58 - x47*x65 - x59*x64);
auto x122 = V{0.0833333333333333}*cell[11];
auto x123 = V{0.0833333333333333}*cell[12];
auto x124 = V{0.0833333333333333}*cell[2];
auto x125 = V{0.0833333333333333}*cell[3];
auto x126 = V{0.166666666666667}*cell[10];
auto x127 = V{0.166666666666667}*cell[1];
auto x128 = V{0.0833333333333334}*cell[13];
auto x129 = V{0.0833333333333334}*cell[14];
auto x130 = V{0.0833333333333334}*cell[15];
auto x131 = V{0.0833333333333334}*cell[4];
auto x132 = V{0.0833333333333334}*cell[5];
auto x133 = V{0.0833333333333334}*cell[6];
auto x134 = V{0.166666666666667}*cell[7];
auto x135 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[7] + V{3.46944695195361e-18};
auto x136 = x26*(V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[1] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[9] + x135);
auto x137 = x24*(V{3.46944695195361e-18}*cell[12] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] + V{6.93889390390723e-18}*cell[6] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x135);
auto x138 = -V{0.0833333333333333}*x24*x29 + V{0.0833333333333333}*x26*x28;
auto x139 = -V{0.00115740740740741}*x24*x29 + V{0.00115740740740741}*x26*x28;
auto x140 = x139*x52;
auto x141 = x139*x47;
auto x142 = -V{0.0416666666666667}*x24*x29 + V{0.0416666666666667}*x26*x28;
auto x143 = x142*x37;
auto x144 = x142*x41;
auto x145 = -V{0.00231481481481481}*x24*x29 + V{0.00231481481481481}*x26*x28;
auto x146 = x46 + x50;
auto x147 = -x136;
auto x148 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x137 + x142*x39 + x147 + V{-0.0555555555555555};
auto x149 = -x123 - x125 + x128 + x129 + x131 + x132 + x144;
auto x150 = V{0.166666666666667}*cell[11] - V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[2] - V{0.166666666666667}*cell[6] - V{0.333333333333333}*cell[7] - x138*x37 - x145*x146 + x145*x52 + x148 + x149;
auto x151 = -x122 - x124 + x130 + x133 + x134 - x140 + x143;
auto x152 = V{0.166666666666667}*cell[12] - V{0.166666666666667}*cell[13] - V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[4] - V{0.166666666666667}*cell[5] - x138*x41 + x139*x146 + x148 + x151;
auto x153 = -x24*x29 + x26*x28;
auto x154 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x153;
auto x155 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x154;
auto x156 = -V{0.0208333333333333}*x24*x29 + V{0.0208333333333333}*x26*x28;
auto x157 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[11] + V{0.0208333333333333}*cell[2] + V{0.0416666666666667}*cell[7] + V{0.0208333333333333};
auto x158 = x24*(V{0.0208333333333333}*cell[12] + V{0.0208333333333333}*cell[17] + V{0.0208333333333333}*cell[18] + V{0.0416666666666667}*cell[1] + V{0.0208333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + V{0.0416666666666667}*cell[6] + V{0.0208333333333333}*cell[8] + V{0.0208333333333333}*cell[9] + x157);
auto x159 = -x26*(V{0.0208333333333333}*cell[10] + V{0.0416666666666667}*cell[12] + V{0.0208333333333333}*cell[13] + V{0.0208333333333333}*cell[14] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[17] + V{0.0208333333333333}*cell[1] + V{0.0208333333333333}*cell[4] + V{0.0208333333333333}*cell[5] + V{0.0416666666666667}*cell[9] + x157);
auto x160 = -V{0.0416666666666667}*x24*x29 + V{0.0416666666666667}*x26*x28;
auto x161 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x158 + x159 - x160*x39 + V{0.0138888888888889};
auto x162 = -V{0.000578703703703704}*x24*x29 + V{0.000578703703703704}*x26*x28;
auto x163 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0833333333333334}*cell[7] + x146*x162 - x160*x37 - x162*x52;
auto x164 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x156*x41 + x161 + x163;
auto x165 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x155 + x164;
auto x166 = x57 + x81;
auto x167 = -x166*x30;
auto x168 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x155 + x164;
auto x169 = -V{0.00115740740740741}*x24*x29 + V{0.00115740740740741}*x26*x28;
auto x170 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x154;
auto x171 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x160*x41;
auto x172 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x156*x37 + x161 + x171;
auto x173 = V{0.416666666666667}*cell[15] + V{0.416666666666667}*cell[6] - V{0.166666666666667}*cell[7] - x146*x169 + x169*x52 - x170 + x172;
auto x174 = -V{0.00578703703703704}*x24*x29 + V{0.00578703703703704}*x26*x28;
auto x175 = -V{0.0833333333333333}*cell[15] - V{0.0833333333333333}*cell[6] + V{0.833333333333333}*cell[7] + x170 + x172;
auto x176 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x153;
auto x177 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x156*x39 + x158 + x159 + x163 + x171 + V{0.0138888888888889};
auto x178 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x176 + x177;
auto x179 = x116 + x63;
auto x180 = -x179*x30;
auto x181 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x176 + x177;
auto x182 = -x24*(-x109*x146 + x120 - x166*x59 + x180) + x26*(x108 - x146*x65 + x167 - x179*x59);
auto x0 = -x19*(V{0.166666666666667}*x121*x44 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[8] + V{1}*cell[9] + x21 + x22 + x23 + x24*(V{0.25}*cell[12] + V{0.25}*cell[17] + V{0.25}*cell[18] + V{0.25}*cell[3] + V{0.5}*cell[4] + V{0.5}*cell[5] + V{0.5}*cell[6] + V{0.25}*cell[8] + V{0.25}*cell[9] + x22 + x25) - x26*(V{0.25}*cell[10] + V{0.25}*cell[13] + V{0.25}*cell[14] + V{0.5}*cell[15] + V{0.5}*cell[17] + V{0.25}*cell[1] + V{0.25}*cell[4] + V{0.25}*cell[5] + V{0.5}*cell[9] + x21 + x25) + x30*x47 - x30*x52 - x37*x49 - x39*x49 - x41*x49 + x48 + V{0.833333333333333});
auto x1 = -x19*(V{0.0277777777777778}*x121*x98 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x122 + x123 + x124 + x125 - x126 - x127 - x128 - x129 - x130 - x131 - x132 - x133 - x134 + x136 - x137 + x138*x39 + x140 - x141 - x143 - x144 + x99 + V{0.0555555555555555});
auto x2 = -x19*(V{0.0277777777777778}*x101*x121 + V{0.0555555555555556}) - x20*(x102 + x150);
auto x3 = -x19*(V{0.0277777777777778}*x110*x121 + V{0.0555555555555556}) - x20*(x111 + x152);
auto x4 = -x19*(V{0.0138888888888889}*x121*x75 + V{0.0277777777777778}) - x20*(x165 + x76);
auto x5 = -x19*(V{0.0138888888888889}*x121*x58 + V{0.0277777777777778}) - x20*(x167 + x168);
auto x6 = -x19*(V{0.0138888888888889}*x119*x121 + V{0.0277777777777778}) - x20*(-x119*x30 + x173);
auto x7 = -x19*(V{0.0138888888888889}*x121*x47 + V{0.0277777777777778}) - x20*(-x146*(-V{0.00810185185185185}*x24*x29 + V{0.00810185185185185}*x26*x28) - x174*x52 + x175);
auto x8 = -x19*(V{0.0138888888888889}*x113*x121 + V{0.0277777777777778}) - x20*(x114 + x178);
auto x9 = -x19*(V{0.0138888888888889}*x121*x64 + V{0.0277777777777778}) - x20*(x180 + x181);
auto x10 = x19*(V{0.0277777777777778}*x182*x71 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x126 - x127 - x137 - x141 - x147 - x149 - x151 + V{0.0833333333333333}*x153*x39 - x72 + V{0.0555555555555555});
auto x11 = x19*(V{0.0277777777777778}*x103*x182 + V{-0.0555555555555556}) - x20*(x104 + x150);
auto x12 = x19*(V{0.0277777777777778}*x182*x93 + V{-0.0555555555555556}) - x20*(x112 + x152);
auto x13 = x19*(V{0.0138888888888889}*x182*x79 + V{-0.0277777777777778}) - x20*(x165 + x80);
auto x14 = -x19*(V{0.0138888888888889}*x121*x82 + V{0.0277777777777778}) - x20*(x168 + x83);
auto x15 = x19*(V{0.0138888888888889}*x182*x87 + V{-0.0277777777777778}) - x20*(x173 + x30*x87);
auto x16 = -x19*(V{0.0138888888888889}*x121*x52 + V{0.0277777777777778}) - x20*(x146*x174 + x175 - x52*(-V{0.0196759259259259}*x24*x29 + V{0.0196759259259259}*x26*x28));
auto x17 = x19*(V{0.0138888888888889}*x182*x90 + V{-0.0277777777777778}) - x20*(x115 + x178);
auto x18 = -x19*(V{0.0138888888888889}*x117*x121 + V{0.0277777777777778}) - x20*(x118 + x181);
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
return { V{0.5}*x182, x37 + x39 + x41 };
}
};

}

}

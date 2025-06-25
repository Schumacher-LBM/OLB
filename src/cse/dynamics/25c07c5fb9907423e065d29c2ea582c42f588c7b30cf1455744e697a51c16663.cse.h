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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, 1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{0.5}*cell[11];
auto x22 = V{0.5}*cell[3];
auto x23 = V{2}*cell[18];
auto x24 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1});
auto x25 = V{0.25}*cell[0] + V{0.25}*cell[10] + V{0.5}*cell[18] + V{0.25}*cell[1] + V{0.25};
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x27 = cell[0] + cell[10] + cell[1] + x23 + V{1};
auto x28 = V{2}*cell[11] + cell[12] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + cell[3] + V{2}*cell[5] + cell[6] + cell[7] + x27;
auto x29 = cell[11] + cell[13] + cell[14] + V{2}*cell[16] + cell[2] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + V{2}*cell[8] + x27;
auto x30 = -V{0.0138888888888889}*x24*x29 + V{0.0138888888888889}*x26*x28;
auto x31 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x32 = -x31;
auto x33 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x34 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x33;
auto x35 = -V{4.5}*x34*x34;
auto x36 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x37 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x38 = V{1.5}*x37;
auto x39 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x40 = V{1.5}*x39;
auto x41 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x42 = V{1.5}*x41;
auto x43 = x40 + x42 + V{-1};
auto x44 = x38 + x43;
auto x45 = x36 + x44;
auto x46 = x32 + x35 + x45;
auto x47 = x44*(-V{0.166666666666667}*x24*x29 + V{0.166666666666667}*x26*x28);
auto x48 = -V{0.25}*x24*x29 + V{0.25}*x26*x28;
auto x49 = -x34;
auto x50 = -x36 + x44;
auto x51 = x31 + x50;
auto x52 = x51 - V{4.5}*x49*x49;
auto x53 = V{1.11022302462516e-16}*cell[3];
auto x54 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x55 = V{4.5}*(x54*x54);
auto x56 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x57 = -x38;
auto x58 = V{1} - x42;
auto x59 = x57 + x58;
auto x60 = x56 + x59;
auto x61 = -x40;
auto x62 = x36 + x61;
auto x63 = x55 + x60 + x62;
auto x64 = -V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28;
auto x65 = V{3}*x37;
auto x66 = x58 + x62 + x65;
auto x67 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x68 = V{4.5}*(x67*x67);
auto x69 = x31 + x61;
auto x70 = x60 + x68 + x69;
auto x71 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x72 = V{4.5}*(x71*x71);
auto x73 = x31 + x59 + x62 + x72;
auto x74 = -V{0.0555555555555556}*x24*x29 + V{0.0555555555555556}*x26*x28;
auto x75 = V{3}*x41;
auto x76 = x57 + x69 + x75 + V{1};
auto x77 = -V{5.55111512312578e-17}*x24*x29 + V{5.55111512312578e-17}*x26*x28;
auto x78 = -V{1.11022302462516e-16}*x24*x29 + V{1.11022302462516e-16}*x26*x28;
auto x79 = x45 - x55 + x56;
auto x80 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x33;
auto x81 = -x80;
auto x82 = x50 + x56;
auto x83 = x82 - V{4.5}*x81*x81;
auto x84 = -x56;
auto x85 = -V{4.5}*x80*x80;
auto x86 = x45 + x84 + x85;
auto x87 = x36 + x43 - x65;
auto x88 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x89 = -x88;
auto x90 = x44 + x56;
auto x91 = x32 + x90;
auto x92 = x91 - V{4.5}*x89*x89;
auto x93 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[10] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[1] + V{5.55111512312578e-17};
auto x94 = V{1.11022302462516e-16}*cell[11];
auto x95 = V{3}*x39;
auto x96 = x60 + x95;
auto x97 = x38 + V{-1};
auto x98 = x42 + x56 - x95 + x97;
auto x99 = V{1.66533453693773e-16}*cell[10] + V{4.44089209850063e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + x24*(V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[16] + V{5.55111512312578e-17}*cell[2] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[6] + V{1.11022302462516e-16}*cell[8] + x53 + x93) - x26*(V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{1.11022302462516e-16}*cell[17] + V{5.55111512312578e-17}*cell[3] + V{1.11022302462516e-16}*cell[5] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + x93 + x94) - x39*(-V{8.32667268468867e-17}*x24*x29 + V{8.32667268468867e-17}*x26*x28) - x46*(-V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28) - x47 - x52*(-V{3.08395284618099e-18}*x24*x29 + V{3.08395284618099e-18}*x26*x28) + x64*x96 - x64*x98 + V{-2.22044604925031e-16};
auto x100 = x26*(V{2.22044604925031e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[2] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + x30*x63 - x30*x79 - x30*x83 - x30*x86 - x37*x77 - x41*x78 + x53 + x64*x66 + x64*x70 + x64*x73 - x64*x87 - x64*x92 + x74*x76 + x99);
auto x101 = x31 - x68 + x90;
auto x102 = -V{4.5}*x88*x88;
auto x103 = x102 + x31 + x44 + x84;
auto x104 = x31 + x40 - x75 + x97;
auto x105 = x31 + x45 - x72;
auto x106 = x24*(V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] - x101*x30 - x103*x30 - x104*x64 - x105*x64 + x30*x70 - x30*x92 - x37*x78 - x41*x77 + x64*x76 - x64*x79 - x64*x86 - x74*x87 + x94 + x99);
auto x107 = x100 - x106;
auto x108 = V{0.0277777777777778}*x24*x29 - V{0.0277777777777778}*x26*x28;
auto x109 = V{0.00231481481481481}*x24*x29 - V{0.00231481481481481}*x26*x28;
auto x110 = V{0.0833333333333333}*x24*x29 - V{0.0833333333333333}*x26*x28;
auto x111 = x35 + x51;
auto x112 = V{0.0833333333333333}*cell[12];
auto x113 = V{0.0833333333333333}*cell[3];
auto x114 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[1] + V{3.46944695195361e-18};
auto x115 = V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[3] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[8] + x114;
auto x116 = x115*x24;
auto x117 = x26*(V{6.93889390390723e-18}*cell[11] + V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[5] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x114);
auto x118 = -x117;
auto x119 = V{0.0416666666666667}*x24*x29 - V{0.0416666666666667}*x26*x28;
auto x120 = x119*x37;
auto x121 = V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] - x112 - x113 + x116 + x118 - x120 + V{-0.0555555555555555};
auto x122 = V{0.0833333333333333}*cell[11];
auto x123 = V{0.0833333333333333}*cell[2];
auto x124 = x119*x41;
auto x125 = V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] - x122 - x123 - x124;
auto x126 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] - V{0.333333333333333}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - x109*x111 + x109*x46 + x110*x39 + x121 + x125;
auto x127 = V{0.166666666666667}*cell[15];
auto x128 = V{0.166666666666667}*cell[16];
auto x129 = V{0.166666666666667}*cell[6];
auto x130 = V{0.166666666666667}*cell[7];
auto x131 = V{0.0833333333333333}*cell[10];
auto x132 = V{0.0833333333333333}*cell[1];
auto x133 = V{0.00115740740740741}*x24*x29 - V{0.00115740740740741}*x26*x28;
auto x134 = x133*x46;
auto x135 = x119*x39;
auto x136 = V{0.0833333333333334}*cell[17] + V{0.166666666666667}*cell[18] + V{0.0833333333333334}*cell[8] + x111*x133 - x131 - x132 - x134 - x135;
auto x137 = V{0.166666666666667}*cell[13];
auto x138 = V{0.166666666666667}*cell[14];
auto x139 = V{0.166666666666667}*cell[4];
auto x140 = V{0.166666666666667}*cell[5];
auto x141 = V{0.0138888888888889}*x24*x29 - V{0.0138888888888889}*x26*x28;
auto x142 = x24*x29 - x26*x28;
auto x143 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x142;
auto x144 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x143;
auto x145 = V{0.0208333333333333}*x24*x29 - V{0.0208333333333333}*x26*x28;
auto x146 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[10] + V{0.0416666666666667}*cell[18] + V{0.0208333333333333}*cell[1] + V{0.0208333333333333};
auto x147 = x24*(V{0.0208333333333333}*cell[11] + V{0.0208333333333333}*cell[13] + V{0.0208333333333333}*cell[14] + V{0.0416666666666667}*cell[16] + V{0.0208333333333333}*cell[2] + V{0.0416666666666667}*cell[3] + V{0.0208333333333333}*cell[4] + V{0.0208333333333333}*cell[5] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[8] + x146);
auto x148 = -x26*(V{0.0416666666666667}*cell[11] + V{0.0208333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0208333333333333}*cell[15] + V{0.0208333333333333}*cell[16] + V{0.0416666666666667}*cell[17] + V{0.0208333333333333}*cell[3] + V{0.0416666666666667}*cell[5] + V{0.0208333333333333}*cell[6] + V{0.0208333333333333}*cell[7] + x146);
auto x149 = V{0.0416666666666667}*x24*x29 - V{0.0416666666666667}*x26*x28;
auto x150 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x147 + x148 + x149*x41 + V{0.0138888888888889};
auto x151 = V{0.000578703703703704}*x24*x29 - V{0.000578703703703704}*x26*x28;
auto x152 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + x111*x151 + x149*x39 - x151*x46;
auto x153 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x145*x37 + x150 + x152;
auto x154 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x144 + x153;
auto x155 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x144 + x153;
auto x156 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x143;
auto x157 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + x149*x37;
auto x158 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x145*x41 + x147 + x148 + x152 + x157 + V{0.0138888888888889};
auto x159 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x156 + x158;
auto x160 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x156 + x158;
auto x161 = V{0.00115740740740741}*x24*x29 - V{0.00115740740740741}*x26*x28;
auto x162 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x142;
auto x163 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x145*x39 + x150 + x157;
auto x164 = V{0.416666666666667}*cell[17] - V{0.166666666666667}*cell[18] + V{0.416666666666667}*cell[8] - x111*x161 + x161*x46 + x162 + x163;
auto x165 = V{0.00578703703703704}*x24*x29 - V{0.00578703703703704}*x26*x28;
auto x166 = -V{0.0833333333333333}*cell[17] + V{0.833333333333333}*cell[18] - V{0.0833333333333333}*cell[8] - x162 + x163;
auto x167 = -x100 + x106;
auto x168 = -V{0.0833333333333334}*cell[17] - V{0.166666666666667}*cell[18] - V{0.0833333333333334}*cell[8] - V{0.00115740740740741}*x111*x142 - x115*x24 + x117 + x131 + x132 + x134 + x135 + V{0.0555555555555555};
auto x0 = -x19*(V{0.166666666666667}*x107*x44 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + x21 + x22 + x23 + x24*(V{0.25}*cell[11] + V{0.25}*cell[13] + V{0.25}*cell[14] + V{0.5}*cell[16] + V{0.25}*cell[2] + V{0.25}*cell[4] + V{0.25}*cell[5] + V{0.5}*cell[6] + V{0.5}*cell[8] + x22 + x25) - x26*(V{0.25}*cell[12] + V{0.5}*cell[13] + V{0.25}*cell[15] + V{0.25}*cell[16] + V{0.5}*cell[17] + V{0.25}*cell[3] + V{0.5}*cell[5] + V{0.25}*cell[6] + V{0.25}*cell[7] + x21 + x25) + x30*x46 - x30*x52 - x37*x48 - x39*x48 - x41*x48 + x47 + V{0.833333333333333});
auto x1 = -x19*(V{0.0277777777777778}*x107*x98 + V{0.0555555555555556}) - x20*(x108*x98 + x126);
auto x2 = -x19*(V{0.0277777777777778}*x104*x107 + V{0.0555555555555556}) - x20*(V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + x104*x108 + x110*x41 + x121 - x127 - x128 - x129 - x130 + x136);
auto x3 = -x19*(V{0.0277777777777778}*x107*x87 + V{0.0555555555555556}) - x20*(V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + x108*x87 + x110*x37 + x116 + x118 + x125 + x136 - x137 - x138 - x139 - x140 + V{-0.0555555555555555});
auto x4 = -x19*(V{0.0138888888888889}*x101*x107 + V{0.0277777777777778}) - x20*(x101*x141 + x154);
auto x5 = -x19*(V{0.0138888888888889}*x107*x92 + V{0.0277777777777778}) - x20*(x141*(x102 + x91) + x155);
auto x6 = -x19*(V{0.0138888888888889}*x107*x79 + V{0.0277777777777778}) - x20*(x141*x79 + x159);
auto x7 = -x19*(V{0.0138888888888889}*x107*x83 + V{0.0277777777777778}) - x20*(x141*(x82 + x85) + x160);
auto x8 = -x19*(V{0.0138888888888889}*x105*x107 + V{0.0277777777777778}) - x20*(x105*x141 + x164);
auto x9 = -x19*(V{0.0138888888888889}*x107*x52 + V{0.0277777777777778}) - x20*(x111*(V{0.0196759259259259}*x24*x29 - V{0.0196759259259259}*x26*x28) - x165*x46 + x166);
auto x10 = -x19*(V{0.0277777777777778}*x167*x96 + V{0.0555555555555556}) - x20*(-x108*x96 + x126);
auto x11 = -x19*(V{0.0277777777777778}*x167*x76 + V{0.0555555555555556}) - x20*(V{0.166666666666667}*cell[11] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] + V{0.166666666666667}*cell[2] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] - x108*x76 - x112 - x113 - x120 - x127 - x128 - x129 - x130 + V{0.0833333333333333}*x142*x41 - x168);
auto x12 = -x19*(V{0.0277777777777778}*x167*x66 + V{0.0555555555555556}) - x20*(V{0.166666666666667}*cell[12] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] + V{0.166666666666667}*cell[3] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] - x108*x66 - x122 - x123 - x124 - x137 - x138 - x139 - x140 + V{0.0833333333333333}*x142*x37 - x168);
auto x13 = -x19*(V{0.0138888888888889}*x167*x70 + V{0.0277777777777778}) - x20*(-x141*x70 + x154);
auto x14 = -x19*(V{0.0138888888888889}*x103*x107 + V{0.0277777777777778}) - x20*(x103*x141 + x155);
auto x15 = -x19*(V{0.0138888888888889}*x167*x63 + V{0.0277777777777778}) - x20*(-x141*x63 + x159);
auto x16 = -x19*(V{0.0138888888888889}*x107*x86 + V{0.0277777777777778}) - x20*(x141*x86 + x160);
auto x17 = -x19*(V{0.0138888888888889}*x167*x73 + V{0.0277777777777778}) - x20*(-x141*x73 + x164);
auto x18 = -x19*(V{0.0138888888888889}*x107*x46 + V{0.0277777777777778}) - x20*(x111*x165 + x166 + x46*(V{0.00810185185185185}*x24*x29 - V{0.00810185185185185}*x26*x28));
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
return { V{0.5}*x107, x37 + x39 + x41 };
}
};

}

}

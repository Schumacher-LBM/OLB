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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<2, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<2, 1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{0.5}*cell[10];
auto x22 = V{0.5}*cell[2];
auto x23 = V{2}*cell[14];
auto x24 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1});
auto x25 = V{0.25}*cell[0] + V{0.25}*cell[12] + V{0.5}*cell[14] + V{0.25}*cell[3] + V{0.25};
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x27 = cell[0] + cell[12] + cell[3] + x23 + V{1};
auto x28 = V{2}*cell[10] + cell[11] + V{2}*cell[13] + V{2}*cell[15] + V{2}*cell[16] + cell[17] + cell[18] + cell[2] + cell[8] + cell[9] + x27;
auto x29 = cell[10] + cell[15] + cell[16] + cell[1] + V{2}*cell[2] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[8] + V{2}*cell[9] + x27;
auto x30 = -V{0.0138888888888889}*x24*x29 + V{0.0138888888888889}*x26*x28;
auto x31 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x32 = -x31;
auto x33 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x34 = -V{4.5}*x33*x33;
auto x35 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x36 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x37 = V{1.5}*x36;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x39 = V{1.5}*x38;
auto x40 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x41 = V{1.5}*x40;
auto x42 = x39 + x41 + V{-1};
auto x43 = x37 + x42;
auto x44 = x35 + x43;
auto x45 = x32 + x34 + x44;
auto x46 = x43*(-V{0.166666666666667}*x24*x29 + V{0.166666666666667}*x26*x28);
auto x47 = -V{0.25}*x24*x29 + V{0.25}*x26*x28;
auto x48 = -x33;
auto x49 = -x35 + x43;
auto x50 = x31 + x49;
auto x51 = x50 - V{4.5}*x48*x48;
auto x52 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x53 = V{4.5}*(x52*x52);
auto x54 = -x37;
auto x55 = V{1} - x39;
auto x56 = x54 + x55;
auto x57 = x35 + x56;
auto x58 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x59 = -x41;
auto x60 = x58 + x59;
auto x61 = x53 + x57 + x60;
auto x62 = -V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28;
auto x63 = V{3}*x40;
auto x64 = x57 + x63;
auto x65 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x66 = V{4.5}*(x65*x65);
auto x67 = x31 + x59;
auto x68 = x57 + x66 + x67;
auto x69 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x70 = V{4.5}*(x69*x69);
auto x71 = x31 + x56 + x60 + x70;
auto x72 = -V{0.0555555555555556}*x24*x29 + V{0.0555555555555556}*x26*x28;
auto x73 = V{3}*x38;
auto x74 = x54 + x67 + x73 + V{1};
auto x75 = -V{1.11022302462516e-16}*x24*x29 + V{1.11022302462516e-16}*x26*x28;
auto x76 = -V{8.32667268468867e-17}*x24*x29 + V{8.32667268468867e-17}*x26*x28;
auto x77 = x44 - x53 + x58;
auto x78 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x79 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x78;
auto x80 = -x79;
auto x81 = -x58;
auto x82 = x44 + x81;
auto x83 = x82 - V{4.5}*x80*x80;
auto x84 = -V{4.5}*x79*x79;
auto x85 = x49 + x58 + x84;
auto x86 = x37 + V{-1};
auto x87 = x35 + x39 - x63 + x86;
auto x88 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x78;
auto x89 = -V{4.5}*x88*x88;
auto x90 = x43 + x58;
auto x91 = x32 + x89 + x90;
auto x92 = V{1.11022302462516e-16}*cell[4];
auto x93 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[14] + V{5.55111512312578e-17}*cell[3] + V{5.55111512312578e-17};
auto x94 = V{1.11022302462516e-16}*cell[13];
auto x95 = V{3}*x36;
auto x96 = x55 + x60 + x95;
auto x97 = x42 + x58 - x95;
auto x98 = V{6.66133814775094e-16}*cell[14] + x24*(V{5.55111512312578e-17}*cell[10] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{5.55111512312578e-17}*cell[1] + V{1.11022302462516e-16}*cell[2] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + V{1.11022302462516e-16}*cell[8] + V{1.11022302462516e-16}*cell[9] + x92 + x93) - x26*(V{1.11022302462516e-16}*cell[10] + V{5.55111512312578e-17}*cell[11] + V{1.11022302462516e-16}*cell[15] + V{1.11022302462516e-16}*cell[16] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{5.55111512312578e-17}*cell[2] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x93 + x94) - x46 - x51*(-V{4.62592926927149e-18}*x24*x29 + V{4.62592926927149e-18}*x26*x28) + x62*x96 - x62*x97 + V{-2.22044604925031e-16};
auto x99 = x26*(V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x30*x61 - x30*x77 - x30*x83 - x30*x85 - x36*(-V{5.55111512312578e-17}*x24*x29 + V{5.55111512312578e-17}*x26*x28) - x38*x75 - x40*x76 - x45*(-V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28) + x62*x64 + x62*x68 + x62*x71 - x62*x87 - x62*x91 + x72*x74 + x98);
auto x100 = x62*x74;
auto x101 = x31 - x70 + x90;
auto x102 = -x88;
auto x103 = x31 + x43 + x81;
auto x104 = x103 - V{4.5}*x102*x102;
auto x105 = x31 + x41 - x73 + x86;
auto x106 = x105*x62;
auto x107 = x31 + x44 - x66;
auto x108 = x24*(V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] + x100 - x101*x30 - x104*x30 - x106 - x107*x62 + x30*x71 - x30*x91 - x36*(-V{4.85722573273506e-17}*x24*x29 + V{4.85722573273506e-17}*x26*x28) - x38*x76 - x40*x75 - x45*(-V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28) - x62*x77 - x62*x83 - x72*x87 + x92 + x94 + x98);
auto x109 = -x108 + x99;
auto x110 = V{0.0833333333333333}*cell[11];
auto x111 = V{0.0833333333333333}*cell[2];
auto x112 = V{0.166666666666667}*cell[10];
auto x113 = V{0.166666666666667}*cell[1];
auto x114 = V{0.0833333333333334}*cell[15];
auto x115 = V{0.0833333333333334}*cell[16];
auto x116 = V{0.0833333333333334}*cell[6];
auto x117 = V{0.0833333333333334}*cell[7];
auto x118 = -V{0.00115740740740741}*x24*x29 + V{0.00115740740740741}*x26*x28;
auto x119 = x118*x45;
auto x120 = -V{0.0416666666666667}*x24*x29 + V{0.0416666666666667}*x26*x28;
auto x121 = x120*x40;
auto x122 = x120*x36;
auto x123 = V{0.0833333333333333}*cell[12];
auto x124 = V{0.0833333333333333}*cell[3];
auto x125 = V{0.0833333333333334}*cell[13];
auto x126 = V{0.0833333333333334}*cell[4];
auto x127 = V{0.166666666666667}*cell[14];
auto x128 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[14] + V{3.46944695195361e-18}*cell[3] + V{3.46944695195361e-18};
auto x129 = x26*(V{6.93889390390723e-18}*cell[10] + V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[2] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x128);
auto x130 = x24*(V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[2] + V{6.93889390390723e-18}*cell[4] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] + x128);
auto x131 = x123 + x124 - x125 - x126 - x127 + x129 - x130 + V{0.0555555555555555};
auto x132 = V{0.166666666666667}*cell[15];
auto x133 = V{0.166666666666667}*cell[16];
auto x134 = V{0.166666666666667}*cell[6];
auto x135 = V{0.166666666666667}*cell[7];
auto x136 = V{0.0277777777777778}*x24*x29 - V{0.0277777777777778}*x26*x28;
auto x137 = V{0.0833333333333333}*x24*x29 - V{0.0833333333333333}*x26*x28;
auto x138 = V{0.00115740740740741}*x24*x29 - V{0.00115740740740741}*x26*x28;
auto x139 = x34 + x50;
auto x140 = x138*x45;
auto x141 = V{0.0416666666666667}*x24*x29 - V{0.0416666666666667}*x26*x28;
auto x142 = x141*x36;
auto x143 = -x129;
auto x144 = -x123 - x124 + x125 + x126 + x127 + x130 + x143 + V{-0.0555555555555555};
auto x145 = V{0.0833333333333333}*cell[10];
auto x146 = V{0.0833333333333333}*cell[1];
auto x147 = x141*x38;
auto x148 = V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] - x145 - x146 - x147;
auto x149 = V{0.00231481481481481}*x24*x29 - V{0.00231481481481481}*x26*x28;
auto x150 = -x110 - x111 + x114 + x115 + x116 + x117;
auto x151 = V{0.166666666666667}*cell[12] - V{0.166666666666667}*cell[13] - V{0.333333333333333}*cell[14] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[4] + x130 + x137*x36 - x139*x149 - x141*x40 + x143 + x148 + x149*x45 + x150 + V{-0.0555555555555555};
auto x152 = V{0.0138888888888889}*x24*x29 - V{0.0138888888888889}*x26*x28;
auto x153 = V{0.00115740740740741}*x24*x29 - V{0.00115740740740741}*x26*x28;
auto x154 = x24*x29 - x26*x28;
auto x155 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x154;
auto x156 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x155;
auto x157 = V{0.0208333333333333}*x24*x29 - V{0.0208333333333333}*x26*x28;
auto x158 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[12] + V{0.0416666666666667}*cell[14] + V{0.0208333333333333}*cell[3] + V{0.0208333333333333};
auto x159 = x24*(V{0.0208333333333333}*cell[10] + V{0.0208333333333333}*cell[15] + V{0.0208333333333333}*cell[16] + V{0.0208333333333333}*cell[1] + V{0.0416666666666667}*cell[2] + V{0.0416666666666667}*cell[4] + V{0.0208333333333333}*cell[6] + V{0.0208333333333333}*cell[7] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x158);
auto x160 = -x26*(V{0.0416666666666667}*cell[10] + V{0.0208333333333333}*cell[11] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0208333333333333}*cell[17] + V{0.0208333333333333}*cell[18] + V{0.0208333333333333}*cell[2] + V{0.0208333333333333}*cell[8] + V{0.0208333333333333}*cell[9] + x158);
auto x161 = V{0.0416666666666667}*x24*x29 - V{0.0416666666666667}*x26*x28;
auto x162 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x159 + x160 + x161*x38 + V{0.0138888888888889};
auto x163 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x161*x40;
auto x164 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x157*x36 + x162 + x163;
auto x165 = V{0.416666666666667}*cell[13] - V{0.166666666666667}*cell[14] + V{0.416666666666667}*cell[4] - x139*x153 + x153*x45 + x156 + x164;
auto x166 = V{0.00578703703703704}*x24*x29 - V{0.00578703703703704}*x26*x28;
auto x167 = -V{0.0833333333333333}*cell[13] + V{0.833333333333333}*cell[14] - V{0.0833333333333333}*cell[4] - x156 + x164;
auto x168 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x155;
auto x169 = V{0.000578703703703704}*x24*x29 - V{0.000578703703703704}*x26*x28;
auto x170 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0833333333333334}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + x139*x169 + x161*x36 - x169*x45;
auto x171 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x157*x40 + x162 + x170;
auto x172 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x168 + x171;
auto x173 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x168 + x171;
auto x174 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x154;
auto x175 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x157*x38 + x159 + x160 + x163 + x170 + V{0.0138888888888889};
auto x176 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] + x174 + x175;
auto x177 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] - x174 + x175;
auto x178 = -x154;
auto x179 = x108 - x99;
auto x0 = -x19*(V{0.166666666666667}*x109*x43 + V{0.333333333333333}) + x20*(V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21 + x22 + x23 + x24*(V{0.25}*cell[10] + V{0.25}*cell[15] + V{0.25}*cell[16] + V{0.25}*cell[1] + V{0.5}*cell[4] + V{0.25}*cell[6] + V{0.25}*cell[7] + V{0.5}*cell[8] + V{0.5}*cell[9] + x22 + x25) - x26*(V{0.25}*cell[11] + V{0.5}*cell[13] + V{0.5}*cell[15] + V{0.5}*cell[16] + V{0.25}*cell[17] + V{0.25}*cell[18] + V{0.25}*cell[2] + V{0.25}*cell[8] + V{0.25}*cell[9] + x21 + x25) + x30*x45 - x30*x51 - x36*x47 - x38*x47 - x40*x47 + x46 + V{0.833333333333333});
auto x1 = -x19*(V{0.0277777777777778}*x105*x109 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x106 + x110 + x111 - x112 - x113 - x114 - x115 - x116 - x117 + x118*x51 - x119 - x121 - x122 + x131 + x38*(-V{0.0833333333333333}*x24*x29 + V{0.0833333333333333}*x26*x28));
auto x2 = -x19*(V{0.0277777777777778}*x109*x87 + V{0.0555555555555556}) - x20*(V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] - x132 - x133 - x134 - x135 + x136*x87 + x137*x40 + x138*x139 - x140 - x142 + x144 + x148);
auto x3 = -x19*(V{0.0277777777777778}*x109*x97 + V{0.0555555555555556}) - x20*(x136*x97 + x151);
auto x4 = -x19*(V{0.0138888888888889}*x107*x109 + V{0.0277777777777778}) - x20*(x107*x152 + x165);
auto x5 = -x19*(V{0.0138888888888889}*x109*x51 + V{0.0277777777777778}) - x20*(x139*(V{0.0196759259259259}*x24*x29 - V{0.0196759259259259}*x26*x28) - x166*x45 + x167);
auto x6 = -x19*(V{0.0138888888888889}*x101*x109 + V{0.0277777777777778}) - x20*(x101*x152 + x172);
auto x7 = -x19*(V{0.0138888888888889}*x104*x109 + V{0.0277777777777778}) - x20*(x152*(x103 + x89) + x173);
auto x8 = -x19*(V{0.0138888888888889}*x109*x77 + V{0.0277777777777778}) - x20*(x152*x77 + x176);
auto x9 = -x19*(V{0.0138888888888889}*x109*x83 + V{0.0277777777777778}) - x20*(x152*(x82 + x84) + x177);
auto x10 = -x19*(V{0.0277777777777778}*x179*x74 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x100 - x112 - x113 - x119 - x121 - x122 - x144 - x150 + V{0.0833333333333333}*x178*x38 + V{0.00115740740740741}*x178*x51);
auto x11 = -x19*(V{0.0277777777777778}*x179*x64 + V{0.0555555555555556}) - x20*(V{0.166666666666667}*cell[11] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.166666666666667}*cell[2] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] - x131 - x132 - x133 - x134 - x135 - x136*x64 + V{0.00115740740740741}*x139*x154 - x140 - x142 - x145 - x146 - x147 + V{0.0833333333333333}*x154*x40);
auto x12 = -x19*(V{0.0277777777777778}*x179*x96 + V{0.0555555555555556}) - x20*(-x136*x96 + x151);
auto x13 = -x19*(V{0.0138888888888889}*x179*x68 + V{0.0277777777777778}) - x20*(-x152*x68 + x165);
auto x14 = -x19*(V{0.0138888888888889}*x109*x45 + V{0.0277777777777778}) - x20*(x139*x166 + x167 + x45*(V{0.00810185185185185}*x24*x29 - V{0.00810185185185185}*x26*x28));
auto x15 = -x19*(V{0.0138888888888889}*x179*x71 + V{0.0277777777777778}) - x20*(-x152*x71 + x172);
auto x16 = -x19*(V{0.0138888888888889}*x109*x91 + V{0.0277777777777778}) - x20*(x152*x91 + x173);
auto x17 = -x19*(V{0.0138888888888889}*x179*x61 + V{0.0277777777777778}) - x20*(-x152*x61 + x176);
auto x18 = -x19*(V{0.0138888888888889}*x109*x85 + V{0.0277777777777778}) - x20*(x152*x85 + x177);
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
return { V{0.5}*x109, x36 + x38 + x40 };
}
};

}

}

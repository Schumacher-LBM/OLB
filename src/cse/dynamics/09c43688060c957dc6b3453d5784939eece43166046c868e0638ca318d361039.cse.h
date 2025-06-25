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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<1, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<1, -1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{0.5}*cell[1];
auto x22 = V{0.5}*cell[3];
auto x23 = V{2}*cell[6];
auto x24 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1});
auto x25 = V{0.25}*cell[0] + V{0.25}*cell[11] + V{0.25}*cell[2] + V{0.5}*cell[6] + V{0.25};
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1});
auto x27 = cell[0] + cell[11] + cell[2] + x23 + V{1};
auto x28 = cell[12] + cell[17] + cell[18] + V{2}*cell[1] + cell[3] + V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[7] + cell[8] + cell[9] + x27;
auto x29 = cell[10] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[18] + cell[1] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[8] + x27;
auto x30 = V{0.25}*x24*x28 + V{0.25}*x26*x29;
auto x31 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x33 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x34 = V{0.0138888888888889}*x24*x28 + V{0.0138888888888889}*x26*x29;
auto x35 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x36 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x37 = V{4.5}*(x36*x36);
auto x38 = V{1.5}*x33;
auto x39 = -x38;
auto x40 = V{1.5}*x31;
auto x41 = V{1} - x40;
auto x42 = x39 + x41;
auto x43 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x44 = V{1.5}*x32;
auto x45 = -x44;
auto x46 = x43 + x45;
auto x47 = x35 + x37 + x42 + x46;
auto x48 = x40 + x44 + V{-1};
auto x49 = x38 + x48;
auto x50 = x43 + x49;
auto x51 = x35 - x37 + x50;
auto x52 = x49*(V{0.166666666666667}*x24*x28 + V{0.166666666666667}*x26*x29);
auto x53 = V{1.11022302462516e-16}*cell[3];
auto x54 = V{0.0277777777777778}*x24*x28 + V{0.0277777777777778}*x26*x29;
auto x55 = V{3}*x33;
auto x56 = x43 + x48 - x55;
auto x57 = x54*x56;
auto x58 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x59 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x60 = V{4.5}*(x59*x59);
auto x61 = x50 + x58 - x60;
auto x62 = x34*x61;
auto x63 = x42 + x58;
auto x64 = x46 + x60 + x63;
auto x65 = -x34*x64;
auto x66 = -x58;
auto x67 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x68 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x67;
auto x69 = -V{4.5}*x68*x68;
auto x70 = x50 + x66 + x69;
auto x71 = x34*x70;
auto x72 = V{5.55111512312578e-17}*x24*x28 + V{5.55111512312578e-17}*x26*x29;
auto x73 = V{1.11022302462516e-16}*x24*x28 + V{1.11022302462516e-16}*x26*x29;
auto x74 = V{8.32667268468867e-17}*x24*x28 + V{8.32667268468867e-17}*x26*x29;
auto x75 = -x68;
auto x76 = -x43;
auto x77 = x49 + x58;
auto x78 = x76 + x77;
auto x79 = x78 - V{4.5}*x75*x75;
auto x80 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x81 = V{4.5}*(x80*x80);
auto x82 = x35 + x77 - x81;
auto x83 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x84 = -x83;
auto x85 = x35 + x49;
auto x86 = x66 + x85;
auto x87 = x86 - V{4.5}*x84*x84;
auto x88 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x67;
auto x89 = -x88;
auto x90 = x76 + x85;
auto x91 = x90 - V{4.5}*x89*x89;
auto x92 = V{0.0555555555555556}*x24*x28 + V{0.0555555555555556}*x26*x29;
auto x93 = V{3}*x31;
auto x94 = x38 + V{-1};
auto x95 = x35 + x44 - x93 + x94;
auto x96 = x41 + x46 + x55;
auto x97 = x54*x96;
auto x98 = V{3}*x32;
auto x99 = x40 + x58 + x94 - x98;
auto x100 = x54*x99;
auto x101 = x63 + x98;
auto x102 = -x101*x54;
auto x103 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[2] + V{1.11022302462516e-16}*cell[6] + V{5.55111512312578e-17};
auto x104 = x100 + x102 + x24*(V{5.55111512312578e-17}*cell[12] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{1.11022302462516e-16}*cell[1] + V{5.55111512312578e-17}*cell[3] + V{1.11022302462516e-16}*cell[4] + V{1.11022302462516e-16}*cell[5] + V{1.11022302462516e-16}*cell[7] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x103) + x26*(V{5.55111512312578e-17}*cell[10] + V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[16] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[1] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[8] + x103 + x53) + x51*(V{0.0277777777777778}*x24*x28 + V{0.0277777777777778}*x26*x29) + x52 + V{-2.22044604925031e-16};
auto x105 = x34*x82;
auto x106 = x35 + x45;
auto x107 = x106 + x63 + x81;
auto x108 = -x107*x34;
auto x109 = -x35;
auto x110 = -V{4.5}*x83*x83;
auto x111 = x109 + x110 + x77;
auto x112 = x111*x34;
auto x113 = x54*x95;
auto x114 = -V{4.5}*x88*x88;
auto x115 = x109 + x114 + x50;
auto x116 = x106 + x39 + x93 + V{1};
auto x117 = x116*x54;
auto x118 = x24*(V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + V{8.88178419700125e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x104 + x31*x73 + x32*x74 + x33*x72 + x34*x79 - x47*(V{6.16790569236198e-18}*x24*x28 + V{6.16790569236198e-18}*x26*x29) + x53 + x54*x82 + x54*x87 + x54*x91 + x57 + x62 + x65 + x71 + x92*x95 - x97) + x26*(V{1.66533453693773e-16}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{6.66133814775094e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + V{2.22044604925031e-16}*cell[9] + x104 + x105 + x108 + x112 + x113 + x115*x54 - x117 + x31*x74 + x32*x72 + x33*x73 + x34*x87 - x47*(V{4.62592926927149e-18}*x24*x28 + V{4.62592926927149e-18}*x26*x29) + x54*x61 + x54*x70 + x56*x92);
auto x119 = V{0.0833333333333333}*cell[12];
auto x120 = V{0.0833333333333333}*cell[3];
auto x121 = V{0.0833333333333334}*cell[13];
auto x122 = V{0.0833333333333334}*cell[14];
auto x123 = V{0.0833333333333334}*cell[4];
auto x124 = V{0.0833333333333334}*cell[5];
auto x125 = V{0.0416666666666667}*x24*x28 + V{0.0416666666666667}*x26*x29;
auto x126 = x125*x33;
auto x127 = V{0.0833333333333333}*x24*x28 + V{0.0833333333333333}*x26*x29;
auto x128 = V{0.0833333333333333}*cell[11];
auto x129 = V{0.0833333333333333}*cell[2];
auto x130 = V{0.0833333333333334}*cell[16];
auto x131 = V{0.0833333333333334}*cell[7];
auto x132 = V{0.166666666666667}*cell[6];
auto x133 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[6] + V{3.46944695195361e-18};
auto x134 = x24*(V{3.46944695195361e-18}*cell[12] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] + V{6.93889390390723e-18}*cell[7] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x133);
auto x135 = x26*(V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[16] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[3] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[8] + x133);
auto x136 = V{0.00115740740740741}*x24*x28 + V{0.00115740740740741}*x26*x29;
auto x137 = x136*x47;
auto x138 = x136*x51;
auto x139 = x125*x32;
auto x140 = x128 + x129 - x130 - x131 - x132 - x134 - x135 + x137 + x138 + x139 + V{0.0555555555555555};
auto x141 = -V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x119 + x120 - x121 - x122 - x123 - x124 + x126 - x127*x31 + x140;
auto x142 = V{0.00231481481481481}*x24*x28 + V{0.00231481481481481}*x26*x29;
auto x143 = V{0.0833333333333333}*cell[10];
auto x144 = V{0.0833333333333333}*cell[1];
auto x145 = x125*x31;
auto x146 = V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x134 + x135 - x143 - x144 - x145 + V{-0.0555555555555555};
auto x147 = V{0.166666666666667}*cell[11] - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.166666666666667}*cell[7] - x119 - x120 + x121 + x122 + x123 + x124 - x126 + x127*x32 + x142*x47 + x142*x51 + x146;
auto x148 = V{0.166666666666667}*cell[13];
auto x149 = V{0.166666666666667}*cell[14];
auto x150 = V{0.166666666666667}*cell[4];
auto x151 = V{0.166666666666667}*cell[5];
auto x152 = x24*x28 + x26*x29;
auto x153 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x152;
auto x154 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x153;
auto x155 = V{0.0208333333333333}*x24*x28 + V{0.0208333333333333}*x26*x29;
auto x156 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[11] + V{0.0208333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0208333333333333};
auto x157 = x24*(V{0.0208333333333333}*cell[12] + V{0.0208333333333333}*cell[17] + V{0.0208333333333333}*cell[18] + V{0.0416666666666667}*cell[1] + V{0.0208333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + V{0.0416666666666667}*cell[7] + V{0.0208333333333333}*cell[8] + V{0.0208333333333333}*cell[9] + x156);
auto x158 = x26*(V{0.0208333333333333}*cell[10] + V{0.0208333333333333}*cell[13] + V{0.0208333333333333}*cell[14] + V{0.0416666666666667}*cell[16] + V{0.0416666666666667}*cell[18] + V{0.0208333333333333}*cell[1] + V{0.0416666666666667}*cell[3] + V{0.0208333333333333}*cell[4] + V{0.0208333333333333}*cell[5] + V{0.0416666666666667}*cell[8] + x156);
auto x159 = V{0.0416666666666667}*x24*x28 + V{0.0416666666666667}*x26*x29;
auto x160 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x157 + x158 + x159*x31 + V{0.0138888888888889};
auto x161 = V{0.000578703703703704}*x24*x28 + V{0.000578703703703704}*x26*x29;
auto x162 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0416666666666667}*cell[7] + x159*x32 - x161*x47 - x161*x51;
auto x163 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x155*x33 + x160 + x162;
auto x164 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x154 + x163;
auto x165 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x154 + x163;
auto x166 = V{0.00578703703703704}*x24*x28 + V{0.00578703703703704}*x26*x29;
auto x167 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x153;
auto x168 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + x159*x33;
auto x169 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x155*x32 + x160 + x168;
auto x170 = -V{0.0833333333333333}*cell[16] + V{0.833333333333333}*cell[6] - V{0.0833333333333333}*cell[7] + x167 + x169;
auto x171 = V{0.00115740740740741}*x24*x28 + V{0.00115740740740741}*x26*x29;
auto x172 = V{0.416666666666667}*cell[16] - V{0.166666666666667}*cell[6] + V{0.416666666666667}*cell[7] - x167 + x169 + x171*x47 + x171*x51;
auto x173 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x152;
auto x174 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x155*x31 + x157 + x158 + x162 + x168 + V{0.0138888888888889};
auto x175 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] + x173 + x174;
auto x176 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] - x173 + x174;
auto x0 = -x19*(-V{0.166666666666667}*x118*x49 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[2] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21 + x22 + x23 + x24*(V{0.25}*cell[12] + V{0.25}*cell[17] + V{0.25}*cell[18] + V{0.25}*cell[3] + V{0.5}*cell[4] + V{0.5}*cell[5] + V{0.5}*cell[7] + V{0.25}*cell[8] + V{0.25}*cell[9] + x21 + x25) + x26*(V{0.25}*cell[10] + V{0.25}*cell[13] + V{0.25}*cell[14] + V{0.5}*cell[16] + V{0.5}*cell[18] + V{0.25}*cell[1] + V{0.25}*cell[4] + V{0.25}*cell[5] + V{0.5}*cell[8] + x22 + x25) + x30*x31 + x30*x32 + x30*x33 - x34*x47 - x34*x51 - x52 + V{0.833333333333333});
auto x1 = -x19*(-V{0.0277777777777778}*x118*x95 + V{0.0555555555555556}) + x20*(-x113 + x141);
auto x2 = -x19*(-V{0.0277777777777778}*x118*x99 + V{0.0555555555555556}) - x20*(x100 + x147);
auto x3 = -x19*(-V{0.0277777777777778}*x118*x56 + V{0.0555555555555556}) - x20*(V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + x127*x33 - x128 - x129 + x130 + x131 + x132 - x137 - x138 - x139 + x146 - x148 - x149 - x150 - x151 + x57);
auto x4 = -x19*(-V{0.0138888888888889}*x118*x82 + V{0.0277777777777778}) - x20*(x105 + x164);
auto x5 = -x19*(-V{0.0138888888888889}*x118*x87 + V{0.0277777777777778}) - x20*(x165 + x34*(x110 + x86));
auto x6 = -x19*(-V{0.0138888888888889}*x118*x51 + V{0.0277777777777778}) - x20*(-x166*x47 + x170 + x51*(V{0.00810185185185185}*x24*x28 + V{0.00810185185185185}*x26*x29));
auto x7 = -x19*(-V{0.0138888888888889}*x118*x91 + V{0.0277777777777778}) - x20*(x172 + x34*(x114 + x90));
auto x8 = -x19*(-V{0.0138888888888889}*x118*x61 + V{0.0277777777777778}) - x20*(x175 + x62);
auto x9 = -x19*(-V{0.0138888888888889}*x118*x79 + V{0.0277777777777778}) - x20*(x176 + x34*(x69 + x78));
auto x10 = -x19*(V{0.0277777777777778}*x116*x118 + V{0.0555555555555556}) + x20*(x117 + x141);
auto x11 = -x19*(V{0.0277777777777778}*x101*x118 + V{0.0555555555555556}) - x20*(x102 + x147);
auto x12 = -x19*(V{0.0277777777777778}*x118*x96 + V{0.0555555555555556}) - x20*(V{0.166666666666667}*cell[12] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] + V{0.166666666666667}*cell[3] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] - x140 - x143 - x144 - x145 - x148 - x149 - x150 - x151 + V{0.0833333333333333}*x152*x33 - x97);
auto x13 = -x19*(V{0.0138888888888889}*x107*x118 + V{0.0277777777777778}) - x20*(x108 + x164);
auto x14 = -x19*(-V{0.0138888888888889}*x111*x118 + V{0.0277777777777778}) - x20*(x112 + x165);
auto x15 = -x19*(V{0.0138888888888889}*x118*x47 + V{0.0277777777777778}) - x20*(-x166*x51 + x170 - x47*(V{0.0196759259259259}*x24*x28 + V{0.0196759259259259}*x26*x29));
auto x16 = -x19*(-V{0.0138888888888889}*x115*x118 + V{0.0277777777777778}) - x20*(x115*x34 + x172);
auto x17 = -x19*(V{0.0138888888888889}*x118*x64 + V{0.0277777777777778}) - x20*(x175 + x65);
auto x18 = -x19*(-V{0.0138888888888889}*x118*x70 + V{0.0277777777777778}) - x20*(x176 + x71);
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
return { -V{0.5}*x118, x31 + x32 + x33 };
}
};

}

}

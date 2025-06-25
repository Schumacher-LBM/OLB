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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<2, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<2, -1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{0.5}*cell[11];
auto x22 = V{0.5}*cell[1];
auto x23 = V{2}*cell[5];
auto x24 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1});
auto x25 = V{0.25}*cell[0] + V{0.25}*cell[12] + V{0.25}*cell[3] + V{0.5}*cell[5] + V{0.25};
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x27 = cell[0] + cell[12] + cell[3] + x23 + V{1};
auto x28 = cell[10] + V{2}*cell[11] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[18] + cell[1] + cell[6] + cell[7] + x27;
auto x29 = cell[11] + cell[17] + cell[18] + V{2}*cell[1] + cell[2] + V{2}*cell[4] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + cell[9] + x27;
auto x30 = -V{0.0138888888888889}*x24*x29 + V{0.0138888888888889}*x26*x28;
auto x31 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x32 = -x31;
auto x33 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x34 = -x33;
auto x35 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x36 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x37 = V{1.5}*x36;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x39 = V{1.5}*x38;
auto x40 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x41 = V{1.5}*x40;
auto x42 = x39 + x41 + V{-1};
auto x43 = x37 + x42;
auto x44 = x35 + x43;
auto x45 = x34 + x44;
auto x46 = x45 - V{4.5}*x32*x32;
auto x47 = x43*(-V{0.166666666666667}*x24*x29 + V{0.166666666666667}*x26*x28);
auto x48 = -V{0.25}*x24*x29 + V{0.25}*x26*x28;
auto x49 = -V{4.5}*x31*x31;
auto x50 = -x35 + x43;
auto x51 = x33 + x49 + x50;
auto x52 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x53 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x52;
auto x54 = -x53;
auto x55 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x56 = -x55;
auto x57 = x44 + x56;
auto x58 = x57 - V{4.5}*x54*x54;
auto x59 = -V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28;
auto x60 = -V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28;
auto x61 = V{3}*x36;
auto x62 = x42 + x55 - x61;
auto x63 = -x60*x62;
auto x64 = V{1} - x41;
auto x65 = -x39;
auto x66 = x55 + x65;
auto x67 = x61 + x64 + x66;
auto x68 = x60*x67;
auto x69 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[12] + V{5.55111512312578e-17}*cell[3] + V{1.11022302462516e-16}*cell[5] + V{5.55111512312578e-17};
auto x70 = V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{6.66133814775094e-16}*cell[5] + x24*(V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{1.11022302462516e-16}*cell[1] + V{5.55111512312578e-17}*cell[2] + V{1.11022302462516e-16}*cell[4] + V{1.11022302462516e-16}*cell[6] + V{1.11022302462516e-16}*cell[7] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x69) - x26*(V{5.55111512312578e-17}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{1.11022302462516e-16}*cell[13] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{1.11022302462516e-16}*cell[17] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[1] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + x69) - x36*(-V{5.55111512312578e-17}*x24*x29 + V{5.55111512312578e-17}*x26*x28) - x47 - x51*(-V{4.62592926927149e-18}*x24*x29 + V{4.62592926927149e-18}*x26*x28) + x63 + x68 + V{-2.22044604925031e-16};
auto x71 = -x46*x59 + x70;
auto x72 = V{3}*x38;
auto x73 = -x37;
auto x74 = x64 + x73;
auto x75 = x35 + x74;
auto x76 = x72 + x75;
auto x77 = x60*x76;
auto x78 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x79 = V{4.5}*(x78*x78);
auto x80 = x44 + x55 - x79;
auto x81 = -x30*x80;
auto x82 = x66 + x75 + x79;
auto x83 = x30*x82;
auto x84 = -V{4.5}*x53*x53;
auto x85 = x50 + x55 + x84;
auto x86 = -x30*x85;
auto x87 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x88 = V{4.5}*(x87*x87);
auto x89 = x33 + x65;
auto x90 = x75 + x88 + x89;
auto x91 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x92 = V{4.5}*(x91*x91);
auto x93 = x33 + x66 + x74 + x92;
auto x94 = -V{0.0555555555555556}*x24*x29 + V{0.0555555555555556}*x26*x28;
auto x95 = V{3}*x40;
auto x96 = x73 + x89 + x95 + V{1};
auto x97 = -V{1.11022302462516e-16}*x24*x29 + V{1.11022302462516e-16}*x26*x28;
auto x98 = -V{8.32667268468867e-17}*x24*x29 + V{8.32667268468867e-17}*x26*x28;
auto x99 = x37 + V{-1};
auto x100 = x35 + x41 - x72 + x99;
auto x101 = x100*x60;
auto x102 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x52;
auto x103 = -V{4.5}*x102*x102;
auto x104 = x43 + x55;
auto x105 = x103 + x104 + x34;
auto x106 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] - x101 - x105*x60 - x38*x98 - x40*x97 + x60*x90 + x60*x93 + x77 + x81 + x83 + x86 + x94*x96;
auto x107 = -x102;
auto x108 = x33 + x43 + x56;
auto x109 = x108 - V{4.5}*x107*x107;
auto x110 = x33 + x39 - x95 + x99;
auto x111 = -x110*x60;
auto x112 = x60*x96;
auto x113 = x104 + x33 - x92;
auto x114 = -x113*x30;
auto x115 = x30*x93;
auto x116 = -x105*x30;
auto x117 = x33 + x44 - x88;
auto x118 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] - x100*x94 + x111 + x112 + x114 + x115 + x116 - x117*x60 - x38*x97 - x40*x98 - x60*x80;
auto x119 = -x24*(-x109*x30 + x118 - x58*x60 + x71) + x26*(x106 - x30*x58 + x71);
auto x120 = V{0.0833333333333333}*cell[11];
auto x121 = V{0.0833333333333333}*cell[12];
auto x122 = V{0.0833333333333333}*cell[2];
auto x123 = V{0.0833333333333333}*cell[3];
auto x124 = V{0.166666666666667}*cell[10];
auto x125 = V{0.166666666666667}*cell[1];
auto x126 = V{0.0833333333333334}*cell[13];
auto x127 = V{0.0833333333333334}*cell[15];
auto x128 = V{0.0833333333333334}*cell[16];
auto x129 = V{0.0833333333333334}*cell[4];
auto x130 = V{0.0833333333333334}*cell[6];
auto x131 = V{0.0833333333333334}*cell[7];
auto x132 = V{0.166666666666667}*cell[5];
auto x133 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[12] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[5] + V{3.46944695195361e-18};
auto x134 = x26*(V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[11] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[1] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x133);
auto x135 = x24*(V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x133);
auto x136 = -V{0.0833333333333333}*x24*x29 + V{0.0833333333333333}*x26*x28;
auto x137 = -V{0.00115740740740741}*x24*x29 + V{0.00115740740740741}*x26*x28;
auto x138 = x137*x51;
auto x139 = x137*x46;
auto x140 = -V{0.0416666666666667}*x24*x29 + V{0.0416666666666667}*x26*x28;
auto x141 = x140*x40;
auto x142 = x140*x36;
auto x143 = x45 + x49;
auto x144 = -x134;
auto x145 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] + x135 + x140*x38 + x144 + V{-0.0555555555555555};
auto x146 = -x121 - x123 + x126 + x129 + x132 - x138 + x142;
auto x147 = V{0.166666666666667}*cell[11] - V{0.166666666666667}*cell[15] - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[2] - V{0.166666666666667}*cell[6] - V{0.166666666666667}*cell[7] - x136*x40 + x137*x143 + x145 + x146;
auto x148 = -V{0.00231481481481481}*x24*x29 + V{0.00231481481481481}*x26*x28;
auto x149 = -x120 - x122 + x127 + x128 + x130 + x131 + x141;
auto x150 = V{0.166666666666667}*cell[12] - V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[4] - V{0.333333333333333}*cell[5] - x136*x36 - x143*x148 + x145 + x148*x51 + x149;
auto x151 = -V{0.00115740740740741}*x24*x29 + V{0.00115740740740741}*x26*x28;
auto x152 = -x24*x29 + x26*x28;
auto x153 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x152;
auto x154 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x153;
auto x155 = -V{0.0208333333333333}*x24*x29 + V{0.0208333333333333}*x26*x28;
auto x156 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[12] + V{0.0208333333333333}*cell[3] + V{0.0416666666666667}*cell[5] + V{0.0208333333333333};
auto x157 = x24*(V{0.0208333333333333}*cell[11] + V{0.0208333333333333}*cell[17] + V{0.0208333333333333}*cell[18] + V{0.0416666666666667}*cell[1] + V{0.0208333333333333}*cell[2] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + V{0.0208333333333333}*cell[8] + V{0.0208333333333333}*cell[9] + x156);
auto x158 = -x26*(V{0.0208333333333333}*cell[10] + V{0.0416666666666667}*cell[11] + V{0.0416666666666667}*cell[13] + V{0.0208333333333333}*cell[15] + V{0.0208333333333333}*cell[16] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0208333333333333}*cell[1] + V{0.0208333333333333}*cell[6] + V{0.0208333333333333}*cell[7] + x156);
auto x159 = -V{0.0416666666666667}*x24*x29 + V{0.0416666666666667}*x26*x28;
auto x160 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x157 + x158 - x159*x38 + V{0.0138888888888889};
auto x161 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x159*x40;
auto x162 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x155*x36 + x160 + x161;
auto x163 = V{0.416666666666667}*cell[13] + V{0.416666666666667}*cell[4] - V{0.166666666666667}*cell[5] - x143*x151 + x151*x51 - x154 + x162;
auto x164 = -V{0.00578703703703704}*x24*x29 + V{0.00578703703703704}*x26*x28;
auto x165 = -V{0.0833333333333333}*cell[13] - V{0.0833333333333333}*cell[4] + V{0.833333333333333}*cell[5] + x154 + x162;
auto x166 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x153;
auto x167 = -V{0.000578703703703704}*x24*x29 + V{0.000578703703703704}*x26*x28;
auto x168 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0833333333333334}*cell[5] + x143*x167 - x159*x36 - x167*x51;
auto x169 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x155*x40 + x160 + x168;
auto x170 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x166 + x169;
auto x171 = x57 + x84;
auto x172 = -x171*x30;
auto x173 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x166 + x169;
auto x174 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x152;
auto x175 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x155*x38 + x157 + x158 + x161 + x168 + V{0.0138888888888889};
auto x176 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x174 + x175;
auto x177 = -x30*(x103 + x108);
auto x178 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x174 + x175;
auto x179 = -x143*x59 + x70;
auto x180 = -x24*(x118 - x171*x60 + x177 + x179) + x26*(x106 + x172 + x179);
auto x0 = -x19*(V{0.166666666666667}*x119*x43 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21 + x22 + x23 + x24*(V{0.25}*cell[11] + V{0.25}*cell[17] + V{0.25}*cell[18] + V{0.25}*cell[2] + V{0.5}*cell[4] + V{0.5}*cell[6] + V{0.5}*cell[7] + V{0.25}*cell[8] + V{0.25}*cell[9] + x22 + x25) - x26*(V{0.25}*cell[10] + V{0.5}*cell[13] + V{0.25}*cell[15] + V{0.25}*cell[16] + V{0.5}*cell[17] + V{0.5}*cell[18] + V{0.25}*cell[1] + V{0.25}*cell[6] + V{0.25}*cell[7] + x21 + x25) + x30*x46 - x30*x51 - x36*x48 - x38*x48 - x40*x48 + x47 + V{0.833333333333333});
auto x1 = -x19*(V{0.0277777777777778}*x100*x119 + V{0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x101 + x120 + x121 + x122 + x123 - x124 - x125 - x126 - x127 - x128 - x129 - x130 - x131 - x132 + x134 - x135 + x136*x38 + x138 - x139 - x141 - x142 + V{0.0555555555555555});
auto x2 = -x19*(V{0.0277777777777778}*x110*x119 + V{0.0555555555555556}) - x20*(x111 + x147);
auto x3 = -x19*(V{0.0277777777777778}*x119*x62 + V{0.0555555555555556}) - x20*(x150 + x63);
auto x4 = -x19*(V{0.0138888888888889}*x117*x119 + V{0.0277777777777778}) - x20*(-x117*x30 + x163);
auto x5 = -x19*(V{0.0138888888888889}*x119*x46 + V{0.0277777777777778}) - x20*(-x143*(-V{0.00810185185185185}*x24*x29 + V{0.00810185185185185}*x26*x28) - x164*x51 + x165);
auto x6 = -x19*(V{0.0138888888888889}*x119*x80 + V{0.0277777777777778}) - x20*(x170 + x81);
auto x7 = -x19*(V{0.0138888888888889}*x119*x58 + V{0.0277777777777778}) - x20*(x172 + x173);
auto x8 = -x19*(V{0.0138888888888889}*x113*x119 + V{0.0277777777777778}) - x20*(x114 + x176);
auto x9 = -x19*(V{0.0138888888888889}*x109*x119 + V{0.0277777777777778}) - x20*(x177 + x178);
auto x10 = x19*(V{0.0277777777777778}*x180*x76 + V{-0.0555555555555556}) + x20*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] - x124 - x125 - x135 - x139 - x144 - x146 - x149 + V{0.0833333333333333}*x152*x38 - x77 + V{0.0555555555555555});
auto x11 = x19*(V{0.0277777777777778}*x180*x96 + V{-0.0555555555555556}) - x20*(x112 + x147);
auto x12 = x19*(V{0.0277777777777778}*x180*x67 + V{-0.0555555555555556}) - x20*(x150 + x68);
auto x13 = x19*(V{0.0138888888888889}*x180*x90 + V{-0.0277777777777778}) - x20*(x163 + x30*x90);
auto x14 = -x19*(V{0.0138888888888889}*x119*x51 + V{0.0277777777777778}) - x20*(x143*x164 + x165 - x51*(-V{0.0196759259259259}*x24*x29 + V{0.0196759259259259}*x26*x28));
auto x15 = x19*(V{0.0138888888888889}*x180*x82 + V{-0.0277777777777778}) - x20*(x170 + x83);
auto x16 = -x19*(V{0.0138888888888889}*x119*x85 + V{0.0277777777777778}) - x20*(x173 + x86);
auto x17 = x19*(V{0.0138888888888889}*x180*x93 + V{-0.0277777777777778}) - x20*(x115 + x176);
auto x18 = -x19*(V{0.0138888888888889}*x105*x119 + V{0.0277777777777778}) - x20*(x116 + x178);
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
return { V{0.5}*x180, x36 + x38 + x40 };
}
};

}

}

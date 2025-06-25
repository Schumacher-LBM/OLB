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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<-1, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<-1, -1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1});
auto x22 = V{0.166666666666667}*cell[0];
auto x23 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + V{0.333333333333333}*cell[4] + x22 + V{0.166666666666667};
auto x24 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + V{0.333333333333333}*cell[6];
auto x25 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1});
auto x26 = V{0.166666666666667}*cell[15];
auto x27 = V{0.166666666666667}*cell[16];
auto x28 = V{0.166666666666667}*cell[6];
auto x29 = V{0.166666666666667}*cell[7];
auto x30 = V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] + V{0.333333333333333}*cell[8];
auto x31 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1});
auto x32 = V{0.166666666666667}*cell[13];
auto x33 = V{0.166666666666667}*cell[14];
auto x34 = V{0.166666666666667}*cell[4];
auto x35 = V{0.166666666666667}*cell[5];
auto x36 = cell[0] + cell[12] + cell[3] + V{2}*cell[4] + V{1};
auto x37 = cell[11] + cell[2] + V{2}*cell[6];
auto x38 = cell[17] + cell[18] + V{2}*cell[1] + V{2}*cell[5] + V{2}*cell[7] + cell[8] + cell[9] + x36 + x37;
auto x39 = cell[10] + cell[1] + V{2}*cell[8];
auto x40 = V{2}*cell[14] + cell[15] + cell[16] + V{2}*cell[2] + cell[6] + cell[7] + V{2}*cell[9] + x36 + x39;
auto x41 = cell[0] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[18] + V{2}*cell[3] + cell[4] + cell[5] + x37 + x39 + V{1};
auto x42 = V{0.166666666666667}*x21*x38 + V{0.166666666666667}*x25*x40 + V{0.166666666666667}*x31*x41;
auto x43 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x44 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x45 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x46 = V{1.5}*x45;
auto x47 = V{1.5}*x43;
auto x48 = V{1.5}*x44;
auto x49 = x47 + x48 + V{-1};
auto x50 = x46 + x49;
auto x51 = x50*(V{0.111111111111111}*x21*x38 + V{0.111111111111111}*x25*x40 + V{0.111111111111111}*x31*x41);
auto x52 = V{1.11022302462516e-16}*cell[3];
auto x53 = V{0.00925925925925926}*x21*x38 + V{0.00925925925925926}*x25*x40 + V{0.00925925925925926}*x31*x41;
auto x54 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x55 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x56 = V{4.5}*(x55*x55);
auto x57 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x58 = x50 + x57;
auto x59 = x54 - x56 + x58;
auto x60 = x53*x59;
auto x61 = -x46;
auto x62 = V{1} - x47;
auto x63 = x61 + x62;
auto x64 = x57 + x63;
auto x65 = -x48;
auto x66 = x54 + x65;
auto x67 = x56 + x64 + x66;
auto x68 = -x53*x67;
auto x69 = -x57;
auto x70 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x71 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x70;
auto x72 = -V{4.5}*x71*x71;
auto x73 = x50 + x54;
auto x74 = x69 + x72 + x73;
auto x75 = x53*x74;
auto x76 = V{5.55111512312578e-17}*x21*x38 + V{5.55111512312578e-17}*x25*x40 + V{5.55111512312578e-17}*x31*x41;
auto x77 = -x71;
auto x78 = -x54;
auto x79 = x58 + x78;
auto x80 = x79 - V{4.5}*x77*x77;
auto x81 = V{0.0185185185185185}*x21*x38 + V{0.0185185185185185}*x25*x40 + V{0.0185185185185185}*x31*x41;
auto x82 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x83 = -x82;
auto x84 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x85 = x50 + x84;
auto x86 = x69 + x85;
auto x87 = x86 - V{4.5}*x83*x83;
auto x88 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x70;
auto x89 = -x88;
auto x90 = x78 + x85;
auto x91 = x90 - V{4.5}*x89*x89;
auto x92 = V{0.037037037037037}*x21*x38 + V{0.037037037037037}*x25*x40 + V{0.037037037037037}*x31*x41;
auto x93 = V{3}*x43;
auto x94 = x46 + V{-1};
auto x95 = x48 + x84 - x93 + x94;
auto x96 = V{3.70074341541719e-17}*x21*x38 + V{3.70074341541719e-17}*x25*x40 + V{3.70074341541719e-17}*x31*x41;
auto x97 = V{7.40148683083438e-17}*x21*x38 + V{7.40148683083438e-17}*x25*x40 + V{7.40148683083438e-17}*x31*x41;
auto x98 = V{3}*x44;
auto x99 = x47 + x57 + x94 - x98;
auto x100 = x81*x99;
auto x101 = x64 + x98;
auto x102 = -x101*x81;
auto x103 = V{5.55111512312578e-17}*cell[0];
auto x104 = V{1.11022302462516e-16}*cell[4];
auto x105 = V{5.55111512312578e-17}*cell[12] + V{5.55111512312578e-17}*cell[3] + x103 + x104 + V{5.55111512312578e-17};
auto x106 = V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[2] + V{1.11022302462516e-16}*cell[6];
auto x107 = x21*(V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{1.11022302462516e-16}*cell[1] + V{1.11022302462516e-16}*cell[5] + V{1.11022302462516e-16}*cell[7] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x105 + x106);
auto x108 = V{1.11022302462516e-16}*cell[2];
auto x109 = V{5.55111512312578e-17}*cell[10] + V{5.55111512312578e-17}*cell[1] + V{1.11022302462516e-16}*cell[8];
auto x110 = x25*(V{1.11022302462516e-16}*cell[14] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + V{1.11022302462516e-16}*cell[9] + x105 + x108 + x109);
auto x111 = x31*(V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[16] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + x103 + x106 + x109 + x52 + V{5.55111512312578e-17});
auto x112 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x113 = V{4.5}*(x112*x112);
auto x114 = -x113 + x73 + x84;
auto x115 = x100 + x102 + x107 + x110 + x111 + x114*x81 + x51 + V{-2.22044604925031e-16};
auto x116 = V{3}*x45;
auto x117 = -x116 + x49 + x54;
auto x118 = x117*x81;
auto x119 = x116 + x62 + x66;
auto x120 = -x119*x81;
auto x121 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x122 = V{4.5}*(x121*x121);
auto x123 = -x122 + x58 + x84;
auto x124 = V{3.33066907387547e-16}*cell[14] + V{3.33066907387547e-16}*cell[5] + x118 + x120 + x123*x81;
auto x125 = x114*x53;
auto x126 = x113 + x63 + x66 + x84;
auto x127 = -x126*x53;
auto x128 = -x84;
auto x129 = -V{4.5}*x88*x88;
auto x130 = x128 + x129 + x73;
auto x131 = x130*x53;
auto x132 = -V{4.5}*x82*x82;
auto x133 = x128 + x132 + x58;
auto x134 = x81*x95;
auto x135 = x65 + x84;
auto x136 = x135 + x61 + x93 + V{1};
auto x137 = -x136*x81;
auto x138 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[9] + x134 + x137 + x43*x76 + x59*x81;
auto x139 = x123*x53;
auto x140 = x122 + x135 + x64;
auto x141 = -x140*x53;
auto x142 = x133*x53;
auto x143 = x21*(V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{3.33066907387547e-16}*cell[4] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x115 + x124 + x43*x97 + x44*x76 + x45*x96 + x52 + x53*x80 + x60 + x68 + x75 + x81*x87 + x81*x91 + x92*x95) + x25*(V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[17] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + x104 + x107 + x110 + x111 + x124 + x125 + x127 + x131 + x133*x81 + x138 + x44*x97 + x45*(V{3.23815048849004e-17}*x21*x38 + V{3.23815048849004e-17}*x25*x40 + V{3.23815048849004e-17}*x31*x41) + x51 + x53*x91 + x80*x81 + x92*x99 + V{-2.22044604925031e-16}) + x31*(V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x108 + x115 + x117*x92 + x130*x81 + x138 + x139 + x141 + x142 + x44*x96 + x45*x97 + x53*x87 + x74*x81);
auto x144 = V{0.0555555555555556}*x21*x38 + V{0.0555555555555556}*x25*x40 + V{0.0555555555555556}*x31*x41;
auto x145 = V{3.46944695195361e-18}*cell[0];
auto x146 = V{3.46944695195361e-18}*cell[12] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[4] + x145 + V{3.46944695195361e-18};
auto x147 = V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[6];
auto x148 = x21*(V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{6.93889390390723e-18}*cell[5] + V{6.93889390390723e-18}*cell[7] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x146 + x147);
auto x149 = V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[8];
auto x150 = x25*(V{6.93889390390723e-18}*cell[14] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + V{6.93889390390723e-18}*cell[9] + x146 + x149);
auto x151 = x31*(V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[16] + V{6.93889390390723e-18}*cell[18] + V{6.93889390390723e-18}*cell[3] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + x145 + x147 + x149 + V{3.46944695195361e-18});
auto x152 = V{0.0277777777777778}*x21*x38 + V{0.0277777777777778}*x25*x40 + V{0.0277777777777778}*x31*x41;
auto x153 = -V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] - V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x148 + x150 + x151 - x152*x45 + V{-0.0555555555555555};
auto x154 = -V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] - V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] - x152*x44;
auto x155 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - V{0.166666666666667}*cell[9] + x144*x43 + x153 + x154;
auto x156 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] - x152*x43;
auto x157 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + x144*x44 + x153 + x156 - x26 - x27 - x28 - x29;
auto x158 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + x144*x45 + x148 + x150 + x151 + x154 + x156 - x32 - x33 - x34 - x35 + V{-0.0555555555555555};
auto x159 = x21*x38 + x25*x40 + x31*x41;
auto x160 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x159;
auto x161 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x160;
auto x162 = V{0.0138888888888889}*x21*x38 + V{0.0138888888888889}*x25*x40 + V{0.0138888888888889}*x31*x41;
auto x163 = V{0.0138888888888889}*cell[0];
auto x164 = V{0.0138888888888889}*cell[12] + V{0.0138888888888889}*cell[3] + V{0.0277777777777778}*cell[4] + x163 + V{0.0138888888888889};
auto x165 = V{0.0138888888888889}*cell[11] + V{0.0138888888888889}*cell[2] + V{0.0277777777777778}*cell[6];
auto x166 = x21*(V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[7] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + x164 + x165);
auto x167 = V{0.0138888888888889}*cell[10] + V{0.0138888888888889}*cell[1] + V{0.0277777777777778}*cell[8];
auto x168 = x25*(V{0.0277777777777778}*cell[14] + V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0277777777777778}*cell[2] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + V{0.0277777777777778}*cell[9] + x164 + x167);
auto x169 = x31*(V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[3] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + x163 + x165 + x167 + V{0.0138888888888889});
auto x170 = V{0.0277777777777778}*x21*x38 + V{0.0277777777777778}*x25*x40 + V{0.0277777777777778}*x31*x41;
auto x171 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x166 + x168 + x169 + x170*x43 + V{0.0138888888888889};
auto x172 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x170*x44;
auto x173 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x162*x45 + x171 + x172;
auto x174 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x161 + x173;
auto x175 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x161 + x173;
auto x176 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x160;
auto x177 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + x170*x45;
auto x178 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x162*x44 + x171 + x177;
auto x179 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x176 + x178;
auto x180 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x176 + x178;
auto x181 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x159;
auto x182 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x162*x43 + x166 + x168 + x169 + x172 + x177 + V{0.0138888888888889};
auto x183 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] + x181 + x182;
auto x184 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] - x181 + x182;
auto x0 = -x19*(-V{0.111111111111111}*x143*x50 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[7] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x23 + x24) + x25*(V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[9] + x23 + x26 + x27 + x28 + x29 + x30) + x31*(V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[3] + x22 + x24 + x30 + x32 + x33 + x34 + x35 + V{0.166666666666667}) + x42*x43 + x42*x44 + x42*x45 - x51 + V{0.833333333333333});
auto x1 = -x19*(-V{0.0185185185185185}*x143*x95 + V{0.0555555555555556}) + x20*(-x134 - x155);
auto x2 = -x19*(-V{0.0185185185185185}*x143*x99 + V{0.0555555555555556}) - x20*(x100 + x157);
auto x3 = -x19*(-V{0.0185185185185185}*x117*x143 + V{0.0555555555555556}) - x20*(x118 + x158);
auto x4 = -x19*(-V{0.00925925925925926}*x123*x143 + V{0.0277777777777778}) - x20*(x139 + x174);
auto x5 = -x19*(-V{0.00925925925925926}*x143*x87 + V{0.0277777777777778}) - x20*(x175 + x53*(x132 + x86));
auto x6 = -x19*(-V{0.00925925925925926}*x114*x143 + V{0.0277777777777778}) - x20*(x125 + x179);
auto x7 = -x19*(-V{0.00925925925925926}*x143*x91 + V{0.0277777777777778}) - x20*(x180 + x53*(x129 + x90));
auto x8 = -x19*(-V{0.00925925925925926}*x143*x59 + V{0.0277777777777778}) - x20*(x183 + x60);
auto x9 = -x19*(-V{0.00925925925925926}*x143*x80 + V{0.0277777777777778}) - x20*(x184 + x53*(x72 + x79));
auto x10 = -x19*(V{0.0185185185185185}*x136*x143 + V{0.0555555555555556}) + x20*(-x137 - x155);
auto x11 = -x19*(V{0.0185185185185185}*x101*x143 + V{0.0555555555555556}) - x20*(x102 + x157);
auto x12 = -x19*(V{0.0185185185185185}*x119*x143 + V{0.0555555555555556}) - x20*(x120 + x158);
auto x13 = -x19*(V{0.00925925925925926}*x140*x143 + V{0.0277777777777778}) - x20*(x141 + x174);
auto x14 = -x19*(-V{0.00925925925925926}*x133*x143 + V{0.0277777777777778}) - x20*(x142 + x175);
auto x15 = -x19*(V{0.00925925925925926}*x126*x143 + V{0.0277777777777778}) - x20*(x127 + x179);
auto x16 = -x19*(-V{0.00925925925925926}*x130*x143 + V{0.0277777777777778}) - x20*(x131 + x180);
auto x17 = -x19*(V{0.00925925925925926}*x143*x67 + V{0.0277777777777778}) - x20*(x183 + x68);
auto x18 = -x19*(-V{0.00925925925925926}*x143*x74 + V{0.0277777777777778}) - x20*(x184 + x75);
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
return { -V{0.333333333333333}*x143, x43 + x44 + x45 };
}
};

}

}

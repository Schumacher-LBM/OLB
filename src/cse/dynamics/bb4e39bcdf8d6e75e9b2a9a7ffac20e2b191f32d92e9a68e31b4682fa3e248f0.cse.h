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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<1, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<1, -1, -1>, momenta::DefineSeparately> >> {
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
auto x27 = V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] + V{0.333333333333333}*cell[8] + x26 + V{0.166666666666667};
auto x28 = V{0.166666666666667}*cell[12] + V{0.333333333333333}*cell[14] + V{0.166666666666667}*cell[3];
auto x29 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1});
auto x30 = V{0.166666666666667}*cell[13];
auto x31 = V{0.166666666666667}*cell[14];
auto x32 = V{0.166666666666667}*cell[4];
auto x33 = V{0.166666666666667}*cell[5];
auto x34 = V{0.166666666666667}*cell[11] + V{0.333333333333333}*cell[16] + V{0.166666666666667}*cell[2];
auto x35 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x36 = cell[0] + cell[12] + V{2}*cell[14] + cell[3] + V{1};
auto x37 = cell[11] + V{2}*cell[16] + cell[2];
auto x38 = V{2}*cell[10] + V{2}*cell[13] + V{2}*cell[15] + cell[17] + cell[18] + cell[8] + cell[9] + x36 + x37;
auto x39 = cell[10] + cell[1] + V{2}*cell[8];
auto x40 = cell[15] + cell[16] + V{2}*cell[2] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[9] + x36 + x39;
auto x41 = cell[0] + cell[13] + cell[14] + V{2}*cell[18] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + x37 + x39 + V{1};
auto x42 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x43 = V{1.5}*x42;
auto x44 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x45 = V{1.5}*x44;
auto x46 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x47 = V{1.5}*x46;
auto x48 = x45 + x47 + V{-1};
auto x49 = x43 + x48;
auto x50 = x49*(-V{0.111111111111111}*x21*x40 - V{0.111111111111111}*x29*x41 + V{0.111111111111111}*x35*x38);
auto x51 = -V{0.166666666666667}*x21*x40 - V{0.166666666666667}*x29*x41 + V{0.166666666666667}*x35*x38;
auto x52 = V{1.11022302462516e-16}*cell[4];
auto x53 = V{1.11022302462516e-16}*cell[13];
auto x54 = V{1.11022302462516e-16}*cell[2];
auto x55 = V{5.55111512312578e-17}*cell[0];
auto x56 = V{5.55111512312578e-17}*cell[10] + V{5.55111512312578e-17}*cell[1] + V{1.11022302462516e-16}*cell[8] + x55 + V{5.55111512312578e-17};
auto x57 = V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[14] + V{5.55111512312578e-17}*cell[3];
auto x58 = x21*(V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + V{1.11022302462516e-16}*cell[9] + x52 + x54 + x56 + x57);
auto x59 = V{1.11022302462516e-16}*cell[3];
auto x60 = V{5.55111512312578e-17}*cell[11] + V{1.11022302462516e-16}*cell[16] + V{5.55111512312578e-17}*cell[2];
auto x61 = x29*(V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[6] + x56 + x59 + x60);
auto x62 = -x35*(V{1.11022302462516e-16}*cell[10] + V{1.11022302462516e-16}*cell[15] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x53 + x55 + x57 + x60 + V{5.55111512312578e-17});
auto x63 = -x50;
auto x64 = -V{0.00925925925925926}*x21*x40 - V{0.00925925925925926}*x29*x41 + V{0.00925925925925926}*x35*x38;
auto x65 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x66 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x67 = V{4.5}*(x66*x66);
auto x68 = -x43;
auto x69 = V{1} - x45;
auto x70 = x68 + x69;
auto x71 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x72 = -x47;
auto x73 = x71 + x72;
auto x74 = x65 + x67 + x70 + x73;
auto x75 = x49 + x71;
auto x76 = x65 - x67 + x75;
auto x77 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x78 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x77;
auto x79 = -x78;
auto x80 = -x71;
auto x81 = x49 + x65;
auto x82 = x80 + x81;
auto x83 = x82 - V{4.5}*x79*x79;
auto x84 = -x65;
auto x85 = -V{4.5}*x78*x78;
auto x86 = x75 + x84 + x85;
auto x87 = -V{0.0185185185185185}*x21*x40 - V{0.0185185185185185}*x29*x41 + V{0.0185185185185185}*x35*x38;
auto x88 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x89 = V{4.5}*(x88*x88);
auto x90 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x91 = x49 + x90;
auto x92 = x65 - x89 + x91;
auto x93 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x77;
auto x94 = -x93;
auto x95 = x80 + x91;
auto x96 = x95 - V{4.5}*x94*x94;
auto x97 = -V{0.037037037037037}*x21*x40 - V{0.037037037037037}*x29*x41 + V{0.037037037037037}*x35*x38;
auto x98 = V{3}*x46;
auto x99 = x43 + V{-1};
auto x100 = x45 + x90 - x98 + x99;
auto x101 = -V{7.40148683083438e-17}*x21*x40 - V{7.40148683083438e-17}*x29*x41 + V{7.40148683083438e-17}*x35*x38;
auto x102 = V{3}*x44;
auto x103 = -x102 + x47 + x65 + x99;
auto x104 = -x103*x87;
auto x105 = x65 + x72;
auto x106 = x102 + x105 + x68 + V{1};
auto x107 = x106*x87;
auto x108 = -V{5.55111512312578e-17}*x21*x40 - V{5.55111512312578e-17}*x29*x41 + V{5.55111512312578e-17}*x35*x38;
auto x109 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x110 = V{4.5}*(x109*x109);
auto x111 = -x110 + x71 + x91;
auto x112 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[9] + x104 + x107 - x108*x44 - x111*x87;
auto x113 = V{3}*x42;
auto x114 = x113 + x69 + x73;
auto x115 = -x113 + x48 + x71;
auto x116 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x117 = -V{4.5}*x116*x116;
auto x118 = x117 + x84 + x91;
auto x119 = V{3.33066907387547e-16}*cell[14] + V{3.33066907387547e-16}*cell[5] + x114*x87 - x115*x87 - x118*x87;
auto x120 = x70 + x90;
auto x121 = x105 + x120 + x89;
auto x122 = -x116;
auto x123 = -x90;
auto x124 = x123 + x81;
auto x125 = x124 - V{4.5}*x122*x122;
auto x126 = -V{4.5}*x93*x93;
auto x127 = x123 + x126 + x75;
auto x128 = -V{3.70074341541719e-17}*x21*x40 - V{3.70074341541719e-17}*x29*x41 + V{3.70074341541719e-17}*x35*x38;
auto x129 = x120 + x98;
auto x130 = -x100*x87 + x129*x87 + x58 + x61 + x62 + x63 - x86*x87 + V{-2.22044604925031e-16};
auto x131 = x110 + x120 + x73;
auto x132 = x21*(V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[17] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] - x100*x97 - x101*x46 + x112 + x119 - x42*(-V{3.23815048849004e-17}*x21*x40 - V{3.23815048849004e-17}*x29*x41 + V{3.23815048849004e-17}*x35*x38) + x52 + x53 + x58 + x61 + x62 + x63 + x64*x74 - x64*x76 - x64*x83 - x64*x86 - x87*x92 - x87*x96 + V{-2.22044604925031e-16}) + x29*(V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] - x101*x42 + x112 - x115*x97 - x118*x64 + x121*x64 - x125*x64 - x127*x87 - x128*x46 + x130 + x54 - x64*x92 - x76*x87) - x35*(V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{3.33066907387547e-16}*cell[4] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] - x101*x44 + x106*x97 - x108*x46 - x111*x64 + x119 + x121*x87 - x127*x64 - x128*x42 + x130 + x131*x64 + x59 - x64*x96 + x74*x87);
auto x133 = -x132;
auto x134 = -V{0.0277777777777778}*x21*x40 - V{0.0277777777777778}*x29*x41 + V{0.0277777777777778}*x35*x38;
auto x135 = x21*x40 + x29*x41 - x35*x38;
auto x136 = V{3.46944695195361e-18}*cell[0];
auto x137 = V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[14] + V{3.46944695195361e-18}*cell[3] + x136 + V{3.46944695195361e-18};
auto x138 = V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[8];
auto x139 = x21*(V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{6.93889390390723e-18}*cell[4] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + V{6.93889390390723e-18}*cell[9] + x137 + x138);
auto x140 = V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[2];
auto x141 = x29*(V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[18] + V{6.93889390390723e-18}*cell[3] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[6] + x136 + x138 + x140 + V{3.46944695195361e-18});
auto x142 = -x35*(V{6.93889390390723e-18}*cell[10] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[15] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x137 + x140);
auto x143 = -V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] - V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x139 + x141 + x142 + V{-0.0555555555555555};
auto x144 = -V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] - V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7];
auto x145 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - V{0.166666666666667}*cell[9] + x134*x42 + x134*x46 + V{0.0555555555555556}*x135*x44 + x143 + x144;
auto x146 = V{0.0185185185185185}*x21*x40 + V{0.0185185185185185}*x29*x41 - V{0.0185185185185185}*x35*x38;
auto x147 = V{0.0555555555555556}*x21*x40 + V{0.0555555555555556}*x29*x41 - V{0.0555555555555556}*x35*x38;
auto x148 = V{0.0277777777777778}*x21*x40 + V{0.0277777777777778}*x29*x41 - V{0.0277777777777778}*x35*x38;
auto x149 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] - x148*x44;
auto x150 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + x143 + x147*x46 - x148*x42 + x149 - x22 - x23 - x24 - x25;
auto x151 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + x139 + x141 + x142 + x144 + x147*x42 - x148*x46 + x149 - x30 - x31 - x32 - x33 + V{-0.0555555555555555};
auto x152 = V{0.00925925925925926}*x21*x40 + V{0.00925925925925926}*x29*x41 - V{0.00925925925925926}*x35*x38;
auto x153 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x135;
auto x154 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x153;
auto x155 = V{0.0138888888888889}*x21*x40 + V{0.0138888888888889}*x29*x41 - V{0.0138888888888889}*x35*x38;
auto x156 = V{0.0138888888888889}*cell[0];
auto x157 = V{0.0138888888888889}*cell[10] + V{0.0138888888888889}*cell[1] + V{0.0277777777777778}*cell[8] + x156 + V{0.0138888888888889};
auto x158 = V{0.0138888888888889}*cell[12] + V{0.0277777777777778}*cell[14] + V{0.0138888888888889}*cell[3];
auto x159 = x21*(V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[4] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + V{0.0277777777777778}*cell[9] + x157 + x158);
auto x160 = V{0.0138888888888889}*cell[11] + V{0.0277777777777778}*cell[16] + V{0.0138888888888889}*cell[2];
auto x161 = x29*(V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[3] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + V{0.0277777777777778}*cell[6] + x157 + x160);
auto x162 = -x35*(V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[15] + V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + x156 + x158 + x160 + V{0.0138888888888889});
auto x163 = V{0.0277777777777778}*x21*x40 + V{0.0277777777777778}*x29*x41 - V{0.0277777777777778}*x35*x38;
auto x164 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x159 + x161 + x162 + x163*x44 + V{0.0138888888888889};
auto x165 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x163*x46;
auto x166 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x155*x42 + x164 + x165;
auto x167 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x154 + x166;
auto x168 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x154 + x166;
auto x169 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x153;
auto x170 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + x163*x42;
auto x171 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x155*x46 + x164 + x170;
auto x172 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x169 + x171;
auto x173 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x169 + x171;
auto x174 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x135;
auto x175 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x155*x44 + x159 + x161 + x162 + x165 + x170 + V{0.0138888888888889};
auto x176 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] + x174 + x175;
auto x177 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] - x174 + x175;
auto x0 = -x19*(V{0.111111111111111}*x133*x49 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21*(V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[9] + x22 + x23 + x24 + x25 + x27 + x28) + x29*(V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[6] + x27 + x30 + x31 + x32 + x33 + x34) - x35*(V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[15] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x26 + x28 + x34 + V{0.166666666666667}) - x42*x51 - x44*x51 - x46*x51 + x50 + V{0.833333333333333});
auto x1 = -x19*(V{0.0185185185185185}*x103*x133 + V{0.0555555555555556}) + x20*(-x104 - x145);
auto x2 = -x19*(V{0.0185185185185185}*x100*x133 + V{0.0555555555555556}) - x20*(x100*x146 + x150);
auto x3 = -x19*(V{0.0185185185185185}*x115*x133 + V{0.0555555555555556}) - x20*(x115*x146 + x151);
auto x4 = -x19*(V{0.00925925925925926}*x133*x92 + V{0.0277777777777778}) - x20*(x152*x92 + x167);
auto x5 = -x19*(V{0.00925925925925926}*x125*x133 + V{0.0277777777777778}) - x20*(x152*(x117 + x124) + x168);
auto x6 = -x19*(V{0.00925925925925926}*x133*x76 + V{0.0277777777777778}) - x20*(x152*x76 + x172);
auto x7 = -x19*(V{0.00925925925925926}*x133*x83 + V{0.0277777777777778}) - x20*(x152*(x82 + x85) + x173);
auto x8 = -x19*(V{0.00925925925925926}*x111*x133 + V{0.0277777777777778}) - x20*(x111*x152 + x176);
auto x9 = -x19*(V{0.00925925925925926}*x133*x96 + V{0.0277777777777778}) - x20*(x152*(x126 + x95) + x177);
auto x10 = -x19*(V{0.0185185185185185}*x106*x132 + V{0.0555555555555556}) + x20*(-x107 - x145);
auto x11 = -x19*(V{0.0185185185185185}*x129*x132 + V{0.0555555555555556}) - x20*(-x129*x146 + x150);
auto x12 = -x19*(V{0.0185185185185185}*x114*x132 + V{0.0555555555555556}) - x20*(-x114*x146 + x151);
auto x13 = -x19*(V{0.00925925925925926}*x121*x132 + V{0.0277777777777778}) - x20*(-x121*x152 + x167);
auto x14 = -x19*(V{0.00925925925925926}*x118*x133 + V{0.0277777777777778}) - x20*(x118*x152 + x168);
auto x15 = -x19*(V{0.00925925925925926}*x132*x74 + V{0.0277777777777778}) - x20*(-x152*x74 + x172);
auto x16 = -x19*(V{0.00925925925925926}*x133*x86 + V{0.0277777777777778}) - x20*(x152*x86 + x173);
auto x17 = -x19*(V{0.00925925925925926}*x131*x132 + V{0.0277777777777778}) - x20*(-x131*x152 + x176);
auto x18 = -x19*(V{0.00925925925925926}*x127*x133 + V{0.0277777777777778}) - x20*(x127*x152 + x177);
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
return { V{0.333333333333333}*x133, x42 + x44 + x46 };
}
};

}

}

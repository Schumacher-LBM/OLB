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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<-1, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<-1, -1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1});
auto x22 = V{0.166666666666667}*cell[0];
auto x23 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + V{0.333333333333333}*cell[4] + x22 + V{0.166666666666667};
auto x24 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + V{0.333333333333333}*cell[7];
auto x25 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1});
auto x26 = V{0.166666666666667}*cell[15];
auto x27 = V{0.166666666666667}*cell[16];
auto x28 = V{0.166666666666667}*cell[6];
auto x29 = V{0.166666666666667}*cell[7];
auto x30 = V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[1] + V{0.333333333333333}*cell[9];
auto x31 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x32 = V{0.166666666666667}*cell[13];
auto x33 = V{0.166666666666667}*cell[14];
auto x34 = V{0.166666666666667}*cell[4];
auto x35 = V{0.166666666666667}*cell[5];
auto x36 = cell[0] + cell[11] + cell[2] + V{2}*cell[7] + V{1};
auto x37 = cell[10] + cell[1] + V{2}*cell[9];
auto x38 = V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + V{2}*cell[17] + cell[4] + cell[5] + x36 + x37;
auto x39 = cell[12] + cell[3] + V{2}*cell[4];
auto x40 = cell[17] + cell[18] + V{2}*cell[1] + V{2}*cell[5] + V{2}*cell[6] + cell[8] + cell[9] + x36 + x39;
auto x41 = cell[0] + V{2}*cell[14] + cell[15] + cell[16] + V{2}*cell[2] + cell[6] + cell[7] + V{2}*cell[8] + x37 + x39 + V{1};
auto x42 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x43 = V{1.5}*x42;
auto x44 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x45 = V{1.5}*x44;
auto x46 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x47 = V{1.5}*x46;
auto x48 = x45 + x47 + V{-1};
auto x49 = x43 + x48;
auto x50 = x49*(-V{0.111111111111111}*x21*x40 - V{0.111111111111111}*x25*x41 + V{0.111111111111111}*x31*x38);
auto x51 = -V{0.166666666666667}*x21*x40 - V{0.166666666666667}*x25*x41 + V{0.166666666666667}*x31*x38;
auto x52 = V{1.11022302462516e-16}*cell[12];
auto x53 = -V{0.00925925925925926}*x21*x40 - V{0.00925925925925926}*x25*x41 + V{0.00925925925925926}*x31*x38;
auto x54 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x55 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x56 = V{4.5}*(x55*x55);
auto x57 = -x43;
auto x58 = V{1} - x47;
auto x59 = x57 + x58;
auto x60 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x61 = -x45;
auto x62 = x60 + x61;
auto x63 = x54 + x56 + x59 + x62;
auto x64 = -V{5.55111512312578e-17}*x21*x40 - V{5.55111512312578e-17}*x25*x41 + V{5.55111512312578e-17}*x31*x38;
auto x65 = x49 + x60;
auto x66 = x54 - x56 + x65;
auto x67 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x68 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x67;
auto x69 = -x68;
auto x70 = -x54;
auto x71 = x65 + x70;
auto x72 = x71 - V{4.5}*x69*x69;
auto x73 = -x60;
auto x74 = -V{4.5}*x68*x68;
auto x75 = x49 + x54;
auto x76 = x73 + x74 + x75;
auto x77 = -V{0.0185185185185185}*x21*x40 - V{0.0185185185185185}*x25*x41 + V{0.0185185185185185}*x31*x38;
auto x78 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x79 = V{4.5}*(x78*x78);
auto x80 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x81 = x49 + x80;
auto x82 = x54 - x79 + x81;
auto x83 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x84 = -x83;
auto x85 = x73 + x81;
auto x86 = x85 - V{4.5}*x84*x84;
auto x87 = -V{0.037037037037037}*x21*x40 - V{0.037037037037037}*x25*x41 + V{0.037037037037037}*x31*x38;
auto x88 = V{3}*x44;
auto x89 = x43 + V{-1};
auto x90 = x47 + x80 - x88 + x89;
auto x91 = -V{3.70074341541719e-17}*x21*x40 - V{3.70074341541719e-17}*x25*x41 + V{3.70074341541719e-17}*x31*x38;
auto x92 = -V{7.40148683083438e-17}*x21*x40 - V{7.40148683083438e-17}*x25*x41 + V{7.40148683083438e-17}*x31*x38;
auto x93 = V{5.55111512312578e-17}*cell[0];
auto x94 = V{1.11022302462516e-16}*cell[4];
auto x95 = V{5.55111512312578e-17}*cell[12] + V{5.55111512312578e-17}*cell[3] + x93 + x94 + V{5.55111512312578e-17};
auto x96 = V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[2] + V{1.11022302462516e-16}*cell[7];
auto x97 = x21*(V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{1.11022302462516e-16}*cell[1] + V{1.11022302462516e-16}*cell[5] + V{1.11022302462516e-16}*cell[6] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x95 + x96);
auto x98 = V{1.11022302462516e-16}*cell[2];
auto x99 = V{5.55111512312578e-17}*cell[10] + V{5.55111512312578e-17}*cell[1] + V{1.11022302462516e-16}*cell[9];
auto x100 = x25*(V{1.11022302462516e-16}*cell[14] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + V{1.11022302462516e-16}*cell[8] + x95 + x98 + x99);
auto x101 = -x31*(V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[15] + V{1.11022302462516e-16}*cell[17] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + x52 + x93 + x96 + x99 + V{5.55111512312578e-17});
auto x102 = V{3}*x42;
auto x103 = x102 + x58 + x62;
auto x104 = -x102 + x48 + x60;
auto x105 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x67;
auto x106 = -x105;
auto x107 = x70 + x81;
auto x108 = x107 - V{4.5}*x106*x106;
auto x109 = -x50;
auto x110 = x100 + x101 + x103*x77 - x104*x77 - x108*x77 + x109 + x97 + V{-2.22044604925031e-16};
auto x111 = V{3}*x46;
auto x112 = x54 + x61;
auto x113 = x111 + x112 + x57 + V{1};
auto x114 = -x111 + x45 + x54 + x89;
auto x115 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x116 = V{4.5}*(x115*x115);
auto x117 = -x116 + x60 + x81;
auto x118 = V{3.33066907387547e-16}*cell[14] + V{3.33066907387547e-16}*cell[5] + x113*x77 - x114*x77 - x117*x77;
auto x119 = x59 + x80;
auto x120 = x112 + x119 + x79;
auto x121 = -x80;
auto x122 = -V{4.5}*x105*x105;
auto x123 = x121 + x122 + x75;
auto x124 = -V{4.5}*x83*x83;
auto x125 = x121 + x124 + x65;
auto x126 = -x77*x90;
auto x127 = x119 + x88;
auto x128 = x127*x77;
auto x129 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[17] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[8] + x126 + x128 - x44*x64 - x72*x77;
auto x130 = x116 + x119 + x62;
auto x131 = x21*(V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{3.33066907387547e-16}*cell[13] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[4] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x110 + x118 - x42*x64 - x44*x92 - x46*x91 + x52 + x53*x63 - x53*x66 - x53*x72 - x53*x76 - x77*x82 - x77*x86 - x87*x90) + x25*(V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[18] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{2.22044604925031e-16}*cell[9] + x100 + x101 - x104*x87 - x108*x53 + x109 + x118 + x120*x53 - x123*x53 - x125*x77 + x129 - x42*x92 - x46*(-V{3.23815048849004e-17}*x21*x40 - V{3.23815048849004e-17}*x25*x41 + V{3.23815048849004e-17}*x31*x38) - x53*x82 - x66*x77 + x94 + x97 + V{-2.22044604925031e-16}) - x31*(V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x110 + x113*x87 - x117*x53 + x120*x77 - x125*x53 + x129 + x130*x53 - x42*x91 - x46*x92 - x53*x86 + x63*x77 + x98);
auto x132 = -x131;
auto x133 = -V{0.0277777777777778}*x21*x40 - V{0.0277777777777778}*x25*x41 + V{0.0277777777777778}*x31*x38;
auto x134 = x21*x40 + x25*x41 - x31*x38;
auto x135 = V{3.46944695195361e-18}*cell[0];
auto x136 = V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[7] + x135 + V{3.46944695195361e-18};
auto x137 = V{3.46944695195361e-18}*cell[12] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[4];
auto x138 = x21*(V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{6.93889390390723e-18}*cell[5] + V{6.93889390390723e-18}*cell[6] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x136 + x137);
auto x139 = V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[9];
auto x140 = x25*(V{6.93889390390723e-18}*cell[14] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + V{6.93889390390723e-18}*cell[8] + x135 + x137 + x139 + V{3.46944695195361e-18});
auto x141 = -x31*(V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + x136 + x139);
auto x142 = -V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] - V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x138 + x140 + x141 + V{-0.0555555555555555};
auto x143 = -V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] - V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7];
auto x144 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - V{0.166666666666667}*cell[9] + x133*x42 + x133*x46 + V{0.0555555555555556}*x134*x44 + x142 + x143;
auto x145 = V{0.0185185185185185}*x21*x40 + V{0.0185185185185185}*x25*x41 - V{0.0185185185185185}*x31*x38;
auto x146 = V{0.0555555555555556}*x21*x40 + V{0.0555555555555556}*x25*x41 - V{0.0555555555555556}*x31*x38;
auto x147 = V{0.0277777777777778}*x21*x40 + V{0.0277777777777778}*x25*x41 - V{0.0277777777777778}*x31*x38;
auto x148 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] - x147*x44;
auto x149 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + x142 + x146*x42 - x147*x46 + x148 - x26 - x27 - x28 - x29;
auto x150 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + x138 + x140 + x141 + x143 + x146*x46 - x147*x42 + x148 - x32 - x33 - x34 - x35 + V{-0.0555555555555555};
auto x151 = V{0.00925925925925926}*x21*x40 + V{0.00925925925925926}*x25*x41 - V{0.00925925925925926}*x31*x38;
auto x152 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x134;
auto x153 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x152;
auto x154 = V{0.0138888888888889}*x21*x40 + V{0.0138888888888889}*x25*x41 - V{0.0138888888888889}*x31*x38;
auto x155 = V{0.0138888888888889}*cell[0];
auto x156 = V{0.0138888888888889}*cell[12] + V{0.0138888888888889}*cell[3] + V{0.0277777777777778}*cell[4] + x155 + V{0.0138888888888889};
auto x157 = V{0.0138888888888889}*cell[11] + V{0.0138888888888889}*cell[2] + V{0.0277777777777778}*cell[7];
auto x158 = x21*(V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + x156 + x157);
auto x159 = V{0.0138888888888889}*cell[10] + V{0.0138888888888889}*cell[1] + V{0.0277777777777778}*cell[9];
auto x160 = x25*(V{0.0277777777777778}*cell[14] + V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0277777777777778}*cell[2] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + V{0.0277777777777778}*cell[8] + x156 + x159);
auto x161 = -x31*(V{0.0277777777777778}*cell[12] + V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[17] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + x155 + x157 + x159 + V{0.0138888888888889});
auto x162 = V{0.0277777777777778}*x21*x40 + V{0.0277777777777778}*x25*x41 - V{0.0277777777777778}*x31*x38;
auto x163 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x158 + x160 + x161 + x162*x44 + V{0.0138888888888889};
auto x164 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x162*x42;
auto x165 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x154*x46 + x163 + x164;
auto x166 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x153 + x165;
auto x167 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x153 + x165;
auto x168 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x152;
auto x169 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + x162*x46;
auto x170 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x154*x42 + x163 + x169;
auto x171 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x168 + x170;
auto x172 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x168 + x170;
auto x173 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x134;
auto x174 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x154*x44 + x158 + x160 + x161 + x164 + x169 + V{0.0138888888888889};
auto x175 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] + x173 + x174;
auto x176 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] - x173 + x174;
auto x0 = -x19*(V{0.111111111111111}*x132*x49 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x23 + x24) + x25*(V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[8] + x23 + x26 + x27 + x28 + x29 + x30) - x31*(V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[17] + x22 + x24 + x30 + x32 + x33 + x34 + x35 + V{0.166666666666667}) - x42*x51 - x44*x51 - x46*x51 + x50 + V{0.833333333333333});
auto x1 = -x19*(V{0.0185185185185185}*x132*x90 + V{0.0555555555555556}) + x20*(-x126 - x144);
auto x2 = -x19*(V{0.0185185185185185}*x104*x132 + V{0.0555555555555556}) - x20*(x104*x145 + x149);
auto x3 = -x19*(V{0.0185185185185185}*x114*x132 + V{0.0555555555555556}) - x20*(x114*x145 + x150);
auto x4 = -x19*(V{0.00925925925925926}*x117*x132 + V{0.0277777777777778}) - x20*(x117*x151 + x166);
auto x5 = -x19*(V{0.00925925925925926}*x132*x86 + V{0.0277777777777778}) - x20*(x151*(x124 + x85) + x167);
auto x6 = -x19*(V{0.00925925925925926}*x132*x82 + V{0.0277777777777778}) - x20*(x151*x82 + x171);
auto x7 = -x19*(V{0.00925925925925926}*x108*x132 + V{0.0277777777777778}) - x20*(x151*(x107 + x122) + x172);
auto x8 = -x19*(V{0.00925925925925926}*x132*x66 + V{0.0277777777777778}) - x20*(x151*x66 + x175);
auto x9 = -x19*(V{0.00925925925925926}*x132*x72 + V{0.0277777777777778}) - x20*(x151*(x71 + x74) + x176);
auto x10 = -x19*(V{0.0185185185185185}*x127*x131 + V{0.0555555555555556}) + x20*(-x128 - x144);
auto x11 = -x19*(V{0.0185185185185185}*x103*x131 + V{0.0555555555555556}) - x20*(-x103*x145 + x149);
auto x12 = -x19*(V{0.0185185185185185}*x113*x131 + V{0.0555555555555556}) - x20*(-x113*x145 + x150);
auto x13 = -x19*(V{0.00925925925925926}*x130*x131 + V{0.0277777777777778}) - x20*(-x130*x151 + x166);
auto x14 = -x19*(V{0.00925925925925926}*x125*x132 + V{0.0277777777777778}) - x20*(x125*x151 + x167);
auto x15 = -x19*(V{0.00925925925925926}*x120*x131 + V{0.0277777777777778}) - x20*(-x120*x151 + x171);
auto x16 = -x19*(V{0.00925925925925926}*x123*x132 + V{0.0277777777777778}) - x20*(x123*x151 + x172);
auto x17 = -x19*(V{0.00925925925925926}*x131*x63 + V{0.0277777777777778}) - x20*(-x151*x63 + x175);
auto x18 = -x19*(V{0.00925925925925926}*x132*x76 + V{0.0277777777777778}) - x20*(x151*x76 + x176);
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
return { V{0.333333333333333}*x132, x42 + x44 + x46 };
}
};

}

}

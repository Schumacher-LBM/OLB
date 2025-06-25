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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, -1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{0.5}*cell[12];
auto x22 = V{0.5}*cell[2];
auto x23 = V{2}*cell[9];
auto x24 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1});
auto x25 = V{0.25}*cell[0] + V{0.25}*cell[10] + V{0.25}*cell[1] + V{0.5}*cell[9] + V{0.25};
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{1});
auto x27 = cell[0] + cell[10] + cell[1] + x23 + V{1};
auto x28 = cell[11] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + V{2}*cell[17] + cell[2] + cell[4] + cell[5] + V{2}*cell[7] + x27;
auto x29 = cell[12] + V{2}*cell[14] + cell[15] + cell[16] + V{2}*cell[2] + cell[3] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[8] + x27;
auto x30 = -V{0.0138888888888889}*x24*x29 + V{0.0138888888888889}*x26*x28;
auto x31 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x31;
auto x33 = -x32;
auto x34 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x35 = -x34;
auto x36 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
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
auto x55 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x56 = x51 + x55;
auto x57 = x56 - V{4.5}*x54*x54;
auto x58 = -V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28;
auto x59 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x31;
auto x60 = -x59;
auto x61 = x44 + x55;
auto x62 = x35 + x61;
auto x63 = x62 - V{4.5}*x60*x60;
auto x64 = -V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28;
auto x65 = V{1.11022302462516e-16}*cell[2];
auto x66 = V{3}*x37;
auto x67 = x36 + x43 - x66;
auto x68 = -x58*x67;
auto x69 = V{1} - x42;
auto x70 = -x40;
auto x71 = x36 + x70;
auto x72 = x66 + x69 + x71;
auto x73 = x58*x72;
auto x74 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x75 = V{4.5}*(x74*x74);
auto x76 = x45 + x55 - x75;
auto x77 = -x30*x76;
auto x78 = -x38;
auto x79 = x69 + x78;
auto x80 = x55 + x79;
auto x81 = x71 + x75 + x80;
auto x82 = x30*x81;
auto x83 = -x55;
auto x84 = -V{4.5}*x53*x53;
auto x85 = x45 + x83 + x84;
auto x86 = -x30*x85;
auto x87 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x88 = V{4.5}*(x87*x87);
auto x89 = x34 + x70;
auto x90 = x80 + x88 + x89;
auto x91 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x92 = V{4.5}*(x91*x91);
auto x93 = x34 + x71 + x79 + x92;
auto x94 = -V{0.0555555555555556}*x24*x29 + V{0.0555555555555556}*x26*x28;
auto x95 = V{3}*x41;
auto x96 = x78 + x89 + x95 + V{1};
auto x97 = -V{1.11022302462516e-16}*x24*x29 + V{1.11022302462516e-16}*x26*x28;
auto x98 = V{3}*x39;
auto x99 = x38 + V{-1};
auto x100 = x42 + x55 - x98 + x99;
auto x101 = -x100*x58;
auto x102 = x80 + x98;
auto x103 = x102*x58;
auto x104 = V{1.11022302462516e-16}*cell[4];
auto x105 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[10] + V{5.55111512312578e-17}*cell[1] + V{1.11022302462516e-16}*cell[9] + V{5.55111512312578e-17};
auto x106 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[17] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[8] + x101 + x103 + x24*(V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[14] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{5.55111512312578e-17}*cell[3] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + V{1.11022302462516e-16}*cell[8] + x104 + x105 + x65) - x26*(V{5.55111512312578e-17}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[15] + V{1.11022302462516e-16}*cell[17] + V{5.55111512312578e-17}*cell[2] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[7] + x105) - x39*(-V{8.32667268468867e-17}*x24*x29 + V{8.32667268468867e-17}*x26*x28) - x48 + V{-2.22044604925031e-16};
auto x107 = V{1.11022302462516e-16}*cell[11] + V{2.22044604925031e-16}*cell[12] + V{4.44089209850063e-16}*cell[13] + V{4.9960036108132e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{2.22044604925031e-16}*cell[3] + V{4.44089209850063e-16}*cell[4] + V{4.9960036108132e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x106 - x37*(-V{5.55111512312578e-17}*x24*x29 + V{5.55111512312578e-17}*x26*x28) - x41*x97 + x58*x90 + x58*x93 + x65 + x68 + x73 + x77 + x82 + x86 + x94*x96;
auto x108 = -V{0.0277777777777778}*x24*x29 + V{0.0277777777777778}*x26*x28;
auto x109 = x34 + x40 - x95 + x99;
auto x110 = -x109*x58;
auto x111 = x58*x96;
auto x112 = x34 + x61 - x88;
auto x113 = -x112*x30;
auto x114 = x30*x90;
auto x115 = -V{4.5}*x59*x59;
auto x116 = x115 + x34 + x44 + x83;
auto x117 = -x116*x30;
auto x118 = x34 + x45 - x92;
auto x119 = V{2.22044604925031e-16}*cell[11] + V{9.71445146547012e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{2.22044604925031e-16}*cell[2] + V{9.71445146547012e-17}*cell[3] + V{3.33066907387547e-16}*cell[5] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[9] + x104 + x106 + x110 + x111 + x113 + x114 + x117 - x118*x58 - x37*x97 - x41*(-V{4.85722573273506e-17}*x24*x29 + V{4.85722573273506e-17}*x26*x28) - x52*(-V{3.08395284618099e-18}*x24*x29 + V{3.08395284618099e-18}*x26*x28) - x58*x76 - x58*x85 - x67*x94;
auto x120 = -x24*(-x108*x47 + x119 - x30*x63) + x26*(x107 - x30*x57 - x47*x64 - x58*x63);
auto x121 = -V{0.00231481481481481}*x24*x29 + V{0.00231481481481481}*x26*x28;
auto x122 = x46 + x50;
auto x123 = -V{0.0833333333333333}*x24*x29 + V{0.0833333333333333}*x26*x28;
auto x124 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[10] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[9] + V{3.46944695195361e-18};
auto x125 = x24*(V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[14] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[2] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + V{6.93889390390723e-18}*cell[8] + x124);
auto x126 = -x26*(V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[12] + V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[2] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[7] + x124);
auto x127 = -V{0.0416666666666667}*x24*x29 + V{0.0416666666666667}*x26*x28;
auto x128 = -V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] - V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x125 + x126 + x127*x41 + V{-0.0555555555555555};
auto x129 = -V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] - V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7] + x127*x37;
auto x130 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - V{0.333333333333333}*cell[9] - x121*x122 + x121*x52 - x123*x39 + x128 + x129;
auto x131 = -V{0.00115740740740741}*x24*x29 + V{0.00115740740740741}*x26*x28;
auto x132 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.166666666666667}*cell[9] + x122*x131 + x127*x39 - x131*x52;
auto x133 = V{0.166666666666667}*cell[11] - V{0.166666666666667}*cell[15] - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[2] - V{0.166666666666667}*cell[6] - V{0.166666666666667}*cell[7] - x123*x37 + x128 + x132;
auto x134 = V{0.166666666666667}*cell[12] - V{0.166666666666667}*cell[13] - V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[4] - V{0.166666666666667}*cell[5] - x123*x41 + x125 + x126 + x129 + x132 + V{-0.0555555555555555};
auto x135 = -x24*x29 + x26*x28;
auto x136 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x135;
auto x137 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x136;
auto x138 = -V{0.0208333333333333}*x24*x29 + V{0.0208333333333333}*x26*x28;
auto x139 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[10] + V{0.0208333333333333}*cell[1] + V{0.0416666666666667}*cell[9] + V{0.0208333333333333};
auto x140 = x24*(V{0.0208333333333333}*cell[12] + V{0.0416666666666667}*cell[14] + V{0.0208333333333333}*cell[15] + V{0.0208333333333333}*cell[16] + V{0.0416666666666667}*cell[2] + V{0.0208333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0208333333333333}*cell[6] + V{0.0208333333333333}*cell[7] + V{0.0416666666666667}*cell[8] + x139);
auto x141 = -x26*(V{0.0208333333333333}*cell[11] + V{0.0416666666666667}*cell[12] + V{0.0208333333333333}*cell[13] + V{0.0208333333333333}*cell[14] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[17] + V{0.0208333333333333}*cell[2] + V{0.0208333333333333}*cell[4] + V{0.0208333333333333}*cell[5] + V{0.0416666666666667}*cell[7] + x139);
auto x142 = -V{0.0416666666666667}*x24*x29 + V{0.0416666666666667}*x26*x28;
auto x143 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x140 + x141 - x142*x37 + V{0.0138888888888889};
auto x144 = -V{0.000578703703703704}*x24*x29 + V{0.000578703703703704}*x26*x28;
auto x145 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0833333333333334}*cell[9] + x122*x144 - x142*x39 - x144*x52;
auto x146 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x138*x41 + x143 + x145;
auto x147 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] - x137 + x146;
auto x148 = -x30*(x56 + x84);
auto x149 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] + x137 + x146;
auto x150 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x136;
auto x151 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] - x142*x41;
auto x152 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x138*x37 + x140 + x141 + x145 + x151 + V{0.0138888888888889};
auto x153 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x150 + x152;
auto x154 = x115 + x62;
auto x155 = -x154*x30;
auto x156 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x150 + x152;
auto x157 = -V{0.00115740740740741}*x24*x29 + V{0.00115740740740741}*x26*x28;
auto x158 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x135;
auto x159 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x138*x39 + x143 + x151;
auto x160 = V{0.416666666666667}*cell[17] + V{0.416666666666667}*cell[8] - V{0.166666666666667}*cell[9] - x122*x157 + x157*x52 - x158 + x159;
auto x161 = -V{0.00578703703703704}*x24*x29 + V{0.00578703703703704}*x26*x28;
auto x162 = -V{0.0833333333333333}*cell[17] - V{0.0833333333333333}*cell[8] + V{0.833333333333333}*cell[9] + x158 + x159;
auto x163 = -x24*(-x108*x122 + x119 + x155) + x26*(x107 - x122*x64 + x148 - x154*x58);
auto x0 = -x19*(V{0.166666666666667}*x120*x44 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{0.5}*cell[1] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + x21 + x22 + x23 + x24*(V{0.25}*cell[12] + V{0.5}*cell[14] + V{0.25}*cell[15] + V{0.25}*cell[16] + V{0.25}*cell[3] + V{0.5}*cell[4] + V{0.25}*cell[6] + V{0.25}*cell[7] + V{0.5}*cell[8] + x22 + x25) - x26*(V{0.25}*cell[11] + V{0.25}*cell[13] + V{0.25}*cell[14] + V{0.5}*cell[15] + V{0.5}*cell[17] + V{0.25}*cell[2] + V{0.25}*cell[4] + V{0.25}*cell[5] + V{0.5}*cell[7] + x21 + x25) + x30*x47 - x30*x52 - x37*x49 - x39*x49 - x41*x49 + x48 + V{0.833333333333333});
auto x1 = -x19*(V{0.0277777777777778}*x100*x120 + V{0.0555555555555556}) - x20*(x101 + x130);
auto x2 = -x19*(V{0.0277777777777778}*x120*x67 + V{0.0555555555555556}) - x20*(x133 + x68);
auto x3 = -x19*(V{0.0277777777777778}*x109*x120 + V{0.0555555555555556}) - x20*(x110 + x134);
auto x4 = -x19*(V{0.0138888888888889}*x120*x76 + V{0.0277777777777778}) - x20*(x147 + x77);
auto x5 = -x19*(V{0.0138888888888889}*x120*x57 + V{0.0277777777777778}) - x20*(x148 + x149);
auto x6 = -x19*(V{0.0138888888888889}*x112*x120 + V{0.0277777777777778}) - x20*(x113 + x153);
auto x7 = -x19*(V{0.0138888888888889}*x120*x63 + V{0.0277777777777778}) - x20*(x155 + x156);
auto x8 = -x19*(V{0.0138888888888889}*x118*x120 + V{0.0277777777777778}) - x20*(-x118*x30 + x160);
auto x9 = -x19*(V{0.0138888888888889}*x120*x47 + V{0.0277777777777778}) - x20*(-x122*(-V{0.00810185185185185}*x24*x29 + V{0.00810185185185185}*x26*x28) - x161*x52 + x162);
auto x10 = x19*(V{0.0277777777777778}*x102*x163 + V{-0.0555555555555556}) - x20*(x103 + x130);
auto x11 = x19*(V{0.0277777777777778}*x163*x72 + V{-0.0555555555555556}) - x20*(x133 + x73);
auto x12 = x19*(V{0.0277777777777778}*x163*x96 + V{-0.0555555555555556}) - x20*(x111 + x134);
auto x13 = x19*(V{0.0138888888888889}*x163*x81 + V{-0.0277777777777778}) - x20*(x147 + x82);
auto x14 = -x19*(V{0.0138888888888889}*x120*x85 + V{0.0277777777777778}) - x20*(x149 + x86);
auto x15 = x19*(V{0.0138888888888889}*x163*x90 + V{-0.0277777777777778}) - x20*(x114 + x153);
auto x16 = -x19*(V{0.0138888888888889}*x116*x120 + V{0.0277777777777778}) - x20*(x117 + x156);
auto x17 = x19*(V{0.0138888888888889}*x163*x93 + V{-0.0277777777777778}) - x20*(x160 + x30*x93);
auto x18 = -x19*(V{0.0138888888888889}*x120*x52 + V{0.0277777777777778}) - x20*(x122*x161 + x162 - x52*(-V{0.0196759259259259}*x24*x29 + V{0.0196759259259259}*x26*x28));
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
return { V{0.5}*x163, x37 + x39 + x41 };
}
};

}

}

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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::BGK> >, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x27 = parameters.template get<descriptors::OMEGA>();
auto x28 = cell[18] + cell[25] + cell[26];
auto x29 = cell[10] + cell[4];
auto x30 = cell[12] + x29;
auto x31 = cell[22] + cell[6];
auto x32 = cell[15] + cell[21];
auto x33 = cell[16] + cell[9];
auto x34 = cell[19] + cell[7];
auto x35 = cell[0] + cell[11] + cell[13] + cell[14] + cell[17] + cell[1] + cell[20] + cell[23] + cell[24] + cell[2] + cell[3] + cell[5] + cell[8] + x28 + x30 + x31 + x32 + x33 + x34 + V{1};
auto x36 = V{1} / (x35);
auto x37 = V{1}*x36;
auto x38 = -cell[26];
auto x39 = -cell[25];
auto x40 = -cell[18];
auto x41 = -cell[19];
auto x42 = cell[7] + x41;
auto x43 = -cell[23];
auto x44 = cell[11] + x43;
auto x45 = x42 + x44;
auto x46 = -cell[24];
auto x47 = -cell[17];
auto x48 = cell[5] + x47;
auto x49 = cell[13] + x46 + x48;
auto x50 = -cell[20];
auto x51 = -cell[14] + cell[1];
auto x52 = cell[6] + x50 + x51;
auto x53 = x30 + x38 + x39 + x40 + x45 + x49 + x52;
auto x54 = x53*x53;
auto x55 = V{0.666666666666667}*cell[10];
auto x56 = V{0.666666666666667}*cell[11];
auto x57 = V{0.666666666666667}*cell[12];
auto x58 = V{0.666666666666667}*cell[13];
auto x59 = V{0.666666666666667}*cell[23];
auto x60 = V{0.666666666666667}*cell[24];
auto x61 = V{0.666666666666667}*cell[25];
auto x62 = V{0.666666666666667}*cell[26];
auto x63 = -V{0.333333333333333}*cell[0];
auto x64 = -V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5] + x55 + x56 + x57 + x58 + x59 + x60 + x61 + x62 + x63;
auto x65 = -V{0.333333333333333}*cell[15] + V{0.666666666666667}*cell[19] + V{0.666666666666667}*cell[20] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x66 = V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[21] - V{0.333333333333333}*cell[22] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9] - x37*x54 + x64 + x65;
auto x67 = -cell[12];
auto x68 = -cell[13];
auto x69 = cell[8] + x68;
auto x70 = -cell[15] + cell[2];
auto x71 = x46 + x70;
auto x72 = -cell[5] + x47;
auto x73 = -cell[22];
auto x74 = -cell[21];
auto x75 = cell[9] + x74;
auto x76 = x73 + x75;
auto x77 = x28 + x29 + x44 + x67 + x69 + x71 + x72 + x76;
auto x78 = x77*x77;
auto x79 = -V{0.333333333333333}*cell[14] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[21] + V{0.666666666666667}*cell[22] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
auto x80 = V{0.666666666666667}*cell[15] - V{0.333333333333333}*cell[19] - V{0.333333333333333}*cell[20] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7] - x37*x78 + x64 + x79;
auto x81 = cell[26] + x39;
auto x82 = cell[20] + x81;
auto x83 = -cell[11];
auto x84 = cell[12] + x43 + x83;
auto x85 = -cell[16] + cell[3];
auto x86 = cell[10] + cell[24] + x85;
auto x87 = -cell[9];
auto x88 = x74 + x87;
auto x89 = -cell[7] + x41;
auto x90 = x31 + x69 + x82 + x84 + x86 + x88 + x89;
auto x91 = x90*x90;
auto x92 = V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5] - x37*x91 + x55 + x56 + x57 + x58 + x59 + x60 + x61 + x62 + x63 + x65 + x79;
auto x93 = x36*x90;
auto x94 = -cell[8];
auto x95 = cell[22] + x94;
auto x96 = cell[11] + cell[12] + x68;
auto x97 = cell[25] + x38 + x43;
auto x98 = -cell[10];
auto x99 = cell[24] + x98;
auto x100 = x75 + x77*x93 + x95 + x96 + x97 + x99;
auto x101 = -cell[6];
auto x102 = cell[13] + x67;
auto x103 = x101 + x102 + x45 + x53*x93 + x82 + x99;
auto x104 = -cell[4] + x28 + x36*x53*x77 + x49 + x84 + x98;
auto x105 = V{1} / (V{3.00000046417339}*util::sqrt(x36*(cell.template getFieldComponent<collision::LES::SMAGORINSKY>(0)*cell.template getFieldComponent<collision::LES::SMAGORINSKY>(0))*util::sqrt(x100*x100 + x103*x103 + x104*x104 + V{0.5}*(x66*x66) + V{0.5}*(x80*x80) + V{0.5}*(x92*x92)) + V{0.0277777691819762}/((x27)*(x27))) + V{0.5}/x27);
auto x106 = V{1} / ((x35)*(x35));
auto x107 = V{1.5}*x106;
auto x108 = -x53;
auto x109 = x108*x108;
auto x110 = x107*x109;
auto x111 = -x90;
auto x112 = x111*x111;
auto x113 = x107*x112;
auto x114 = -x77;
auto x115 = x114*x114;
auto x116 = x107*x115;
auto x117 = x113 + x116 + V{-1};
auto x118 = x110 + x117;
auto x119 = V{1} - x105;
auto x120 = V{0.0740740740740741}*cell[0] + V{0.0740740740740741}*cell[10] + V{0.0740740740740741}*cell[11] + V{0.0740740740740741}*cell[12] + V{0.0740740740740741}*cell[13] + V{0.0740740740740741}*cell[14] + V{0.0740740740740741}*cell[15] + V{0.0740740740740741}*cell[16] + V{0.0740740740740741}*cell[17] + V{0.0740740740740741}*cell[18] + V{0.0740740740740741}*cell[19] + V{0.0740740740740741}*cell[1] + V{0.0740740740740741}*cell[20] + V{0.0740740740740741}*cell[21] + V{0.0740740740740741}*cell[22] + V{0.0740740740740741}*cell[23] + V{0.0740740740740741}*cell[24] + V{0.0740740740740741}*cell[25] + V{0.0740740740740741}*cell[26] + V{0.0740740740740741}*cell[2] + V{0.0740740740740741}*cell[3] + V{0.0740740740740741}*cell[4] + V{0.0740740740740741}*cell[5] + V{0.0740740740740741}*cell[6] + V{0.0740740740740741}*cell[7] + V{0.0740740740740741}*cell[8] + V{0.0740740740740741}*cell[9] + V{0.0740740740740741};
auto x121 = V{3}*cell[13];
auto x122 = V{3}*cell[5];
auto x123 = V{3}*cell[7];
auto x124 = V{3}*cell[18];
auto x125 = V{3}*cell[20];
auto x126 = V{3}*cell[26];
auto x127 = V{3}*cell[10];
auto x128 = V{3}*cell[11];
auto x129 = -V{3}*cell[23];
auto x130 = V{3}*cell[24];
auto x131 = -V{3}*cell[17] + V{3}*cell[4] + x127 + x128 + x129 - x130;
auto x132 = V{3}*cell[12];
auto x133 = V{3}*cell[25];
auto x134 = -V{3}*cell[19] + V{3}*cell[6] + x132 - x133;
auto x135 = -V{3}*cell[14] + V{3}*cell[1] + x121 + x122 + x123 - x124 - x125 - x126 + x131 + x134;
auto x136 = -x135*x36;
auto x137 = V{3}*x106;
auto x138 = V{3}*cell[9];
auto x139 = V{3}*cell[22];
auto x140 = -V{3}*cell[21] + V{3}*cell[8] - x121 + x126;
auto x141 = -V{3}*cell[15] + V{3}*cell[2] - x122 + x124 + x131 - x132 + x133 + x138 - x139 + x140;
auto x142 = -x141*x36;
auto x143 = x110 + V{-1};
auto x144 = -V{3}*cell[16] + V{3}*cell[3] - x123 + x125 + x127 - x128 + x129 + x130 + x134 - x138 + x139 + x140;
auto x145 = -x144*x36;
auto x146 = V{0.0185185185185185}*cell[0] + V{0.0185185185185185}*cell[10] + V{0.0185185185185185}*cell[11] + V{0.0185185185185185}*cell[12] + V{0.0185185185185185}*cell[13] + V{0.0185185185185185}*cell[14] + V{0.0185185185185185}*cell[15] + V{0.0185185185185185}*cell[16] + V{0.0185185185185185}*cell[17] + V{0.0185185185185185}*cell[18] + V{0.0185185185185185}*cell[19] + V{0.0185185185185185}*cell[1] + V{0.0185185185185185}*cell[20] + V{0.0185185185185185}*cell[21] + V{0.0185185185185185}*cell[22] + V{0.0185185185185185}*cell[23] + V{0.0185185185185185}*cell[24] + V{0.0185185185185185}*cell[25] + V{0.0185185185185185}*cell[26] + V{0.0185185185185185}*cell[2] + V{0.0185185185185185}*cell[3] + V{0.0185185185185185}*cell[4] + V{0.0185185185185185}*cell[5] + V{0.0185185185185185}*cell[6] + V{0.0185185185185185}*cell[7] + V{0.0185185185185185}*cell[8] + V{0.0185185185185185}*cell[9] + V{0.0185185185185185};
auto x147 = V{4.5}*x106;
auto x148 = x42 + x52;
auto x149 = -V{2}*cell[17] + V{2}*cell[4];
auto x150 = -V{2}*cell[23];
auto x151 = V{2}*cell[10];
auto x152 = cell[8] + x150 + x151;
auto x153 = V{2}*cell[11] - V{2}*cell[24];
auto x154 = x153 + x70;
auto x155 = x148 + x149 + x152 + x154 + x76;
auto x156 = -x155;
auto x157 = x118 + x136;
auto x158 = x142 + x157;
auto x159 = V{2}*cell[25];
auto x160 = V{2}*cell[12];
auto x161 = -x159 + x160;
auto x162 = V{2}*cell[26];
auto x163 = V{2}*cell[13];
auto x164 = -x162 + x163;
auto x165 = V{2}*cell[18];
auto x166 = V{2}*cell[5];
auto x167 = -cell[2] - x165 + x166;
auto x168 = x148 + x161 + x164 + x167 + x32 + x87 + x95;
auto x169 = -x142;
auto x170 = x157 + x169;
auto x171 = cell[4] + x85;
auto x172 = -V{2}*cell[19] + V{2}*cell[6] + x51;
auto x173 = x40 + x48;
auto x174 = cell[22] + x152 + x161 + x171 + x172 + x173 + x88;
auto x175 = -x174;
auto x176 = -x145;
auto x177 = -cell[3];
auto x178 = cell[4] + x177;
auto x179 = V{2}*cell[20];
auto x180 = V{2}*cell[7];
auto x181 = -x179 + x180 + x51;
auto x182 = cell[21] + x153 + x164 + x173 + x178 + x181 + x33 + x73 + x94;
auto x183 = cell[18] + x72;
auto x184 = -V{2}*cell[21] + V{2}*cell[8];
auto x185 = x184 + x70;
auto x186 = cell[20] + cell[6] + x150 + x151 + x162 - x163 + x171 + x183 + x185 + x89;
auto x187 = -x186;
auto x188 = x118 + x142;
auto x189 = x145 + x188;
auto x190 = V{2}*cell[22];
auto x191 = V{2}*cell[9];
auto x192 = cell[16] - x190 + x191;
auto x193 = x101 + x154 + x159 - x160 + x178 + x183 + x192 + x34 + x50;
auto x194 = V{0.00462962962962963}*cell[0] + V{0.00462962962962963}*cell[10] + V{0.00462962962962963}*cell[11] + V{0.00462962962962963}*cell[12] + V{0.00462962962962963}*cell[13] + V{0.00462962962962963}*cell[14] + V{0.00462962962962963}*cell[15] + V{0.00462962962962963}*cell[16] + V{0.00462962962962963}*cell[17] + V{0.00462962962962963}*cell[18] + V{0.00462962962962963}*cell[19] + V{0.00462962962962963}*cell[1] + V{0.00462962962962963}*cell[20] + V{0.00462962962962963}*cell[21] + V{0.00462962962962963}*cell[22] + V{0.00462962962962963}*cell[23] + V{0.00462962962962963}*cell[24] + V{0.00462962962962963}*cell[25] + V{0.00462962962962963}*cell[26] + V{0.00462962962962963}*cell[2] + V{0.00462962962962963}*cell[3] + V{0.00462962962962963}*cell[4] + V{0.00462962962962963}*cell[5] + V{0.00462962962962963}*cell[6] + V{0.00462962962962963}*cell[7] + V{0.00462962962962963}*cell[8] + V{0.00462962962962963}*cell[9] + V{0.00462962962962963};
auto x195 = V{3}*cell[10] - V{3}*cell[23] + x149 + x172 + x184 + x71 + x81 + x85 + x96;
auto x196 = -x195;
auto x197 = cell[10] + V{3}*cell[11] - V{3}*cell[24] + x102 + x149 + x177 + x181 + x192 + x70 + x97;
auto x198 = x43 + x83 + x86;
auto x199 = V{3}*cell[12] + cell[13] + cell[15] - V{3}*cell[25] + x167 + x172 + x190 - x191 + x198 + x38;
auto x200 = x135*x36;
auto x201 = -V{3}*cell[13] + cell[14] - cell[1] + cell[25] + V{3}*cell[26] + x165 - x166 + x179 - x180 + x185 + x198 + x67;
auto x202 = -x201;
auto x203 = x144*x36;
auto x204 = x107*x78;
auto x205 = x107*x91 + V{-1};
auto x206 = x204 + x205;
auto x207 = x141*x36;
auto x208 = x107*x54;
auto x209 = x207 + x208;
auto x210 = x203 + x206 + x209;
auto x211 = x200 + x206;
auto x212 = x203 + x208;
auto x213 = x209 + x211;
auto x214 = -x136;
auto x215 = -x168;
auto x216 = x211 + x212;
auto x217 = -x182;
auto x218 = x118 + x145;
auto x219 = -x193;
auto x220 = -x197;
auto x221 = -x199;
auto x0 = V{1}*cell[0]*x119 - x105*(x118*(V{0.296296296296296}*cell[0] + V{0.296296296296296}*cell[10] + V{0.296296296296296}*cell[11] + V{0.296296296296296}*cell[12] + V{0.296296296296296}*cell[13] + V{0.296296296296296}*cell[14] + V{0.296296296296296}*cell[15] + V{0.296296296296296}*cell[16] + V{0.296296296296296}*cell[17] + V{0.296296296296296}*cell[18] + V{0.296296296296296}*cell[19] + V{0.296296296296296}*cell[1] + V{0.296296296296296}*cell[20] + V{0.296296296296296}*cell[21] + V{0.296296296296296}*cell[22] + V{0.296296296296296}*cell[23] + V{0.296296296296296}*cell[24] + V{0.296296296296296}*cell[25] + V{0.296296296296296}*cell[26] + V{0.296296296296296}*cell[2] + V{0.296296296296296}*cell[3] + V{0.296296296296296}*cell[4] + V{0.296296296296296}*cell[5] + V{0.296296296296296}*cell[6] + V{0.296296296296296}*cell[7] + V{0.296296296296296}*cell[8] + V{0.296296296296296}*cell[9] + V{0.296296296296296}) + V{0.296296296296296});
auto x1 = V{1}*cell[1]*x119 - x105*(x120*(-x109*x137 + x117 + x136) + V{0.0740740740740741});
auto x2 = V{1}*cell[2]*x119 - x105*(x120*(x113 - x115*x137 + x142 + x143) + V{0.0740740740740741});
auto x3 = V{1}*cell[3]*x119 - x105*(x120*(-x112*x137 + x116 + x143 + x145) + V{0.0740740740740741});
auto x4 = V{1}*cell[4]*x119 - x105*(x146*(-x147*x156*x156 + x158) + V{0.0185185185185185});
auto x5 = V{1}*cell[5]*x119 - x105*(x146*(-x147*x168*x168 + x170) + V{0.0185185185185185});
auto x6 = V{1}*cell[6]*x119 - x105*(x146*(x145 - x147*x175*x175 + x157) + V{0.0185185185185185});
auto x7 = V{1}*cell[7]*x119 - x105*(x146*(-x147*x182*x182 + x157 + x176) + V{0.0185185185185185});
auto x8 = V{1}*cell[8]*x119 - x105*(x146*(-x147*x187*x187 + x189) + V{0.0185185185185185});
auto x9 = V{1}*cell[9]*x119 - x105*(x146*(-x147*x193*x193 + x176 + x188) + V{0.0185185185185185});
auto x10 = V{1}*cell[10]*x119 - x105*(x194*(x145 - x147*x196*x196 + x158) + V{0.00462962962962963});
auto x11 = V{1}*cell[11]*x119 - x105*(x194*(-x147*x197*x197 + x158 + x176) + V{0.00462962962962963});
auto x12 = V{1}*cell[12]*x119 - x105*(x194*(x145 - x147*x199*x199 + x170) + V{0.00462962962962963});
auto x13 = V{1}*cell[13]*x119 - x105*(x194*(-x147*x202*x202 - x200 + x210) + V{0.00462962962962963});
auto x14 = V{1}*cell[14]*x119 - x105*(x120*(-x137*x54 + x211) + V{0.0740740740740741});
auto x15 = V{1}*cell[15]*x119 - x105*(x120*(-x137*x78 + x205 + x209) + V{0.0740740740740741});
auto x16 = V{1}*cell[16]*x119 - x105*(x120*(-x137*x91 + x204 + x212 + V{-1}) + V{0.0740740740740741});
auto x17 = V{1}*cell[17]*x119 - x105*(x146*(-x147*x155*x155 + x213) + V{0.0185185185185185});
auto x18 = V{1}*cell[18]*x119 - x105*(x146*(-x147*x215*x215 + x188 + x214) + V{0.0185185185185185});
auto x19 = V{1}*cell[19]*x119 - x105*(x146*(-x147*x174*x174 + x216) + V{0.0185185185185185});
auto x20 = V{1}*cell[20]*x119 - x105*(x146*(-x147*x217*x217 + x214 + x218) + V{0.0185185185185185});
auto x21 = V{1}*cell[21]*x119 - x105*(x146*(-x147*x186*x186 + x210) + V{0.0185185185185185});
auto x22 = V{1}*cell[22]*x119 - x105*(x146*(-x147*x219*x219 + x169 + x218) + V{0.0185185185185185});
auto x23 = V{1}*cell[23]*x119 - x105*(x194*(-x147*x195*x195 + x203 + x213) + V{0.00462962962962963});
auto x24 = V{1}*cell[24]*x119 - x105*(x194*(-x147*x220*x220 - x203 + x213) + V{0.00462962962962963});
auto x25 = V{1}*cell[25]*x119 - x105*(x194*(-x147*x221*x221 - x207 + x216) + V{0.00462962962962963});
auto x26 = V{1}*cell[26]*x119 - x105*(x194*(-x147*x201*x201 + x189 + x214) + V{0.00462962962962963});
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
cell[19] = x19;
cell[20] = x20;
cell[21] = x21;
cell[22] = x22;
cell[23] = x23;
cell[24] = x24;
cell[25] = x25;
cell[26] = x26;
return { x35, V{1}*x106*(x54 + x78 + x91) };
}
};

}

}

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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::ForcedWithStress> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x27 = parameters.template get<descriptors::OMEGA>();
auto x28 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x29 = cell[13] + cell[1] + cell[5] + cell[7];
auto x30 = cell[18] + cell[25] + cell[2] + cell[9];
auto x31 = cell[10] + cell[20] + cell[22] + cell[24] + cell[3];
auto x32 = cell[0] + cell[11] + cell[12] + cell[14] + cell[15] + cell[16] + cell[17] + cell[19] + cell[21] + cell[23] + cell[26] + cell[4] + cell[6] + cell[8] + x29 + x30 + x31;
auto x33 = x32 + V{1};
auto x34 = V{1} / (x33);
auto x35 = -cell[23];
auto x36 = cell[10] + cell[11] - cell[17] - cell[24] + cell[4] + x35;
auto x37 = cell[12] - cell[19] - cell[25] + cell[6];
auto x38 = -cell[14] - cell[18] - cell[20] - cell[26] + x29 + x36 + x37;
auto x39 = x34*(x32 + V{1});
auto x40 = V{0.666666666666667}*cell[10];
auto x41 = V{0.666666666666667}*cell[11];
auto x42 = V{0.666666666666667}*cell[12];
auto x43 = V{0.666666666666667}*cell[13];
auto x44 = V{0.666666666666667}*cell[23];
auto x45 = V{0.666666666666667}*cell[24];
auto x46 = V{0.666666666666667}*cell[25];
auto x47 = V{0.666666666666667}*cell[26];
auto x48 = -V{0.333333333333333}*cell[0];
auto x49 = -V{0.333333333333333}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] - V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5] + x40 + x41 + x42 + x43 + x44 + x45 + x46 + x47 + x48;
auto x50 = -V{0.333333333333333}*cell[15] + V{0.666666666666667}*cell[19] + V{0.666666666666667}*cell[20] - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7];
auto x51 = -cell.template getFieldComponent<descriptors::FORCE>(0)*x38*x39 + V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[1] - V{0.333333333333333}*cell[21] - V{0.333333333333333}*cell[22] - V{0.333333333333333}*cell[8] - V{0.333333333333333}*cell[9] - x34*x38*x38 + x49 + x50;
auto x52 = -cell[13] - cell[21] + cell[26] + cell[8];
auto x53 = -cell[12] - cell[15] - cell[22] - cell[5] + x30 + x36 + x52;
auto x54 = -V{0.333333333333333}*cell[14] - V{0.333333333333333}*cell[1] + V{0.666666666666667}*cell[21] + V{0.666666666666667}*cell[22] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9];
auto x55 = -cell.template getFieldComponent<descriptors::FORCE>(1)*x39*x53 + V{0.666666666666667}*cell[15] - V{0.333333333333333}*cell[19] - V{0.333333333333333}*cell[20] + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[7] - x34*x53*x53 + x49 + x54;
auto x56 = -cell[11] - cell[16] - cell[7] - cell[9] + x31 + x35 + x37 + x52;
auto x57 = -cell.template getFieldComponent<descriptors::FORCE>(2)*x39*x56 + V{0.666666666666667}*cell[16] - V{0.333333333333333}*cell[17] - V{0.333333333333333}*cell[18] + V{0.666666666666667}*cell[3] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[5] - x34*x56*x56 + x40 + x41 + x42 + x43 + x44 + x45 + x46 + x47 + x48 + x50 + x54;
auto x58 = V{1}*cell[18];
auto x59 = V{1}*cell[11];
auto x60 = -x59;
auto x61 = V{1}*cell[4];
auto x62 = x39*(cell.template getFieldComponent<descriptors::FORCE>(0)*x53 + cell.template getFieldComponent<descriptors::FORCE>(1)*x38);
auto x63 = x34*x53;
auto x64 = x38*x63;
auto x65 = V{1}*cell[12];
auto x66 = V{1}*cell[25];
auto x67 = V{1}*cell[10];
auto x68 = -x67;
auto x69 = -V{1}*cell[23];
auto x70 = x65 + x66 + x68 + x69;
auto x71 = V{1}*cell[13];
auto x72 = V{1}*cell[26];
auto x73 = x71 + x72;
auto x74 = V{1}*cell[5];
auto x75 = -V{1}*cell[17];
auto x76 = V{1}*cell[24];
auto x77 = -x76;
auto x78 = x74 + x75 + x77;
auto x79 = V{2}*cell[13];
auto x80 = V{2}*cell[26];
auto x81 = V{2}*cell[11];
auto x82 = V{2}*cell[24];
auto x83 = -V{2}*cell[10] + V{2}*cell[12] - V{2}*cell[23] + V{2}*cell[25];
auto x84 = V{1}*cell[22];
auto x85 = -x72;
auto x86 = V{1}*cell[8];
auto x87 = x39*(cell.template getFieldComponent<descriptors::FORCE>(1)*x56 + cell.template getFieldComponent<descriptors::FORCE>(2)*x53);
auto x88 = x56*x63;
auto x89 = x59 + x76;
auto x90 = V{1}*cell[9];
auto x91 = -x71;
auto x92 = -V{1}*cell[21];
auto x93 = x90 + x91 + x92;
auto x94 = V{1}*cell[20];
auto x95 = V{1}*cell[6];
auto x96 = -x65 + x69;
auto x97 = V{1}*cell[7];
auto x98 = -V{1}*cell[19];
auto x99 = -x66;
auto x100 = x97 + x98 + x99;
auto x101 = -x100 - V{1}*x34*x38*x56 - V{0.5}*x39*(cell.template getFieldComponent<descriptors::FORCE>(0)*x56 + cell.template getFieldComponent<descriptors::FORCE>(2)*x38) - x68 - x73 - x89 - x94 + x95 - x96;
auto x102 = V{1} / (V{3.00000046417339}*util::sqrt(x34*(x28*x28)*util::sqrt(V{0.5}*(-x58 - x60 + x61 - V{0.5}*x62 - V{1}*x64 - x70 - x73 - x78)*(V{2}*cell[17] - V{2}*cell[18] + V{2}*cell[4] - V{2}*cell[5] - V{1}*x62 - V{2}*x64 - x79 - x80 + x81 + x82 - x83) + V{0.5}*(-x70 - x84 - x85 + x86 - V{0.5}*x87 - V{1}*x88 - x89 - x93)*(V{2}*cell[21] - V{2}*cell[22] + V{2}*cell[8] - V{2}*cell[9] + x79 + x80 - x81 - x82 - x83 - V{1}*x87 - V{2}*x88) + x101*x101 + V{0.5}*(x51*x51) + V{0.5}*(x55*x55) + V{0.5}*(x57*x57)) + V{0.0277777691819762}/((x27)*(x27))) + V{0.5}/x27);
auto x103 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x104 = x59 + x61 + x67;
auto x105 = x65 + x69 + x95;
auto x106 = -V{1}*cell[14] + V{1}*cell[1] + x100 + x104 + x105 - x58 + x71 + x78 + x85 - x94;
auto x107 = -x106*x34;
auto x108 = x103 + x107;
auto x109 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x110 = x72 + x86;
auto x111 = -V{1}*cell[15] + V{1}*cell[2] + x104 + x110 + x58 + x66 - x74 + x75 + x77 - x84 + x93 + x96;
auto x112 = -x111*x34;
auto x113 = x109 + x112;
auto x114 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x115 = -V{1}*cell[16] + V{1}*cell[3] + x105 + x110 + x60 + x67 + x76 + x84 - x90 + x91 + x92 + x94 - x97 + x98 + x99;
auto x116 = -x115*x34;
auto x117 = x114 + x116;
auto x118 = V{-1} + V{1.5}*(x108*x108) + V{1.5}*(x113*x113) + V{1.5}*(x117*x117);
auto x119 = V{1} - x102;
auto x120 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
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
auto x136 = -x135*x34;
auto x137 = x120 + x136;
auto x138 = cell.template getFieldComponent<descriptors::FORCE>(0)*x137;
auto x139 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x140 = V{3}*cell[9];
auto x141 = V{3}*cell[22];
auto x142 = -V{3}*cell[21] + V{3}*cell[8] - x121 + x126;
auto x143 = -V{3}*cell[15] + V{3}*cell[2] - x122 + x124 + x131 - x132 + x133 + x140 - x141 + x142;
auto x144 = -x143*x34;
auto x145 = x139 + x144;
auto x146 = cell.template getFieldComponent<descriptors::FORCE>(1)*x145;
auto x147 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x148 = -V{3}*cell[16] + V{3}*cell[3] - x123 + x125 + x127 - x128 + x129 + x130 + x134 - x140 + x141 + x142;
auto x149 = -x148*x34;
auto x150 = x147 + x149;
auto x151 = cell.template getFieldComponent<descriptors::FORCE>(2)*x150;
auto x152 = x146 + x151;
auto x153 = V{1} - V{0.5}*x102;
auto x154 = V{0.0740740740740741}*cell[0] + V{0.0740740740740741}*cell[10] + V{0.0740740740740741}*cell[11] + V{0.0740740740740741}*cell[12] + V{0.0740740740740741}*cell[13] + V{0.0740740740740741}*cell[14] + V{0.0740740740740741}*cell[15] + V{0.0740740740740741}*cell[16] + V{0.0740740740740741}*cell[17] + V{0.0740740740740741}*cell[18] + V{0.0740740740740741}*cell[19] + V{0.0740740740740741}*cell[1] + V{0.0740740740740741}*cell[20] + V{0.0740740740740741}*cell[21] + V{0.0740740740740741}*cell[22] + V{0.0740740740740741}*cell[23] + V{0.0740740740740741}*cell[24] + V{0.0740740740740741}*cell[25] + V{0.0740740740740741}*cell[26] + V{0.0740740740740741}*cell[2] + V{0.0740740740740741}*cell[3] + V{0.0740740740740741}*cell[4] + V{0.0740740740740741}*cell[5] + V{0.0740740740740741}*cell[6] + V{0.0740740740740741}*cell[7] + V{0.0740740740740741}*cell[8] + V{0.0740740740740741}*cell[9] + V{0.0740740740740741};
auto x155 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x156 = V{4.5}*cell[13];
auto x157 = V{4.5}*cell[5];
auto x158 = V{4.5}*cell[7];
auto x159 = V{4.5}*cell[18];
auto x160 = V{4.5}*cell[20];
auto x161 = V{4.5}*cell[26];
auto x162 = V{4.5}*cell[10];
auto x163 = V{4.5}*cell[11];
auto x164 = -V{4.5}*cell[23];
auto x165 = V{4.5}*cell[24];
auto x166 = -V{4.5}*cell[17] + V{4.5}*cell[4] + x162 + x163 + x164 - x165;
auto x167 = V{4.5}*cell[12];
auto x168 = V{4.5}*cell[25];
auto x169 = -V{4.5}*cell[19] + V{4.5}*cell[6] + x167 - x168;
auto x170 = -V{4.5}*cell[14] + V{4.5}*cell[1] + x156 + x157 + x158 - x159 - x160 - x161 + x166 + x169;
auto x171 = -x170*x34;
auto x172 = x155 + x171;
auto x173 = x118 + x137;
auto x174 = V{6}*cell[13];
auto x175 = V{6}*cell[5];
auto x176 = V{6}*cell[7];
auto x177 = V{6}*cell[18];
auto x178 = V{6}*cell[20];
auto x179 = V{6}*cell[26];
auto x180 = V{6}*cell[10];
auto x181 = V{6}*cell[11];
auto x182 = -V{6}*cell[23];
auto x183 = V{6}*cell[24];
auto x184 = -V{6}*cell[17] + V{6}*cell[4] + x180 + x181 + x182 - x183;
auto x185 = V{6}*cell[12];
auto x186 = V{6}*cell[25];
auto x187 = -V{6}*cell[19] + V{6}*cell[6] + x185 - x186;
auto x188 = -V{6}*cell[14] + V{6}*cell[1] + x174 + x175 + x176 - x177 - x178 - x179 + x184 + x187;
auto x189 = x188*x34;
auto x190 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x191 = V{3} - x190;
auto x192 = x189 + x191;
auto x193 = V{0.0740740740740741}*x33;
auto x194 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x195 = V{4.5}*cell[9];
auto x196 = V{4.5}*cell[22];
auto x197 = -V{4.5}*cell[21] + V{4.5}*cell[8] - x156 + x161;
auto x198 = -V{4.5}*cell[15] + V{4.5}*cell[2] - x157 + x159 + x166 - x167 + x168 + x195 - x196 + x197;
auto x199 = -x198*x34;
auto x200 = x194 + x199;
auto x201 = x118 + x145;
auto x202 = V{6}*cell[9];
auto x203 = V{6}*cell[22];
auto x204 = -V{6}*cell[21] + V{6}*cell[8] - x174 + x179;
auto x205 = -V{6}*cell[15] + V{6}*cell[2] - x175 + x177 + x184 - x185 + x186 + x202 - x203 + x204;
auto x206 = x205*x34;
auto x207 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x208 = V{3} - x207;
auto x209 = x206 + x208;
auto x210 = x138 + x151;
auto x211 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x212 = -V{4.5}*cell[16] + V{4.5}*cell[3] - x158 + x160 + x162 - x163 + x164 + x165 + x169 - x195 + x196 + x197;
auto x213 = -x212*x34;
auto x214 = x211 + x213;
auto x215 = x118 + x150;
auto x216 = -V{6}*cell[16] + V{6}*cell[3] - x176 + x178 + x180 - x181 + x182 + x183 + x187 - x202 + x203 + x204;
auto x217 = x216*x34;
auto x218 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x219 = V{3} - x218;
auto x220 = x217 + x219;
auto x221 = x138 + x146;
auto x222 = V{0.0185185185185185}*cell[0] + V{0.0185185185185185}*cell[10] + V{0.0185185185185185}*cell[11] + V{0.0185185185185185}*cell[12] + V{0.0185185185185185}*cell[13] + V{0.0185185185185185}*cell[14] + V{0.0185185185185185}*cell[15] + V{0.0185185185185185}*cell[16] + V{0.0185185185185185}*cell[17] + V{0.0185185185185185}*cell[18] + V{0.0185185185185185}*cell[19] + V{0.0185185185185185}*cell[1] + V{0.0185185185185185}*cell[20] + V{0.0185185185185185}*cell[21] + V{0.0185185185185185}*cell[22] + V{0.0185185185185185}*cell[23] + V{0.0185185185185185}*cell[24] + V{0.0185185185185185}*cell[25] + V{0.0185185185185185}*cell[26] + V{0.0185185185185185}*cell[2] + V{0.0185185185185185}*cell[3] + V{0.0185185185185185}*cell[4] + V{0.0185185185185185}*cell[5] + V{0.0185185185185185}*cell[6] + V{0.0185185185185185}*cell[7] + V{0.0185185185185185}*cell[8] + V{0.0185185185185185}*cell[9] + V{0.0185185185185185};
auto x223 = x108 + x113;
auto x224 = x172 + x200;
auto x225 = x145 + x173;
auto x226 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x227 = V{9}*cell[18];
auto x228 = V{9}*cell[25];
auto x229 = V{9}*cell[9];
auto x230 = V{9}*cell[12];
auto x231 = V{9}*cell[22];
auto x232 = V{9}*cell[5];
auto x233 = V{9}*cell[10];
auto x234 = V{9}*cell[11];
auto x235 = -V{9}*cell[23];
auto x236 = V{9}*cell[24];
auto x237 = -V{9}*cell[17] + V{9}*cell[4] + x233 + x234 + x235 - x236;
auto x238 = V{9}*cell[26];
auto x239 = V{9}*cell[13];
auto x240 = -V{9}*cell[21] + V{9}*cell[8] + x238 - x239;
auto x241 = -V{9}*cell[15] + V{9}*cell[2] + x227 + x228 + x229 - x230 - x231 - x232 + x237 + x240;
auto x242 = x241*x34;
auto x243 = -x226 + x242;
auto x244 = x192 + x243;
auto x245 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x246 = V{9}*cell[7];
auto x247 = V{9}*cell[20];
auto x248 = -V{9}*cell[19] + V{9}*cell[6] - x228 + x230;
auto x249 = -V{9}*cell[14] + V{9}*cell[1] - x227 + x232 + x237 - x238 + x239 + x246 - x247 + x248;
auto x250 = x249*x34;
auto x251 = -x245 + x250;
auto x252 = x209 + x251;
auto x253 = x148*x34;
auto x254 = x147 - x253;
auto x255 = cell.template getFieldComponent<descriptors::FORCE>(2)*x254;
auto x256 = V{0.0185185185185185}*x33;
auto x257 = -x109;
auto x258 = x108 - x112 + x257;
auto x259 = -x194;
auto x260 = x172 - x199 + x259;
auto x261 = -x139;
auto x262 = -x144 + x261;
auto x263 = x173 + x262;
auto x264 = x188*x34 + x191;
auto x265 = x226 - x241*x34;
auto x266 = -x206 + x207 + V{3};
auto x267 = x251 + x266;
auto x268 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x269 = -V{9}*cell[16] + V{9}*cell[3] - x229 + x231 + x233 - x234 + x235 + x236 + x240 - x246 + x247 + x248;
auto x270 = x269*x34;
auto x271 = -x268 + x270;
auto x272 = x192 + x271;
auto x273 = x143*x34;
auto x274 = x139 - x273;
auto x275 = cell.template getFieldComponent<descriptors::FORCE>(1)*x274;
auto x276 = x220 + x251;
auto x277 = -x114;
auto x278 = -x116 + x277;
auto x279 = x108 + x278;
auto x280 = -x211;
auto x281 = -x213 + x280;
auto x282 = x172 + x281;
auto x283 = -x147;
auto x284 = -x149 + x283;
auto x285 = x268 - x269*x34;
auto x286 = -x217 + x218 + V{3};
auto x287 = x251 + x286;
auto x288 = x113 + x117;
auto x289 = x200 + x214;
auto x290 = x150 + x201;
auto x291 = x135*x34;
auto x292 = x120 - x291;
auto x293 = cell.template getFieldComponent<descriptors::FORCE>(0)*x292;
auto x294 = x209 + x271;
auto x295 = x220 + x243;
auto x296 = x113 + x278;
auto x297 = x200 + x281;
auto x298 = x205*x34 + x208;
auto x299 = x243 + x286;
auto x300 = V{0.00462962962962963}*cell[0] + V{0.00462962962962963}*cell[10] + V{0.00462962962962963}*cell[11] + V{0.00462962962962963}*cell[12] + V{0.00462962962962963}*cell[13] + V{0.00462962962962963}*cell[14] + V{0.00462962962962963}*cell[15] + V{0.00462962962962963}*cell[16] + V{0.00462962962962963}*cell[17] + V{0.00462962962962963}*cell[18] + V{0.00462962962962963}*cell[19] + V{0.00462962962962963}*cell[1] + V{0.00462962962962963}*cell[20] + V{0.00462962962962963}*cell[21] + V{0.00462962962962963}*cell[22] + V{0.00462962962962963}*cell[23] + V{0.00462962962962963}*cell[24] + V{0.00462962962962963}*cell[25] + V{0.00462962962962963}*cell[26] + V{0.00462962962962963}*cell[2] + V{0.00462962962962963}*cell[3] + V{0.00462962962962963}*cell[4] + V{0.00462962962962963}*cell[5] + V{0.00462962962962963}*cell[6] + V{0.00462962962962963}*cell[7] + V{0.00462962962962963}*cell[8] + V{0.00462962962962963}*cell[9] + V{0.00462962962962963};
auto x301 = V{0.00462962962962963}*x33;
auto x302 = x268 - x270;
auto x303 = x226 - x242;
auto x304 = x103 - x106*x34;
auto x305 = x115*x34;
auto x306 = x277 + x305;
auto x307 = x111*x34;
auto x308 = x257 + x307;
auto x309 = x155 - x170*x34;
auto x310 = x212*x34;
auto x311 = x280 + x310;
auto x312 = x198*x34;
auto x313 = x259 + x312;
auto x314 = x304*x304;
auto x315 = x109 - x307;
auto x316 = x315*x315;
auto x317 = x114 - x305;
auto x318 = x317*x317;
auto x319 = V{1.5}*x314 + V{1.5}*x316 + V{1.5}*x318 + V{-1};
auto x320 = x261 + x273;
auto x321 = x319 + x320;
auto x322 = x253 + x283;
auto x323 = x321 + x322;
auto x324 = -x120;
auto x325 = x291 + x319 + x324;
auto x326 = -x189 + x190 + V{3};
auto x327 = x194 - x312;
auto x328 = x211 - x310;
auto x329 = x304 + x315;
auto x330 = x309 + x327;
auto x331 = x320 + x325;
auto x332 = x303 + x326;
auto x333 = x245 - x250;
auto x334 = x266 + x333;
auto x335 = -x136 + x324;
auto x336 = x245 - x249*x34;
auto x337 = x243 + x326;
auto x338 = x304 + x317;
auto x339 = x309 + x328;
auto x340 = x322 + x325;
auto x341 = x286 + x333;
auto x342 = x216*x34 + x219;
auto x0 = V{1}*cell[0]*x119 - x102*(x118*(V{0.296296296296296}*cell[0] + V{0.296296296296296}*cell[10] + V{0.296296296296296}*cell[11] + V{0.296296296296296}*cell[12] + V{0.296296296296296}*cell[13] + V{0.296296296296296}*cell[14] + V{0.296296296296296}*cell[15] + V{0.296296296296296}*cell[16] + V{0.296296296296296}*cell[17] + V{0.296296296296296}*cell[18] + V{0.296296296296296}*cell[19] + V{0.296296296296296}*cell[1] + V{0.296296296296296}*cell[20] + V{0.296296296296296}*cell[21] + V{0.296296296296296}*cell[22] + V{0.296296296296296}*cell[23] + V{0.296296296296296}*cell[24] + V{0.296296296296296}*cell[25] + V{0.296296296296296}*cell[26] + V{0.296296296296296}*cell[2] + V{0.296296296296296}*cell[3] + V{0.296296296296296}*cell[4] + V{0.296296296296296}*cell[5] + V{0.296296296296296}*cell[6] + V{0.296296296296296}*cell[7] + V{0.296296296296296}*cell[8] + V{0.296296296296296}*cell[9] + V{0.296296296296296}) + V{0.296296296296296}) - V{0.296296296296296}*x153*x33*(x138 + x152);
auto x1 = V{1}*cell[1]*x119 - x102*(x154*(-x108*x172 + x173) + V{0.0740740740740741}) - x153*x193*(cell.template getFieldComponent<descriptors::FORCE>(0)*x192 + x152);
auto x2 = V{1}*cell[2]*x119 - x102*(x154*(-x113*x200 + x201) + V{0.0740740740740741}) - x153*x193*(cell.template getFieldComponent<descriptors::FORCE>(1)*x209 + x210);
auto x3 = V{1}*cell[3]*x119 - x102*(x154*(-x117*x214 + x215) + V{0.0740740740740741}) - x153*x193*(cell.template getFieldComponent<descriptors::FORCE>(2)*x220 + x221);
auto x4 = V{1}*cell[4]*x119 - x102*(x222*(-x223*x224 + x225) + V{0.0185185185185185}) - x153*x256*(cell.template getFieldComponent<descriptors::FORCE>(0)*x244 + cell.template getFieldComponent<descriptors::FORCE>(1)*x252 + x255);
auto x5 = V{1}*cell[5]*x119 - x102*(x222*(-x258*x260 + x263) + V{0.0185185185185185}) - x153*x256*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x264 + x265) - cell.template getFieldComponent<descriptors::FORCE>(1)*x267 + x151);
auto x6 = V{1}*cell[6]*x119 - x102*(x222*(x150 + x173 - (x108 + x117)*(x172 + x214)) + V{0.0185185185185185}) - x153*x256*(cell.template getFieldComponent<descriptors::FORCE>(0)*x272 + cell.template getFieldComponent<descriptors::FORCE>(2)*x276 + x275);
auto x7 = V{1}*cell[7]*x119 - x102*(x222*(x173 - x279*x282 + x284) + V{0.0185185185185185}) - x153*x256*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x264 + x285) - cell.template getFieldComponent<descriptors::FORCE>(2)*x287 + x146);
auto x8 = V{1}*cell[8]*x119 - x102*(x222*(-x288*x289 + x290) + V{0.0185185185185185}) - x153*x256*(cell.template getFieldComponent<descriptors::FORCE>(1)*x294 + cell.template getFieldComponent<descriptors::FORCE>(2)*x295 + x293);
auto x9 = V{1}*cell[9]*x119 - x102*(x222*(x201 + x284 - x296*x297) + V{0.0185185185185185}) - x153*x256*(cell.template getFieldComponent<descriptors::FORCE>(1)*(x285 + x298) - cell.template getFieldComponent<descriptors::FORCE>(2)*x299 + x138);
auto x10 = V{1}*cell[10]*x119 - x102*(x300*(x150 + x225 - (x117 + x223)*(x214 + x224)) + V{0.00462962962962963}) - x153*x301*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x244 + x271) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x252 + x271) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x243 + x276));
auto x11 = V{1}*cell[11]*x119 - x102*(x300*(x225 + x284 - (-x223 - x278)*(-x224 - x281)) + V{0.00462962962962963}) - x153*x301*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x244 + x302) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x252 + x302) + cell.template getFieldComponent<descriptors::FORCE>(2)*(-x243 - x287));
auto x12 = V{1}*cell[12]*x119 - x102*(x300*(x150 + x263 - (-x117 - x258)*(-x214 - x260)) + V{0.00462962962962963}) - x153*x301*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x272 + x303) + cell.template getFieldComponent<descriptors::FORCE>(1)*(-x267 - x271) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x276 + x303));
auto x13 = -x102*(x300*(x292 + x323 - (-x304 - x306 - x308)*(-x309 - x311 - x313)) + V{0.00462962962962963}) + x119*x71 + x153*x301*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x192 + x302 + x303) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x267 + x302) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x287 + x303));
auto x14 = V{1}*cell[14]*x119 - x102*(x154*(-x304*x309 + x325) + V{0.0740740740740741}) - x153*x193*(-cell.template getFieldComponent<descriptors::FORCE>(0)*x326 + x152);
auto x15 = V{1}*cell[15]*x119 - x102*(x154*(-x315*x327 + x321) + V{0.0740740740740741}) - x153*x193*(-cell.template getFieldComponent<descriptors::FORCE>(1)*x266 + x210);
auto x16 = V{1}*cell[16]*x119 - x102*(x154*(-x317*x328 + x319 + x322) + V{0.0740740740740741}) - x153*x193*(-cell.template getFieldComponent<descriptors::FORCE>(2)*x286 + x221);
auto x17 = V{1}*cell[17]*x119 - x102*(x222*(-x329*x330 + x331) + V{0.0185185185185185}) + x153*x256*(cell.template getFieldComponent<descriptors::FORCE>(0)*x332 + cell.template getFieldComponent<descriptors::FORCE>(1)*x334 - x255);
auto x18 = V{1}*cell[18]*x119 - x102*(x222*(x201 - x258*x260 + x335) + V{0.0185185185185185}) - x153*x256*(-cell.template getFieldComponent<descriptors::FORCE>(0)*x337 + cell.template getFieldComponent<descriptors::FORCE>(1)*(x298 + x336) + x151);
auto x19 = V{1}*cell[19]*x119 - x102*(x222*(-x338*x339 + x340) + V{0.0185185185185185}) + x153*x256*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x302 + x326) + cell.template getFieldComponent<descriptors::FORCE>(2)*x341 - x275);
auto x20 = V{1}*cell[20]*x119 - x102*(x222*(x215 - x279*x282 + x335) + V{0.0185185185185185}) - x153*x256*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x271 + x326) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x336 + x342) + x146);
auto x21 = V{1}*cell[21]*x119 - x102*(x222*(x323 - (x315 + x317)*(x327 + x328)) + V{0.0185185185185185}) + x153*x256*(cell.template getFieldComponent<descriptors::FORCE>(1)*(x266 + x302) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x286 + x303) - x293);
auto x22 = V{1}*cell[22]*x119 - x102*(x222*(x215 + x262 - x296*x297) + V{0.0185185185185185}) - x153*x256*(-cell.template getFieldComponent<descriptors::FORCE>(1)*(x266 + x271) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x265 + x342) + x138);
auto x23 = V{1}*cell[23]*x119 - x102*(x300*(x322 + x331 - (x317 + x329)*(x328 + x330)) + V{0.00462962962962963}) + x153*x301*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x302 + x332) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x302 + x334) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x303 + x341));
auto x24 = -x102*(x300*(x254 + x331 - (x306 + x329)*(x311 + x330)) + V{0.00462962962962963}) + x119*x76 + x153*x301*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x271 + x332) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x271 + x334) - cell.template getFieldComponent<descriptors::FORCE>(2)*(x220 + x303 + x333));
auto x25 = -x102*(x300*(x274 + x340 - (x308 + x338)*(x313 + x339)) + V{0.00462962962962963}) + x119*x66 + x153*x301*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x302 + x337) - cell.template getFieldComponent<descriptors::FORCE>(1)*(x209 + x302 + x333) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x299 + x333));
auto x26 = V{1}*cell[26]*x119 - x102*(x300*(x290 + x335 - (x103 + x107 - x288)*(x155 + x171 - x289)) + V{0.00462962962962963}) - x153*x301*(cell.template getFieldComponent<descriptors::FORCE>(0)*(-x271 - x337) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x294 + x333) + cell.template getFieldComponent<descriptors::FORCE>(2)*(x295 + x333));
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
return { x33, x314 + x316 + x318 };
}
};

}

}

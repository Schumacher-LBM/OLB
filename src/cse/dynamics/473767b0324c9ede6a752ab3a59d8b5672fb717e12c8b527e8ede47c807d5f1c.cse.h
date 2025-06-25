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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[10] + cell[14];
auto x22 = cell[12] + cell[7];
auto x23 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x21 + x22;
auto x24 = x23 + V{1};
auto x25 = x23 + V{1};
auto x26 = V{1} / ((x25)*(x25));
auto x27 = V{1.5}*x26;
auto x28 = cell.template getFieldComponent<descriptors::POROSITY>(0)*cell.template getFieldComponent<descriptors::POROSITY>(0);
auto x29 = cell[13] - cell[4];
auto x30 = cell[15] - cell[6];
auto x31 = x29 + x30;
auto x32 = -cell[1];
auto x33 = cell[16] - cell[7];
auto x34 = x32 + x33;
auto x35 = -cell[5] + x21;
auto x36 = x31 + x34 + x35;
auto x37 = x36*x36;
auto x38 = x28*x37;
auto x39 = x27*x38;
auto x40 = cell[17] - cell[8];
auto x41 = x29 + x40;
auto x42 = cell[18] - cell[9];
auto x43 = -cell[2];
auto x44 = cell[11] - cell[14] + cell[5] + x43;
auto x45 = x41 + x42 + x44;
auto x46 = x45*x45;
auto x47 = x28*x46;
auto x48 = x27*x47;
auto x49 = x30 + x40;
auto x50 = -cell[3];
auto x51 = -cell[18] + cell[9];
auto x52 = x50 + x51;
auto x53 = -cell[16] + x22;
auto x54 = x49 + x52 + x53;
auto x55 = x54*x54;
auto x56 = x28*x55;
auto x57 = x27*x56;
auto x58 = x48 + x57 + V{-1};
auto x59 = x39 + x58;
auto x60 = V{0.0555555555555556}*x19;
auto x61 = V{3}*cell.template getFieldComponent<descriptors::POROSITY>(0)/x25;
auto x62 = x36*x61;
auto x63 = V{3}*x26;
auto x64 = x38*x63;
auto x65 = x45*x61;
auto x66 = x47*x63;
auto x67 = x39 + V{-1};
auto x68 = x54*x61;
auto x69 = x56*x63;
auto x70 = V{0.0277777777777778}*x19;
auto x71 = V{4.5}*x26;
auto x72 = cell[10] + x34;
auto x73 = cell[11] + V{2}*cell[13] - V{2}*cell[4] + x42 + x43 + x49 + x72;
auto x74 = x28*x71*(x73*x73);
auto x75 = x59 + x62;
auto x76 = -x65;
auto x77 = -cell[17] + cell[8];
auto x78 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x30 + x51 + x72 + x77;
auto x79 = -x78;
auto x80 = x32 + x35;
auto x81 = cell[12] + V{2}*cell[15] - V{2}*cell[6] + x41 + x52 + x80;
auto x82 = x28*x71*(x81*x81);
auto x83 = -x68;
auto x84 = -cell[12] + cell[3] + x29;
auto x85 = V{2}*cell[16] - V{2}*cell[7] + x42 + x77 + x80 + x84;
auto x86 = -x85;
auto x87 = V{2}*cell[17] - V{2}*cell[8] + x31 + x44 + x50 + x53;
auto x88 = x28*x71*(x87*x87);
auto x89 = x59 + x65;
auto x90 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x33 + x44 + x84;
auto x91 = -x90;
auto x92 = -x48;
auto x93 = V{1} - x57;
auto x94 = x92 + x93;
auto x95 = x62 + x94;
auto x96 = -x39;
auto x97 = x65 + x96;
auto x98 = x68 + x96;
auto x99 = -x62;
auto x100 = x59 + x68;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x24*x59 + V{1});
auto x1 = -cell[1]*x20 - x60*(x24*(x58 + x62 - x64) + V{1});
auto x2 = -cell[2]*x20 - x60*(x24*(x57 + x65 - x66 + x67) + V{1});
auto x3 = -cell[3]*x20 - x60*(x24*(x48 + x67 + x68 - x69) + V{1});
auto x4 = -cell[4]*x20 - x70*(x24*(x65 - x74 + x75) + V{1});
auto x5 = -(cell[5]*x20 + x70*(x24*(-x28*x71*x79*x79 + x75 + x76) + V{1}));
auto x6 = -cell[6]*x20 - x70*(x24*(x68 + x75 - x82) + V{1});
auto x7 = -(cell[7]*x20 + x70*(x24*(-x28*x71*x86*x86 + x75 + x83) + V{1}));
auto x8 = -cell[8]*x20 - x70*(x24*(x68 - x88 + x89) + V{1});
auto x9 = -(cell[9]*x20 + x70*(x24*(-x28*x71*x91*x91 + x83 + x89) + V{1}));
auto x10 = -cell[10]*x20 + x60*(x24*(x64 + x95) + V{-1});
auto x11 = -cell[11]*x20 + x60*(x24*(x66 + x93 + x97) + V{-1});
auto x12 = -cell[12]*x20 + x60*(x24*(x69 + x92 + x98 + V{1}) + V{-1});
auto x13 = -cell[13]*x20 + x70*(x24*(x74 + x95 + x97) + V{-1});
auto x14 = -(cell[14]*x20 + x70*(x24*(-x28*x71*x78*x78 + x89 + x99) + V{1}));
auto x15 = -cell[15]*x20 + x70*(x24*(x82 + x95 + x98) + V{-1});
auto x16 = -(cell[16]*x20 + x70*(x24*(x100 - x28*x71*x85*x85 + x99) + V{1}));
auto x17 = -cell[17]*x20 + x70*(x24*(x68 + x88 + x94 + x97) + V{-1});
auto x18 = -(cell[18]*x20 + x70*(x24*(x100 - x28*x71*x90*x90 + x76) + V{1}));
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
return { x25, V{1}*x26*x28*(x37 + x46 + x55) };
}
};

}

}

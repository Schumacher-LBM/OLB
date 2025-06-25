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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::ShanChen>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x22 = x21 + V{1};
auto x23 = V{1} / (x19);
auto x24 = cell.template getFieldComponent<descriptors::FORCE>(0)*x23;
auto x25 = cell.template getFieldComponent<descriptors::VELOCITY>(0) + x24;
auto x26 = x25*x25;
auto x27 = V{1.5}*x26;
auto x28 = cell.template getFieldComponent<descriptors::FORCE>(1)*x23;
auto x29 = cell.template getFieldComponent<descriptors::VELOCITY>(1) + x28;
auto x30 = x29*x29;
auto x31 = V{1.5}*x30;
auto x32 = cell.template getFieldComponent<descriptors::FORCE>(2)*x23;
auto x33 = cell.template getFieldComponent<descriptors::VELOCITY>(2) + x32;
auto x34 = x33*x33;
auto x35 = V{1.5}*x34;
auto x36 = x31 + x35 + V{-1};
auto x37 = x27 + x36;
auto x38 = V{0.0555555555555556}*x19;
auto x39 = V{3}*x26;
auto x40 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x41 = V{3}*x24;
auto x42 = x40 + x41;
auto x43 = V{3}*x30;
auto x44 = x27 + V{-1};
auto x45 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x46 = V{3}*x28;
auto x47 = x45 + x46;
auto x48 = V{3}*x34;
auto x49 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x50 = V{3}*x32;
auto x51 = x49 + x50;
auto x52 = V{0.0277777777777778}*x19;
auto x53 = x25 + x29;
auto x54 = V{4.5}*(x53*x53);
auto x55 = x37 + x42;
auto x56 = -cell.template getFieldComponent<descriptors::VELOCITY>(1) + x25 - x28;
auto x57 = -x56;
auto x58 = -x45 - x46;
auto x59 = x25 + x33;
auto x60 = V{4.5}*(x59*x59);
auto x61 = -cell.template getFieldComponent<descriptors::FORCE>(2)*x23 - cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x62 = x25 + x61;
auto x63 = -x62;
auto x64 = -x49 - x50;
auto x65 = x29 + x33;
auto x66 = V{4.5}*(x65*x65);
auto x67 = x37 + x47;
auto x68 = x29 + x61;
auto x69 = -x68;
auto x70 = -x31;
auto x71 = V{1} - x35;
auto x72 = x70 + x71;
auto x73 = x42 + x72;
auto x74 = -x27;
auto x75 = x47 + x74;
auto x76 = x51 + x74;
auto x77 = -x40 - x41;
auto x78 = x37 + x51;
auto x79 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x80 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1) + cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x81 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(2) + cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x22*x37 + V{1});
auto x1 = -cell[1]*x20 - x38*(x22*(x36 - x39 + x42) + V{1});
auto x2 = -cell[2]*x20 - x38*(x22*(x35 - x43 + x44 + x47) + V{1});
auto x3 = -cell[3]*x20 - x38*(x22*(x31 + x44 - x48 + x51) + V{1});
auto x4 = -cell[4]*x20 - x52*(x22*(x47 - x54 + x55) + V{1});
auto x5 = -(cell[5]*x20 + x52*(x22*(x55 + x58 - V{4.5}*x57*x57) + V{1}));
auto x6 = -cell[6]*x20 - x52*(x22*(x51 + x55 - x60) + V{1});
auto x7 = -(cell[7]*x20 + x52*(x22*(x55 + x64 - V{4.5}*x63*x63) + V{1}));
auto x8 = -cell[8]*x20 - x52*(x22*(x51 - x66 + x67) + V{1});
auto x9 = -(cell[9]*x20 + x52*(x22*(x64 + x67 - V{4.5}*x69*x69) + V{1}));
auto x10 = -cell[10]*x20 + x38*(x22*(x39 + x73) + V{-1});
auto x11 = -cell[11]*x20 + x38*(x22*(x43 + x71 + x75) + V{-1});
auto x12 = -cell[12]*x20 + x38*(x22*(x48 + x70 + x76 + V{1}) + V{-1});
auto x13 = -cell[13]*x20 + x52*(x22*(x54 + x73 + x75) + V{-1});
auto x14 = -(cell[14]*x20 + x52*(x22*(x67 + x77 - V{4.5}*x56*x56) + V{1}));
auto x15 = -cell[15]*x20 + x52*(x22*(x60 + x73 + x76) + V{-1});
auto x16 = -(cell[16]*x20 + x52*(x22*(x77 + x78 - V{4.5}*x62*x62) + V{1}));
auto x17 = -cell[17]*x20 + x52*(x22*(x51 + x66 + x72 + x75) + V{-1});
auto x18 = -(cell[18]*x20 + x52*(x22*(x58 + x78 - V{4.5}*x68*x68) + V{1}));
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
return { x21 + V{1}, x79*x79 + x80*x80 + x81*x81 };
}
};

}

}

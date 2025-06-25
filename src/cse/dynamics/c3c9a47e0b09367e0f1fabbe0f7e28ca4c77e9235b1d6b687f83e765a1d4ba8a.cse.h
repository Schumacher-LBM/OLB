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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::ForcedWithStress> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x11 = cell[2] + cell[3];
auto x12 = cell[7] + cell[8];
auto x13 = cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + x11 + x12;
auto x14 = x13 + V{1};
auto x15 = V{1} / (x14);
auto x16 = cell[1] - cell[5];
auto x17 = -cell[6] - cell[7] + x11 + x16;
auto x18 = x15*(x13 + V{1});
auto x19 = -V{0.333333333333333}*cell[0] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[7];
auto x20 = -cell.template getFieldComponent<descriptors::FORCE>(0)*x17*x18 + V{0.666666666666667}*cell[2] - V{0.333333333333333}*cell[4] + V{0.666666666666667}*cell[6] - V{0.333333333333333}*cell[8] - x15*x17*x17 + x19;
auto x21 = -cell[3] - cell[4] + x12 + x16;
auto x22 = cell.template getFieldComponent<descriptors::FORCE>(1)*x18*x21 - V{0.333333333333333}*cell[2] + V{0.666666666666667}*cell[4] - V{0.333333333333333}*cell[6] + V{0.666666666666667}*cell[8] - x15*x21*x21 + x19;
auto x23 = V{1}*cell[7];
auto x24 = V{1}*cell[1];
auto x25 = x18*(cell.template getFieldComponent<descriptors::FORCE>(0)*x21 - cell.template getFieldComponent<descriptors::FORCE>(1)*x17);
auto x26 = x15*x17*x21;
auto x27 = V{1}*cell[3];
auto x28 = -V{1}*cell[5];
auto x29 = x27 + x28;
auto x30 = V{1} / (V{2.52268963608289}*util::sqrt(x15*(x10*x10)*util::sqrt((-x23 + x24 - V{0.5}*x25 - V{1}*x26 - x29)*(V{2}*cell[1] - V{2}*cell[3] + V{2}*cell[5] - V{2}*cell[7] - V{1}*x25 - V{2}*x26) + V{1}*(x20*x20) + V{1}*(x22*x22)) + V{0.0392836979096202}/((x9)*(x9))) + V{0.5}/x9);
auto x31 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x32 = V{1}*cell[2] - V{1}*cell[6] - x23 + x24 + x29;
auto x33 = -x15*x32 + x31;
auto x34 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x35 = V{1}*cell[8];
auto x36 = x15*(-V{1}*cell[4] + x23 + x24 - x27 + x28 + x35);
auto x37 = x34 + x36;
auto x38 = x37*x37;
auto x39 = V{1.5}*x38;
auto x40 = x39 + V{-1};
auto x41 = x40 + V{1.5}*(x33*x33);
auto x42 = V{1} - x30;
auto x43 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x44 = V{3}*cell[3];
auto x45 = V{3}*cell[1] - V{3}*cell[5];
auto x46 = V{3}*cell[2] - V{3}*cell[6] - V{3}*cell[7] + x44 + x45;
auto x47 = -x15*x46 + x43;
auto x48 = cell.template getFieldComponent<descriptors::FORCE>(0)*x47;
auto x49 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x50 = x15*(-V{3}*cell[4] + V{3}*cell[7] + V{3}*cell[8] - x44 + x45);
auto x51 = x49 + x50;
auto x52 = cell.template getFieldComponent<descriptors::FORCE>(1)*x51;
auto x53 = V{1} - V{0.5}*x30;
auto x54 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x55 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x56 = V{4.5}*cell[3];
auto x57 = V{4.5}*cell[1] - V{4.5}*cell[5];
auto x58 = x15*(-V{4.5}*cell[4] + V{4.5}*cell[7] + V{4.5}*cell[8] - x56 + x57);
auto x59 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x60 = V{4.5}*cell[2] - V{4.5}*cell[6] - V{4.5}*cell[7] + x56 + x57;
auto x61 = -x15*x60 + x59;
auto x62 = x41 + x47;
auto x63 = V{6}*cell[3];
auto x64 = V{6}*cell[1] - V{6}*cell[5];
auto x65 = V{6}*cell[2] - V{6}*cell[6] - V{6}*cell[7] + x63 + x64;
auto x66 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x67 = -x66;
auto x68 = x67 + V{3};
auto x69 = V{9}*cell[7];
auto x70 = V{9}*cell[3];
auto x71 = V{9}*cell[1] - V{9}*cell[5];
auto x72 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1) + x15*(-V{9}*cell[4] + V{9}*cell[8] + x69 - x70 + x71);
auto x73 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1) + x15*(-V{6}*cell[4] + V{6}*cell[7] + V{6}*cell[8] - x63 + x64);
auto x74 = x73 + V{3};
auto x75 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x76 = x15*(V{9}*cell[2] - V{9}*cell[6] - x69 + x70 + x71);
auto x77 = -x75 + x76;
auto x78 = V{0.0277777777777778}*x14;
auto x79 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x80 = x15*x65;
auto x81 = V{0.111111111111111}*x14;
auto x82 = x55 + x58;
auto x83 = x66 - x80;
auto x84 = x72 + V{-3};
auto x85 = x73 + V{-3};
auto x86 = x75 - x76;
auto x87 = x37*x82;
auto x88 = x15*x32;
auto x89 = x31 - x88;
auto x90 = x89*x89;
auto x91 = V{1.5}*x90;
auto x92 = x15*x46;
auto x93 = x15*x60;
auto x94 = x59 - x93;
auto x95 = -x39 - x91 + V{1};
auto x96 = x43 - x92 + x95;
auto x97 = x83 + V{3};
auto x0 = V{1}*cell[0]*x42 - V{0.444444444444444}*x14*x53*(x48 + x52) - x30*(x41*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) + V{0.444444444444444});
auto x1 = V{1}*cell[1]*x42 - x30*(x54*(-x49 - x50 + x62 - (-x33 + x34 + x36)*(x55 + x58 - x61)) + V{0.0277777777777778}) - x53*x78*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x15*x65 + x68 + x72) - cell.template getFieldComponent<descriptors::FORCE>(1)*(x74 + x77));
auto x2 = V{1}*cell[2]*x42 - x30*(x79*(-x33*x61 + x62) + V{0.111111111111111}) - x53*x81*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x68 + x80) + x52);
auto x3 = x27*x42 - x30*(x54*(x51 + x62 - (x33 + x37)*(x61 + x82)) + V{0.0277777777777778}) + x53*x78*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x83 + x84) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x85 + x86));
auto x4 = V{1}*cell[4]*x42 - x30*(x79*(x41 + x51 - x87) + V{0.111111111111111}) - x53*x81*(-cell.template getFieldComponent<descriptors::FORCE>(1)*x85 + x48);
auto x5 = V{1}*cell[5]*x42 - x30*(x54*(x40 - x43 + x51 + x91 + x92 - (x31 - x37 - x88)*(x59 - x82 - x93)) + V{0.0277777777777778}) - x53*x78*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x67 + x80 + x84) + cell.template getFieldComponent<descriptors::FORCE>(1)*(-x77 - x85));
auto x6 = V{1}*cell[6]*x42 + x30*(x79*(x89*x94 + x96) + V{-0.111111111111111}) + x53*x81*(cell.template getFieldComponent<descriptors::FORCE>(0)*x97 - x52);
auto x7 = x23*x42 + x30*(x54*(x51 + x96 + (x37 + x89)*(x82 + x94)) + V{-0.0277777777777778}) + x53*x78*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x72 + x97) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x74 + x86));
auto x8 = x30*(x79*(x51 + x87 + x95) + V{-0.111111111111111}) + x35*x42 - x53*x81*(-cell.template getFieldComponent<descriptors::FORCE>(1)*x74 + x48);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x14, x38 + x90 };
}
};

}

}

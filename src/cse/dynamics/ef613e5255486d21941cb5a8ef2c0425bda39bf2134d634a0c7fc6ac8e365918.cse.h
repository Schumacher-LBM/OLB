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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, guoZhao::GuoZhaoSecondOrder, collision::BGK, guoZhao::GuoZhaoForcing<momenta::GuoZhaoForced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x12 = x11 + V{1};
auto x13 = V{1} / (cell.template getFieldComponent<descriptors::EPSILON>(0));
auto x14 = V{0.5}*cell.template getFieldComponent<descriptors::EPSILON>(0);
auto x15 = cell.template getFieldComponent<descriptors::NU>(0)/cell.template getFieldComponent<descriptors::K>(0);
auto x16 = x14*x15 + V{1};
auto x17 = util::sqrt(x16*x16);
auto x18 = V{1} / ((x16 + x17)*(x16 + x17));
auto x19 = x11 + V{1};
auto x20 = V{1} / (x19);
auto x21 = V{1}*cell[7];
auto x22 = V{1}*cell[3];
auto x23 = V{1}*cell[1] - V{1}*cell[5];
auto x24 = -V{1}*cell[4] + V{1}*cell[8] + x21 - x22 + x23;
auto x25 = x20*x24;
auto x26 = cell.template getFieldComponent<descriptors::BODY_FORCE>(1)*x14;
auto x27 = x25 + x26;
auto x28 = x27*x27;
auto x29 = V{6}*x28;
auto x30 = V{1}*cell[2] - V{1}*cell[6] - x21 + x22 + x23;
auto x31 = x20*x30;
auto x32 = cell.template getFieldComponent<descriptors::BODY_FORCE>(0)*x14;
auto x33 = -x31 + x32;
auto x34 = x33*x33;
auto x35 = -x13*x18*(x29 + V{6}*x34) + V{1};
auto x36 = x19*(V{0.5}*x9 + V{-1});
auto x37 = V{1} / (V{0.25}*cell.template getFieldComponent<descriptors::EPSILON>(0)*x15 + V{0.5}*x17 + V{0.5});
auto x38 = V{1}*x15;
auto x39 = cell.template getFieldComponent<descriptors::BODY_FORCE>(1) - x27*x37*x38;
auto x40 = cell.template getFieldComponent<descriptors::BODY_FORCE>(0) - x33*x37*x38;
auto x41 = V{1.5}*cell.template getFieldComponent<descriptors::EPSILON>(0);
auto x42 = cell.template getFieldComponent<descriptors::BODY_FORCE>(0)*x41;
auto x43 = -V{3}*x31 + x42;
auto x44 = x37*x43;
auto x45 = V{18}*x13;
auto x46 = -x20*x30;
auto x47 = x32 + x46;
auto x48 = -x25 - x26 + x47;
auto x49 = cell.template getFieldComponent<descriptors::BODY_FORCE>(1)*x41 + V{3}*x25;
auto x50 = x37*x49;
auto x51 = x35 + x50;
auto x52 = V{0.0277777777777778}*x36;
auto x53 = V{4.5}*cell.template getFieldComponent<descriptors::EPSILON>(0);
auto x54 = cell.template getFieldComponent<descriptors::BODY_FORCE>(1)*x53;
auto x55 = V{9}*x25;
auto x56 = cell.template getFieldComponent<descriptors::BODY_FORCE>(0)*x53;
auto x57 = V{9}*x46 + x56;
auto x58 = x37*(-x54 - x55 + x57);
auto x59 = V{3}*cell.template getFieldComponent<descriptors::EPSILON>(0);
auto x60 = x42 + V{3}*x46;
auto x61 = x37*x60 + x59;
auto x62 = x37*x49;
auto x63 = V{0.111111111111111}*x9;
auto x64 = x47*x47;
auto x65 = x13*x18*(x29 + V{6}*x64) + V{-1};
auto x66 = x37*x60 + x65;
auto x67 = V{9}*x31;
auto x68 = x56 - x67;
auto x69 = -x37*x39*(V{0.166666666666667}*cell.template getFieldComponent<descriptors::BODY_FORCE>(1)*cell.template getFieldComponent<descriptors::EPSILON>(0) + V{0.333333333333333}*x20*x24);
auto x70 = x27 + x47;
auto x71 = x54 + x55;
auto x72 = -x37*(x57 + x71);
auto x73 = x59 + x62;
auto x74 = x18*x28*x45;
auto x75 = -x37*x40*(V{0.166666666666667}*cell.template getFieldComponent<descriptors::BODY_FORCE>(0)*cell.template getFieldComponent<descriptors::EPSILON>(0) - V{0.333333333333333}*x20*x30);
auto x76 = V{0.5}*cell.template getFieldComponent<descriptors::BODY_FORCE>(0)*cell.template getFieldComponent<descriptors::EPSILON>(0) - x27 - x31;
auto x77 = x35 + x44;
auto x78 = x37*(V{4.5}*cell.template getFieldComponent<descriptors::BODY_FORCE>(0)*cell.template getFieldComponent<descriptors::EPSILON>(0) - x67 - x71);
auto x79 = -x37*x43 + x59;
auto x80 = x27 + x33;
auto x81 = x37*(x68 + x71);
auto x82 = x59 - x62;
auto x0 = -cell[0]*x10 + V{1.33333333333333}*x36*x37*(x27*x39 + x33*x40) + V{0.444444444444444}*x9*(x12*x35 + V{-1});
auto x1 = -(cell[1]*x10 + x52*(x39*(-x58 + x59 - x62) + x40*(x58 - x61)) - V{0.0277777777777778}*x9*(x12*(x18*x45*(x48*x48) - x44 + x51) + V{-1}));
auto x2 = -cell[2]*x10 - x36*(V{0.111111111111111}*x40*(x37*x68 - x61) + x69) - x63*(x12*(-x18*x45*x64 + x66) + V{1});
auto x3 = -(cell[3]*x10 + x52*(x39*(-x72 - x73) + x40*(-x61 - x72)) + V{0.0277777777777778}*x9*(x12*(-x18*x45*x70*x70 + x50 + x66) + V{1}));
auto x4 = -cell[4]*x10 - x36*(V{0.111111111111111}*x39*(x37*x71 - x73) + x75) - x63*(x12*(x50 + x65 - x74) + V{1});
auto x5 = -(cell[5]*x10 + x52*(-x39*(x73 + x78) + x40*(x78 + x79)) - V{0.0277777777777778}*x9*(x12*(x18*x45*(x76*x76) - x50 + x77) + V{-1}));
auto x6 = -cell[6]*x10 - x36*(V{0.111111111111111}*x40*(x37*x68 + x79) + x69) + V{0.111111111111111}*x9*(x12*(x18*x34*x45 + x77) + V{-1});
auto x7 = -(cell[7]*x10 + x52*(x39*(x81 + x82) + x40*(x79 + x81)) - V{0.0277777777777778}*x9*(x12*(x18*x45*(x80*x80) + x50 + x77) + V{-1}));
auto x8 = -cell[8]*x10 - x36*(V{0.111111111111111}*x39*(x37*x71 + x82) + x75) + V{0.111111111111111}*x9*(x12*(x51 + x74) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x19, V{4}*x18*(x28 + x34) };
}
};

}

}

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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::Guo<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x12 = x11 + V{1};
auto x13 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x14 = x11 + V{1};
auto x15 = V{1} / (x14);
auto x16 = V{1}*cell[3];
auto x17 = V{1}*cell[1] - V{1}*cell[5];
auto x18 = V{1}*cell[2] - V{1}*cell[6] - V{1}*cell[7] + x16 + x17;
auto x19 = x13 - x15*x18;
auto x20 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x21 = x15*(-V{1}*cell[4] + V{1}*cell[7] + V{1}*cell[8] - x16 + x17);
auto x22 = x20 + x21;
auto x23 = x22*x22;
auto x24 = V{1.5}*x23;
auto x25 = x24 + V{-1};
auto x26 = x25 + V{1.5}*(x19*x19);
auto x27 = V{0.5}*x9 + V{-1};
auto x28 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x29 = V{3}*cell[3];
auto x30 = V{3}*cell[1] - V{3}*cell[5];
auto x31 = V{3}*cell[2] - V{3}*cell[6] - V{3}*cell[7] + x29 + x30;
auto x32 = -x15*x31 + x28;
auto x33 = cell.template getFieldComponent<descriptors::FORCE>(0)*x32;
auto x34 = V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x35 = x15*(-V{3}*cell[4] + V{3}*cell[7] + V{3}*cell[8] - x29 + x30);
auto x36 = x34 + x35;
auto x37 = cell.template getFieldComponent<descriptors::FORCE>(1)*x36;
auto x38 = V{0.0277777777777778}*x9;
auto x39 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x40 = V{4.5}*cell[3];
auto x41 = V{4.5}*cell[1] - V{4.5}*cell[5];
auto x42 = x15*(-V{4.5}*cell[4] + V{4.5}*cell[7] + V{4.5}*cell[8] - x40 + x41);
auto x43 = V{2.25}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x44 = V{4.5}*cell[2] - V{4.5}*cell[6] - V{4.5}*cell[7] + x40 + x41;
auto x45 = -x15*x44 + x43;
auto x46 = x26 + x32;
auto x47 = V{6}*cell[3];
auto x48 = V{6}*cell[1] - V{6}*cell[5];
auto x49 = V{6}*cell[2] - V{6}*cell[6] - V{6}*cell[7] + x47 + x48;
auto x50 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x51 = -x50;
auto x52 = x51 + V{3};
auto x53 = V{9}*cell[7];
auto x54 = V{9}*cell[3];
auto x55 = V{9}*cell[1] - V{9}*cell[5];
auto x56 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(1) + x15*(-V{9}*cell[4] + V{9}*cell[8] + x53 - x54 + x55);
auto x57 = V{4.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x58 = V{9}*cell[2] - V{9}*cell[6] - x53 + x54 + x55;
auto x59 = x15*x58;
auto x60 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x61 = x15*(-V{6}*cell[4] + V{6}*cell[7] + V{6}*cell[8] - x47 + x48);
auto x62 = x60 + x61;
auto x63 = x62 + V{3};
auto x64 = x14*x27;
auto x65 = V{0.0277777777777778}*x64;
auto x66 = V{0.111111111111111}*x9;
auto x67 = x15*x49;
auto x68 = -x37;
auto x69 = V{0.111111111111111}*x64;
auto x70 = x39 + x42;
auto x71 = x50 - x67;
auto x72 = x56 + V{-3};
auto x73 = x62 + V{-3};
auto x74 = x57 - x59;
auto x75 = x22*x70;
auto x76 = x15*x18;
auto x77 = x13 - x76;
auto x78 = x77*x77;
auto x79 = V{1.5}*x78;
auto x80 = x15*x31;
auto x81 = x15*x44;
auto x82 = x43 - x81;
auto x83 = -x24 - x79 + V{1};
auto x84 = x28 - x80 + x83;
auto x85 = x71 + V{3};
auto x0 = -cell[0]*x10 + V{0.444444444444444}*x14*x27*(x33 + x37) - V{0.444444444444444}*x9*(x12*x26 + V{1});
auto x1 = -cell[1]*x10 - x38*(x12*(-x34 - x35 + x46 - (-x19 + x20 + x21)*(x39 + x42 - x45)) + V{1}) - x65*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x15*x49 + x52 + x56) + cell.template getFieldComponent<descriptors::FORCE>(1)*(-x57 + x59 + x63));
auto x2 = -cell[2]*x10 - x66*(x12*(-x19*x45 + x46) + V{1}) - x69*(cell.template getFieldComponent<descriptors::FORCE>(0)*(-x52 - x67) + x68);
auto x3 = -cell[3]*x10 - x38*(x12*(x36 + x46 - (x19 + x22)*(x45 + x70)) + V{1}) - x65*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x71 + x72) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x73 + x74));
auto x4 = -cell[4]*x10 - x66*(x12*(x26 + x36 - x75) + V{1}) - x69*(cell.template getFieldComponent<descriptors::FORCE>(1)*x73 - x33);
auto x5 = -cell[5]*x10 - x65*(cell.template getFieldComponent<descriptors::FORCE>(0)*(-x51 - x67 - x72) - cell.template getFieldComponent<descriptors::FORCE>(1)*(-x15*x58 + x57 - x60 - x61 + V{3})) + V{0.0277777777777778}*x9*(x12*(V{1.5}*cell.template getFieldComponent<descriptors::FORCE>(0) - x25 - x36 - x79 - x80 + (x13 - x22 - x76)*(x43 - x70 - x81)) + V{-1});
auto x6 = -cell[6]*x10 - x69*(cell.template getFieldComponent<descriptors::FORCE>(0)*x85 + x68) + V{0.111111111111111}*x9*(x12*(x77*x82 + x84) + V{-1});
auto x7 = -cell[7]*x10 - x65*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x56 + x85) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x63 + x74)) + V{0.0277777777777778}*x9*(x12*(x36 + x84 + (x22 + x77)*(x70 + x82)) + V{-1});
auto x8 = -cell[8]*x10 - x69*(cell.template getFieldComponent<descriptors::FORCE>(1)*x63 - x33) + V{0.111111111111111}*x9*(x12*(x36 + x75 + x83) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x14, x23 + x78 };
}
};

}

}

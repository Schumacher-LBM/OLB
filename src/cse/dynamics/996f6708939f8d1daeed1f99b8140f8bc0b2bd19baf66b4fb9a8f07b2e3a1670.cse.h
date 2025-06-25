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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>, forcing::Wagner>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-1};
auto x10 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x11 = x10 + V{1};
auto x12 = V{0.666666} - V{0.1666665}*cell.template getFieldComponent<descriptors::OMEGA>(0);
auto x13 = cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x14 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x15 = V{1} / (x11);
auto x16 = cell.template getFieldComponent<descriptors::OMEGA>(0)*cell.template getFieldComponent<descriptors::SCALAR>(0)*x15;
auto x17 = cell[1] - cell[5];
auto x18 = cell[2] - cell[6];
auto x19 = cell[3] - cell[7] + x17 + x18;
auto x20 = -x19;
auto x21 = V{1.33333333333333}*x15;
auto x22 = -cell[4] + cell[8];
auto x23 = -cell[3] + cell[7] + x17 + x22;
auto x24 = cell.template getFieldComponent<descriptors::FORCE>(1)*x23;
auto x25 = x10 + V{1};
auto x26 = V{1} / ((x11)*(x11));
auto x27 = V{1.5}*x26;
auto x28 = x23*x23;
auto x29 = x27*x28;
auto x30 = x27*(x20*x20) + V{-1};
auto x31 = x29 + x30;
auto x32 = V{6}*cell[3];
auto x33 = V{6}*cell[1] - V{6}*cell[5];
auto x34 = V{6}*cell[2] - V{6}*cell[6] - V{6}*cell[7] + x32 + x33;
auto x35 = V{9}*cell[7];
auto x36 = V{9}*cell[3];
auto x37 = V{9}*cell[1] - V{9}*cell[5];
auto x38 = x15*(-V{9}*cell[4] + V{9}*cell[8] + x35 - x36 + x37);
auto x39 = x38 + V{3};
auto x40 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x41 = V{9}*cell[2] - V{9}*cell[6] - x35 + x36 + x37;
auto x42 = x15*x41;
auto x43 = x15*(-V{6}*cell[4] + V{6}*cell[7] + V{6}*cell[8] - x32 + x33);
auto x44 = x43 + V{3};
auto x45 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x46 = V{0.25}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-1};
auto x47 = V{0.25}*cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::FORCE>(1)*x46;
auto x48 = V{0.02083334375}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-0.083333375};
auto x49 = x13*x48 + x14*x48 - V{0.0138888958333333}*x16;
auto x50 = -x47 + x49;
auto x51 = V{3}*cell[3];
auto x52 = V{3}*cell[7];
auto x53 = V{3}*cell[1] - V{3}*cell[5];
auto x54 = V{3}*cell[2] - V{3}*cell[6] + x51 - x52 + x53;
auto x55 = x15*x54;
auto x56 = V{4.5}*x26;
auto x57 = V{2}*cell[1] - V{2}*cell[5] + x18 + x22;
auto x58 = x56*(x57*x57);
auto x59 = V{1} - x29;
auto x60 = x15*(-V{3}*cell[4] + V{3}*cell[8] - x51 + x52 + x53);
auto x61 = x19*x19;
auto x62 = x27*x61;
auto x63 = x60 - x62;
auto x64 = x59 + x63;
auto x65 = x15*x34;
auto x66 = x65 + V{3};
auto x67 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x68 = V{0.083333375}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-0.3333335};
auto x69 = V{0.333333333333333}*x15;
auto x70 = -V{0.0138889166666667}*cell.template getFieldComponent<descriptors::OMEGA>(0)*cell.template getFieldComponent<descriptors::SCALAR>(0)*x15 + x13*x68 - V{0.1666665}*x14*x46 + x24*x69;
auto x71 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::OMEGA>(0);
auto x72 = V{3}*x26;
auto x73 = x59 + x61*x72;
auto x74 = V{3} - x43;
auto x75 = x47 + x49;
auto x76 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::OMEGA>(0);
auto x77 = V{2}*cell[3] + cell[4] - V{2}*cell[7] - cell[8] + x18;
auto x78 = -x77;
auto x79 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x80 = cell.template getFieldComponent<descriptors::FORCE>(0)*x19*x69 + x13*(V{0.041666625}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-0.1666665}) - x14*x68 + V{0.0138889166666667}*x16;
auto x81 = x28*x72;
auto x82 = x65 + V{-3};
auto x83 = -x55;
auto x0 = -V{0.444444444444444}*cell.template getFieldComponent<descriptors::OMEGA>(0)*(x25*x31 + V{1}) - cell[0]*x9 - x11*(cell.template getFieldComponent<descriptors::FORCE>(0)*x20*x21 + x12*x13 + x12*x14 + V{0.111111}*x16 + x21*x24);
auto x1 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::OMEGA>(0)*(x25*(x55 + x58 + x64) + V{-1}) - cell[1]*x9 - x11*(x40*(x15*x34 + x39) - x45*(x42 + x44) + x50);
auto x2 = -cell[2]*x9 + x11*(-x66*x67 - x70) + x71*(x25*(x55 + x73) + V{-1});
auto x3 = -(cell[3]*x9 + x11*(-x40*(x38 - x66) - x45*(-x42 - x74) + x75) + x76*(x25*(-x15*x54 + x31 - x56*x78*x78 + x60) + V{1}));
auto x4 = -cell[4]*x9 + x11*(x79*(x43 + V{-3}) + x80) - x71*(x25*(x30 + x60 - x81) + V{1});
auto x5 = -cell[5]*x9 - x11*(-x40*(-x38 - x82) + x45*(-x15*x41 + x74) + x50) - x76*(x25*(x29 + x55 - x58 + x60 + x62 + V{-1}) + V{1});
auto x6 = -cell[6]*x9 + x11*(-x67*x82 - x70) + x71*(x25*(x73 + x83) + V{-1});
auto x7 = -(-V{0.0277777777777778}*cell.template getFieldComponent<descriptors::OMEGA>(0)*(x25*(x56*(x77*x77) + x64 + x83) + V{-1}) + cell[7]*x9 + x11*(-x40*(x39 - x65) - x45*(-x42 + x44) + x75));
auto x8 = -cell[8]*x9 + x11*(x44*x79 + x80) + x71*(x25*(x63 + x81 + V{1}) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x11, V{1}*x26*(x28 + x61) };
}
};

}

}

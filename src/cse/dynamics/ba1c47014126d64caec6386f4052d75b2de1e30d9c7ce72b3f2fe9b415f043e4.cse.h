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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>, forcing::Wagner>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = V{0.666666} - V{0.1666665}*cell.template getFieldComponent<descriptors::OMEGA>(0);
auto x10 = cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x11 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x12 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x13 = V{1} / (cell.template getFieldComponent<momenta::FixedDensity::RHO>(0));
auto x14 = V{0.111111111111111}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x15 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x16 = V{3}*x15 + V{1};
auto x17 = V{3}*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x18 = cell[0] + cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + V{1};
auto x19 = -x13*x18 + V{1};
auto x20 = -x19;
auto x21 = x20*x20;
auto x22 = V{1.5}*x21;
auto x23 = -x22;
auto x24 = x17 + x23;
auto x25 = x16 + x24;
auto x26 = -x17;
auto x27 = x16 + x26;
auto x28 = x23 + x27;
auto x29 = V{0.444444444444444}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x30 = V{1.5}*x15;
auto x31 = x30 + V{-1};
auto x32 = x22 + x31;
auto x33 = -x32;
auto x34 = V{0.222222222222222}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x35 = V{3}*x21;
auto x36 = x13*(V{3}*cell[0] + V{3}*cell[4] + V{6}*cell[5] + V{6}*cell[6] + V{6}*cell[7] + V{3}*cell[8] + V{3});
auto x37 = -x30 + x36 + V{-2};
auto x38 = x35 + x37;
auto x39 = V{0.0555555555555555}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x40 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1) + x13*x18 + V{-1};
auto x41 = V{4.5}*(x40*x40);
auto x42 = x24 + x37 + x41;
auto x43 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1) + x19;
auto x44 = -x43;
auto x45 = x30 + V{2};
auto x46 = -x36 + x45;
auto x47 = x22 + x46;
auto x48 = x17 + x47 - V{4.5}*x44*x44;
auto x49 = -x48;
auto x50 = V{1.66533453693773e-16}*cell[4] + V{1.33226762955019e-15}*cell[5] + V{6.66133814775094e-16}*cell[6] + V{8.88178419700125e-16}*cell[7] + V{1.66533453693773e-16}*cell[8] + V{1.11022302462516e-16};
auto x51 = x13*(x14*x25 + x14*x28 + x29*x33 + x34*x38 + x39*x42 + x39*x49 + x50);
auto x52 = V{1} - x51;
auto x53 = cell.template getFieldComponent<descriptors::OMEGA>(0)*cell.template getFieldComponent<descriptors::SCALAR>(0)*x13;
auto x54 = cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-1};
auto x55 = V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x56 = -x41;
auto x57 = x26 + x56;
auto x58 = x47 + x57;
auto x59 = -x35 + x46;
auto x60 = x29*x32;
auto x61 = x19*x19;
auto x62 = V{1.5}*x61;
auto x63 = -x62;
auto x64 = x16 + x17;
auto x65 = x63 + x64;
auto x66 = x27 + x63;
auto x67 = x51 + V{-1};
auto x68 = x67*x67;
auto x69 = V{1.5}*x68;
auto x70 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x71 = V{9}*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x72 = V{0.333333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x73 = V{1.33333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x74 = V{0.666666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x75 = V{2.66666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x76 = V{9.99200722162641e-16}*cell[4] + V{7.99360577730113e-15}*cell[5] + V{3.99680288865056e-15}*cell[6] + V{5.32907051820075e-15}*cell[7] + V{9.99200722162641e-16}*cell[8] + V{6.66133814775094e-16};
auto x77 = -x13*(x25*x74 + x28*x74 + x33*x75 + x38*x73 + x42*x72 + x49*x72 + x76);
auto x78 = x77 + V{9};
auto x79 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x80 = V{6}*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x81 = V{1}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x82 = V{0.5}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x83 = x13*(-V{4}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x32 - V{2}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x59 + V{1.49880108324396e-15}*cell[4] + V{1.19904086659517e-14}*cell[5] + V{5.99520433297585e-15}*cell[6] + V{7.99360577730113e-15}*cell[7] + V{1.49880108324396e-15}*cell[8] - x48*x82 - x58*x82 + x65*x81 + x66*x81 + V{9.99200722162641e-16});
auto x84 = x80 - x83;
auto x85 = V{0.25}*cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::FORCE>(1)*(V{0.25}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-1});
auto x86 = V{0.02083334375}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-0.083333375};
auto x87 = x10*x86 + x11*x86 - V{0.0138888958333333}*x53;
auto x88 = -x85 + x87;
auto x89 = V{0.0231481481481481}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x90 = x62 - V{4.5}*x43*x43;
auto x91 = x17 + x46 + x90;
auto x92 = V{0.0277777777777778}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x93 = x30 + V{-4};
auto x94 = x36 + x93;
auto x95 = V{0.00462962962962963}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x96 = x46 + x57 + x62;
auto x97 = -V{3}*x61;
auto x98 = x46 + x97;
auto x99 = V{0.00925925925925926}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x100 = V{0.0185185185185185}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x98 + V{0.0833333333333333}*cell[4] + V{0.166666666666667}*cell[6] + V{0.0833333333333333}*cell[8] - x65*x99 - x66*x99 + V{0.0277777777777778};
auto x101 = V{0.833333333333333}*cell[5] - V{0.166666666666667}*cell[7] + x100 - x95*x96;
auto x102 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::OMEGA>(0);
auto x103 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1) + x52;
auto x104 = V{0.166666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x105 = V{0.333333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x106 = x13*(V{1.33333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x33 + V{0.666666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x38 + V{4.9960036108132e-16}*cell[4] + V{3.99680288865056e-15}*cell[5] + V{1.99840144432528e-15}*cell[6] + V{2.66453525910038e-15}*cell[7] + V{4.9960036108132e-16}*cell[8] + x104*x42 + x104*x49 + x105*x25 + x105*x28 + V{3.33066907387547e-16});
auto x107 = x106 + x69 + x93;
auto x108 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x109 = V{0.083333375}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-0.3333335};
auto x110 = V{0.041666625}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-0.1666665};
auto x111 = V{0.0138889166666667}*x53;
auto x112 = x10*x109 - x11*x110 - x111 + V{0.333333333333333}*x12;
auto x113 = V{0.0185185185185185}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x114 = V{0.00925925925925926}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x115 = V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[7] + x114*x91 + x114*x96 + V{-0.0555555555555555};
auto x116 = -V{0.166666666666667}*cell[4] + V{0.666666666666667}*cell[6] - V{0.166666666666667}*cell[8] + x113*x65 + x113*x66 + x115;
auto x117 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::OMEGA>(0);
auto x118 = x13*(-x32*x75 - x48*x72 - x58*x72 - x59*x73 + x65*x74 + x66*x74 + x76) + x71;
auto x119 = x80 + x83;
auto x120 = x85 + x87;
auto x121 = -V{0.166666666666667}*cell[5] + V{0.833333333333333}*cell[7] + x100 - x91*x95;
auto x122 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1) + x67;
auto x123 = -x122;
auto x124 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x125 = x14*x65 + x14*x66 + x50;
auto x126 = V{0.333333333333333}*cell.template getFieldComponent<descriptors::FORCE>(0)*(-x13*(x125 - x29*(x31 + x62) - x34*x98 - x39*x91 - x39*x96) + V{1}) + x10*x110 - x109*x11 + x111;
auto x127 = V{0.037037037037037}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x128 = -V{0.037037037037037}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x98 + V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[8] + x115;
auto x129 = -x13*(x125 - x34*x59 - x39*x48 - x39*x58 - x60) + V{1};
auto x130 = x129*x129;
auto x131 = -V{1.5}*x130;
auto x132 = x77 + V{3};
auto x133 = V{0.00462962962962963}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x134 = -x103;
auto x135 = -x106 + x45;
auto x136 = x135 + x69;
auto x0 = -V{0.444444444444444}*cell.template getFieldComponent<descriptors::OMEGA>(0)*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x31 + x69) + V{1}) - cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(-V{1.33333333333333}*cell.template getFieldComponent<descriptors::FORCE>(0)*x52 + x10*x9 + x11*x9 + V{1.33333333333333}*x12 + V{0.111111}*x53) - x54*(-V{0.148148148148148}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x59 + V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x65 + V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x66 - V{0.666666666666667}*cell[4] - V{2.66666666666667}*cell[5] - V{1.33333333333333}*cell[6] - V{2.66666666666667}*cell[7] - V{0.666666666666667}*cell[8] - x48*x55 - x55*x58 - x60 + V{-0.888888888888889});
auto x1 = -(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x70*(x71 + x78) - x79*(x84 + V{12}) + x88) + x102*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x107 + x26 - V{4.5}*x103*x103) + V{1}) + x54*(x101 + x89*x91 - x92*(x26 + x90 + x94)));
auto x2 = -(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x108*x78 + x112) + x117*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x107 - V{4.5}*x52*x52) + V{1}) + x54*(V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x98 + x116 - x14*(x94 + x97)));
auto x3 = -(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x120 - x70*(x118 + V{-9}) - x79*(x119 + V{-12})) + x102*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x107 + x17 - V{4.5}*x123*x123) + V{1}) + x54*(x121 + x89*x96 - x92*(x17 + x56 + x62 + x94)));
auto x4 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x124*(x80 + V{-3}) + x126) + x117*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x131 + x27) + V{-1}) - x54*(-x127*x65 + x128 + x55*x66);
auto x5 = -(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x70*(x132 + x71) - x79*(x84 + V{6}) + x88) + x102*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x136 + x17 - V{4.5}*x134*x134) + V{1}) + x54*(x101 - x133*x91));
auto x6 = -cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x108*x132 + x112) - x117*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x135 - V{3}*x68) + V{1}) - x54*(-V{0.037037037037037}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x98 + x116);
auto x7 = -(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x120 - x70*(x118 + V{-3}) - x79*(x119 + V{-6})) + x102*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x136 + x26 - V{4.5}*x122*x122) + V{1}) + x54*(x121 - x133*x96));
auto x8 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x124*(x80 + V{3}) + x126) + x117*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x131 + x64) + V{-1}) - x54*(-x127*x66 + x128 + x55*x65);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { cell.template getFieldComponent<momenta::FixedDensity::RHO>(0), V{1}*x130 + x15 };
}
};

}

}

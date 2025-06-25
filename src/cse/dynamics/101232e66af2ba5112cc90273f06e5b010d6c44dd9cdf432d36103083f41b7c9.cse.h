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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>, forcing::Wagner>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, -1>, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = V{0.666666} - V{0.1666665}*cell.template getFieldComponent<descriptors::OMEGA>(0);
auto x10 = cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x11 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x12 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1);
auto x13 = V{1} / (cell.template getFieldComponent<momenta::FixedDensity::RHO>(0));
auto x14 = V{0.444444444444444}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x15 = x13*(cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + cell[4] + cell[8] + V{1});
auto x16 = V{1} - x15;
auto x17 = x16*x16;
auto x18 = V{1.5}*x17;
auto x19 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1);
auto x20 = V{1.5}*x19;
auto x21 = x20 + V{-1};
auto x22 = x18 + x21;
auto x23 = -x22;
auto x24 = V{0.222222222222222}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x25 = -x16;
auto x26 = x13*(V{3}*cell[0] + V{6}*cell[1] + V{6}*cell[2] + V{6}*cell[3] + V{3}*cell[4] + V{3}*cell[8] + V{3});
auto x27 = V{2} - x26;
auto x28 = x18 + x20;
auto x29 = x27 + x28;
auto x30 = x29 - V{4.5}*x25*x25;
auto x31 = -x30;
auto x32 = V{0.0555555555555555}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x33 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1) + x15 + V{-1};
auto x34 = V{4.5}*(x33*x33);
auto x35 = V{3}*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1);
auto x36 = -x18;
auto x37 = x35 + x36;
auto x38 = -x20 + x26 + x34 + x37 + V{-2};
auto x39 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1) + x16;
auto x40 = -x39;
auto x41 = x29 + x35;
auto x42 = x41 - V{4.5}*x40*x40;
auto x43 = -x42;
auto x44 = V{0.111111111111111}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x45 = V{3}*x19 + V{1};
auto x46 = x37 + x45;
auto x47 = -x35;
auto x48 = x45 + x47;
auto x49 = x36 + x48;
auto x50 = V{1.33226762955019e-15}*cell[1] + V{6.66133814775094e-16}*cell[2] + V{8.88178419700125e-16}*cell[3] + V{1.66533453693773e-16}*cell[4] + V{1.66533453693773e-16}*cell[8] + x44*x46 + x44*x49 + V{1.11022302462516e-16};
auto x51 = x13*(x14*x23 + x24*x31 + x32*x38 + x32*x43 + x50);
auto x52 = x51 + V{-1};
auto x53 = cell.template getFieldComponent<descriptors::OMEGA>(0)*cell.template getFieldComponent<descriptors::SCALAR>(0)*x13;
auto x54 = cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-1};
auto x55 = V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x56 = -x34;
auto x57 = x29 + x47 + x56;
auto x58 = x14*x22;
auto x59 = V{1} - x51;
auto x60 = x59*x59;
auto x61 = V{1.5}*x60;
auto x62 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x63 = V{9}*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1);
auto x64 = V{0.333333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x65 = V{1.33333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x66 = V{2.66666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x67 = V{0.666666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x68 = V{7.99360577730113e-15}*cell[1] + V{3.99680288865056e-15}*cell[2] + V{5.32907051820075e-15}*cell[3] + V{9.99200722162641e-16}*cell[4] + V{9.99200722162641e-16}*cell[8] + x46*x67 + x49*x67 + V{6.66133814775094e-16};
auto x69 = x13*(x23*x66 + x31*x65 + x38*x64 + x43*x64 + x68) + x63;
auto x70 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x71 = V{6}*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1);
auto x72 = V{1}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x73 = V{0.5}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x74 = x13*(-V{4}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x22 - V{2}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x30 + V{1.19904086659517e-14}*cell[1] + V{5.99520433297585e-15}*cell[2] + V{7.99360577730113e-15}*cell[3] + V{1.49880108324396e-15}*cell[4] + V{1.49880108324396e-15}*cell[8] - x42*x73 + x46*x72 + x49*x72 - x57*x73 + V{9.99200722162641e-16});
auto x75 = x71 + x74;
auto x76 = V{0.25}*cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::FORCE>(1)*(V{0.25}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-1});
auto x77 = V{0.02083334375}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-0.083333375};
auto x78 = x10*x77 + x11*x77 - V{0.0138888958333333}*x53;
auto x79 = -x76 + x78;
auto x80 = V{0.00462962962962963}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x81 = V{0.00462962962962963}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x82 = -V{4.5}*x39*x39;
auto x83 = x41 + x82;
auto x84 = -V{3}*x17 + x20;
auto x85 = x27 + x84;
auto x86 = V{0.00925925925925926}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x87 = V{0.0185185185185185}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x85 + V{0.166666666666667}*cell[2] + V{0.0833333333333333}*cell[4] + V{0.0833333333333333}*cell[8] - x46*x86 - x49*x86 + V{0.0277777777777778};
auto x88 = V{0.833333333333333}*cell[1] - V{0.166666666666667}*cell[3] - x81*x83 + x87;
auto x89 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::OMEGA>(0);
auto x90 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1) + x52;
auto x91 = V{0.166666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x92 = V{0.333333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x93 = x13*(V{1.33333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x23 + V{0.666666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x31 + V{3.99680288865056e-15}*cell[1] + V{1.99840144432528e-15}*cell[2] + V{2.66453525910038e-15}*cell[3] + V{4.9960036108132e-16}*cell[4] + V{4.9960036108132e-16}*cell[8] + x38*x91 + x43*x91 + x46*x92 + x49*x92 + V{3.33066907387547e-16});
auto x94 = x20 + x61 - x93 + V{2};
auto x95 = V{0.0185185185185185}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x96 = V{0.00925925925925926}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x97 = V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[3] + x57*x96 + x83*x96 + V{-0.0555555555555555};
auto x98 = V{0.666666666666667}*cell[2] - V{0.166666666666667}*cell[4] - V{0.166666666666667}*cell[8] + x46*x95 + x49*x95 + x97;
auto x99 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x100 = -x22*x66 - x57*x64 + x68;
auto x101 = -x13*(x100 - x64*x83 - x65*x85);
auto x102 = V{0.041666625}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-0.1666665};
auto x103 = V{0.083333375}*cell.template getFieldComponent<descriptors::OMEGA>(0) + V{-0.3333335};
auto x104 = V{0.0138889166666667}*x53;
auto x105 = -x10*x103 + x102*x11 + x104 - V{0.333333333333333}*x12;
auto x106 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::OMEGA>(0);
auto x107 = -x13*(x100 - x30*x65 - x42*x64) + x63;
auto x108 = x71 - x74;
auto x109 = x76 + x78;
auto x110 = -V{0.166666666666667}*cell[1] + V{0.833333333333333}*cell[3] - x57*x81 + x87;
auto x111 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1) + x59;
auto x112 = -x111;
auto x113 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x114 = -x32*x57 + x50 - x58;
auto x115 = -V{0.333333333333333}*cell.template getFieldComponent<descriptors::FORCE>(0)*(-x13*(x114 - x24*x30 - x32*x42) + V{1}) + x10*x102 - x103*x11 + x104;
auto x116 = V{0.037037037037037}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x117 = -V{0.037037037037037}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x85 - V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[8] + x97;
auto x118 = -x13*(x114 - x24*x85 - x32*x83) + V{1};
auto x119 = x118*x118;
auto x120 = -V{1.5}*x119;
auto x121 = V{0.0231481481481481}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x122 = V{0.0277777777777778}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x123 = x26 + V{-4};
auto x124 = x123 + x28;
auto x125 = -x90;
auto x126 = x20 + x93 + V{-4};
auto x127 = x126 + x61;
auto x0 = -V{0.444444444444444}*cell.template getFieldComponent<descriptors::OMEGA>(0)*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x21 + x61) + V{1}) - cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(-V{1.33333333333333}*cell.template getFieldComponent<descriptors::FORCE>(0)*x52 + x10*x9 + x11*x9 + V{1.33333333333333}*x12 + V{0.111111}*x53) - x54*(-V{0.148148148148148}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x30 + V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x46 + V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x49 - V{2.66666666666667}*cell[1] - V{1.33333333333333}*cell[2] - V{2.66666666666667}*cell[3] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[8] - x42*x55 - x55*x57 - x58 + V{-0.888888888888889});
auto x1 = -(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x62*(x69 + V{-3}) - x70*(x75 + V{-6}) + x79) + x54*(-x57*x80 + x88) + x89*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x47 + x94 - V{4.5}*x90*x90) + V{1}));
auto x2 = -(-cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x105 + x99*(x101 + V{3})) + x106*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x94 - V{4.5}*x52*x52) + V{1}) + x54*(-V{0.037037037037037}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x85 + x98));
auto x3 = -(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x109 - x62*(x107 + V{3}) - x70*(x108 + V{6})) + x54*(x110 - x80*x83) + x89*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x35 + x94 - V{4.5}*x112*x112) + V{1}));
auto x4 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x113*(x71 + V{-3}) + x115) + x106*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x120 + x48) + V{-1}) - x54*(-x116*x46 + x117 + x49*x55);
auto x5 = -(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x62*(x69 + V{-9}) - x70*(x75 + V{-12}) + x79) + x54*(x121*x57 - x122*(x124 + x35 + x56) + x88) + x89*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x127 + x35 - V{4.5}*x125*x125) + V{1}));
auto x6 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x105 + x99*(x101 + V{9})) - x106*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x126 - V{3}*x60) + V{1}) - x54*(V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x85 - x44*(x123 + x84) + x98);
auto x7 = -(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x109 - x62*(x107 + V{9}) - x70*(x108 + V{12})) + x54*(x110 + x121*x83 - x122*(x124 + x47 + x82)) + x89*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x127 + x47 - V{4.5}*x111*x111) + V{1}));
auto x8 = cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x113*(x71 + V{3}) + x115) + x106*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x120 + x35 + x45) + V{-1}) - x54*(-x116*x49 + x117 + x46*x55);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { cell.template getFieldComponent<momenta::FixedDensity::RHO>(0), V{1}*x119 + x19 };
}
};

}

}

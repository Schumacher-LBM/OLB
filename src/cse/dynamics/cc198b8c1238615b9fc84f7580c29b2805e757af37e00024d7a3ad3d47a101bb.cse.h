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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x12 = V{3}*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x13 = V{1} / (cell.template getFieldComponent<momenta::FixedDensity::RHO>(0));
auto x14 = V{0.111111111111111}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x15 = x13*(cell[0] + cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + V{1});
auto x16 = V{1} - x15;
auto x17 = -x16;
auto x18 = x17*x17;
auto x19 = V{1.5}*x18;
auto x20 = -x19;
auto x21 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x22 = V{3}*x21 + V{1};
auto x23 = x20 + x22;
auto x24 = x12 + x23;
auto x25 = -x12;
auto x26 = x23 + x25;
auto x27 = V{0.444444444444444}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x28 = V{1.5}*x21;
auto x29 = x28 + V{-1};
auto x30 = x19 + x29;
auto x31 = -x30;
auto x32 = V{0.222222222222222}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x33 = V{3}*x18;
auto x34 = x13*(V{3}*cell[0] + V{3}*cell[4] + V{6}*cell[5] + V{6}*cell[6] + V{6}*cell[7] + V{3}*cell[8] + V{3});
auto x35 = -x28 + V{-2};
auto x36 = x34 + x35;
auto x37 = x33 + x36;
auto x38 = V{0.0555555555555555}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x39 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1) + x15 + V{-1};
auto x40 = V{4.5}*(x39*x39);
auto x41 = x12 + x20 + x36 + x40;
auto x42 = -cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1) - x16;
auto x43 = x28 + V{2};
auto x44 = -x34 + x43;
auto x45 = x19 + x44;
auto x46 = x12 + x45 - V{4.5}*x42*x42;
auto x47 = -x46;
auto x48 = V{1.66533453693773e-16}*cell[4] + V{1.33226762955019e-15}*cell[5] + V{6.66133814775094e-16}*cell[6] + V{8.88178419700125e-16}*cell[7] + V{1.66533453693773e-16}*cell[8] + V{1.11022302462516e-16};
auto x49 = x13*(x14*x24 + x14*x26 + x27*x31 + x32*x37 + x38*x41 + x38*x47 + x48);
auto x50 = V{1} - x49;
auto x51 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1) + x50;
auto x52 = -x51;
auto x53 = x49 + V{-1};
auto x54 = x53*x53;
auto x55 = V{1.5}*x54;
auto x56 = V{0.166666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x57 = V{0.333333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x58 = x13*(V{1.33333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x31 + V{0.666666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x37 + V{4.9960036108132e-16}*cell[4] + V{3.99680288865056e-15}*cell[5] + V{1.99840144432528e-15}*cell[6] + V{2.66453525910038e-15}*cell[7] + V{4.9960036108132e-16}*cell[8] + x24*x57 + x26*x57 + x41*x56 + x47*x56 + V{3.33066907387547e-16});
auto x59 = x43 - x58;
auto x60 = x55 + x59;
auto x61 = x12 + x60 - V{4.5}*x52*x52;
auto x62 = cell.template getFieldComponent<momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1) + x53;
auto x63 = V{4.5}*(x62*x62);
auto x64 = x25 + x60 - x63;
auto x65 = V{3}*x54;
auto x66 = x59 - x65;
auto x67 = V{1.15648231731787e-17}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x68 = x25 - x40 + x45;
auto x69 = -x33 + x44;
auto x70 = x22 - V{1.5}*x16*x16;
auto x71 = x12 + x70;
auto x72 = x25 + x70;
auto x73 = x50*x50;
auto x74 = x22 - V{1.5}*x73;
auto x75 = x12 + x74;
auto x76 = x25 + x74;
auto x77 = V{0.0277777777777778}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x78 = x28 + x55 + x58 + V{-4};
auto x79 = V{1.01192202765314e-18}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x80 = V{1.87928376564154e-18}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x81 = V{0.0231481481481481}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x82 = -x61;
auto x83 = V{0.00462962962962963}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x84 = -x55;
auto x85 = x35 + x58;
auto x86 = x12 + x63 + x84 + x85;
auto x87 = V{1.73472347597681e-18}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x88 = V{0.00925925925925926}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x89 = x22 + x84;
auto x90 = x12 + x89;
auto x91 = x25 + x89;
auto x92 = x65 + x85;
auto x93 = V{4.62592926927149e-18}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x37 + V{0.0185185185185185}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x92 - V{0.0833333333333334}*cell[4] - V{0.166666666666667}*cell[6] - V{0.0833333333333334}*cell[8] + x24*x87 + x26*x87 + x88*x90 + x88*x91 + V{-0.0555555555555556};
auto x94 = -x10*(-V{0.833333333333333}*cell[5] + V{0.166666666666667}*cell[7] + x41*x80 + x47*x79 + x81*x82 - x83*x86 + x93) + V{0.0277777777777778};
auto x95 = V{0.333333333333333}*cell[7];
auto x96 = V{0.333333333333334}*cell[5];
auto x97 = V{3.85185988877447e-34}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x98 = V{0.00925925925925926}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x99 = V{1.44560289664734e-18}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x100 = V{0.0185185185185185}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x101 = -x10*(V{1.15648231731787e-17}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x37 + V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x92 + V{0.166666666666667}*cell[4] - V{0.666666666666667}*cell[6] + V{0.166666666666667}*cell[8] - x100*x90 - x100*x91 + x24*x97 + x26*x97 + x41*x99 + x47*x99 + x82*x98 + x86*x98 - x95 - x96 + V{-0.0555555555555556}) + V{0.111111111111111};
auto x102 = -x62;
auto x103 = -x10*(V{0.166666666666666}*cell[5] - V{0.833333333333333}*cell[7] + x41*x79 + x47*x80 + x81*x86 - x82*x83 + x93) + V{0.0277777777777778};
auto x104 = -x13*(x14*x71 + x14*x72 - x27*x30 - x32*x69 - x38*x46 - x38*x68 + x48) + V{1};
auto x105 = x104*x104;
auto x106 = -V{1.5}*x105 + x22;
auto x107 = x10*(V{0.037037037037037}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x66 + V{2.31296463463574e-18}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x69 + V{3.46944695195361e-18}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x71 + V{3.46944695195361e-18}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x72 + V{0.037037037037037}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x75 + V{0.037037037037037}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x76 - V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[6] - V{0.333333333333333}*cell[8] - x46*x99 - x61*x98 - x64*x98 - x68*x99 - x95 - x96 + V{-0.0555555555555556}) + V{-0.111111111111111};
auto x0 = -x10*(-V{0.148148148148148}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x66 - V{3.70074341541719e-17}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x69 + V{1.38777878078145e-17}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x71 + V{1.38777878078145e-17}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x72 + V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x75 + V{0.0740740740740741}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x76 - V{0.666666666666667}*cell[4] - V{2.66666666666667}*cell[5] - V{1.33333333333333}*cell[6] - V{2.66666666666667}*cell[7] - V{0.666666666666667}*cell[8] - x11*x61 - x11*x64 - x46*x67 - x67*x68 + V{-0.444444444444444}) - x27*(x29 + x55) + V{-0.444444444444444};
auto x1 = -(x77*(x25 + x78 - V{4.5}*x51*x51) + x94);
auto x2 = -x101 - x14*(-V{4.5}*x73 + x78);
auto x3 = -(x103 + x77*(x12 + x78 - V{4.5}*x102*x102));
auto x4 = x107 + x14*(x106 + x25);
auto x5 = -x61*x77 - x94;
auto x6 = -x101 - x14*x66;
auto x7 = -x103 - x64*x77;
auto x8 = x107 + x14*(x106 + x12);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { cell.template getFieldComponent<momenta::FixedDensity::RHO>(0), V{1}*x105 + x21 };
}
};

}

}

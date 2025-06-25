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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1};
auto x12 = V{1} / (x11);
auto x13 = cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + cell[4] + cell[8] + V{1};
auto x14 = x12*x13;
auto x15 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x16 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x17 = V{3}*x16;
auto x18 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x19 = V{1.5}*x18;
auto x20 = -x19;
auto x21 = x17 + x20 + V{1};
auto x22 = x15 + x21;
auto x23 = V{0.0555555555555555}*x14;
auto x24 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x25 = -x24;
auto x26 = V{4.5}*(x25*x25);
auto x27 = -x15;
auto x28 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x29 = V{1.5}*x16;
auto x30 = x19 + V{-1};
auto x31 = x29 + x30;
auto x32 = x28 + x31;
auto x33 = x27 + x32;
auto x34 = -x26 + x33;
auto x35 = V{0.111111111111111}*x12;
auto x36 = x13*x35;
auto x37 = x15 - x17 + x30;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x39 = V{4.5}*(x38*x38);
auto x40 = x15 + x32 - x39;
auto x41 = -x28;
auto x42 = V{1} - x29;
auto x43 = V{3}*x18 + x42;
auto x44 = x41 + x43;
auto x45 = V{1.33226762955019e-15}*cell[1] + V{6.66133814775094e-16}*cell[2] + V{8.88178419700125e-16}*cell[3] + V{1.66533453693773e-16}*cell[4] + V{1.66533453693773e-16}*cell[8] + V{1.11022302462516e-16};
auto x46 = V{0.444444444444444}*x12*x13*x31 - V{0.222222222222222}*x14*x44 - x22*x36 + x23*x40 + x36*x37 + x45;
auto x47 = x23*x34 + x46;
auto x48 = x22*x47;
auto x49 = x44*x47;
auto x50 = V{0.0277777777777778}*x12;
auto x51 = -V{4.5}*x24*x24;
auto x52 = x33 + x51;
auto x53 = x23*x52 + x46;
auto x54 = V{0.00462962962962963}*x12;
auto x55 = x40*x47;
auto x56 = V{0.0833333333333334}*cell[4];
auto x57 = V{0.0833333333333334}*cell[8];
auto x58 = V{0.166666666666667}*cell[2];
auto x59 = -V{1.73472347597681e-18}*x12*x13*x37 - V{0.00925925925925926}*x12*x37*x47 + V{0.00925925925925926}*x12*x48 + V{0.0185185185185185}*x12*x49 + V{1.73472347597681e-18}*x14*x22 + V{4.62592926927149e-18}*x14*x44 + x56 + x57 + x58 + V{0.0555555555555556};
auto x60 = x10*(-V{0.833333333333333}*cell[1] + V{0.166666666666667}*cell[3] + V{1.01192202765314e-18}*x12*x13*x34 + V{1.87928376564154e-18}*x12*x13*x40 + V{0.0231481481481481}*x12*x34*x47 - x54*x55 - x59) + V{-0.0277777777777778};
auto x61 = -1/x11;
auto x62 = x13*x61;
auto x63 = V{3.85185988877447e-34}*x62;
auto x64 = x21 + x27;
auto x65 = V{0.111111111111111}*x62;
auto x66 = V{0.0555555555555555}*x62;
auto x67 = x15 + x41;
auto x68 = x20 + x42;
auto x69 = x26 + x67 + x68;
auto x70 = -x40;
auto x71 = x22*x65 - V{0.444444444444444}*x31*x62 + V{0.222222222222222}*x44*x62 + x45 + x64*x65 + x66*x69 + x66*x70;
auto x72 = x44*x71;
auto x73 = V{0.0185185185185185}*x61;
auto x74 = x22*x71;
auto x75 = x64*x71;
auto x76 = V{0.00925925925925926}*x61;
auto x77 = x69*x71;
auto x78 = x70*x71;
auto x79 = V{1.44560289664734e-18}*x62;
auto x80 = -V{0.333333333333334}*cell[1] - V{0.333333333333333}*cell[3] + V{-0.0555555555555556};
auto x81 = x69*x79 + x70*x79 + x76*x77 + x76*x78 + x80;
auto x82 = -x10*(-V{0.666666666666667}*cell[2] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[8] + x22*x63 + V{1.15648231731787e-17}*x44*x62 + V{0.0740740740740741}*x61*x72 + x63*x64 - x73*x74 - x73*x75 + x81) + V{0.111111111111111};
auto x83 = V{0.833333333333333}*cell[3];
auto x84 = x34*x47;
auto x85 = V{0.00925925925925926}*x12;
auto x86 = V{1.44560289664734e-18}*x14;
auto x87 = V{3.46944695195361e-18}*x14;
auto x88 = V{0.037037037037037}*x12;
auto x89 = V{0.333333333333333}*cell[2] - V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[8];
auto x90 = V{1.73472347597681e-18}*x62;
auto x91 = V{0.00925925925925926}*x61;
auto x92 = V{3.46944695195361e-18}*x62;
auto x93 = V{0.037037037037037}*x61;
auto x0 = -x10*(-V{2.66666666666667}*cell[1] - V{1.33333333333333}*cell[2] - V{2.66666666666667}*cell[3] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[8] + V{1.15648231731787e-17}*x12*x13*x34 + V{1.38777878078145e-17}*x12*x13*x37 + V{1.15648231731787e-17}*x12*x13*x40 + V{0.0740740740740741}*x12*x34*x47 + V{0.0740740740740741}*x12*x37*x47 + V{0.0740740740740741}*x12*x40*x47 - V{0.0740740740740741}*x12*x48 - V{0.148148148148148}*x12*x49 - V{1.38777878078145e-17}*x14*x22 - V{3.70074341541719e-17}*x14*x44 + V{-0.444444444444444}) + V{0.444444444444444}*x12*x31*x47 + V{-0.444444444444444};
auto x1 = x50*x52*x53 + x60;
auto x2 = -x35*x49 - x82;
auto x3 = x10*(V{0.166666666666666}*cell[1] + V{1.87928376564154e-18}*x12*x13*x34 + V{1.01192202765314e-18}*x12*x13*x40 + V{0.0231481481481481}*x12*x40*x47 - x54*x84 - x59 - x83) + x40*x50*x53 + V{-0.0277777777777778};
auto x4 = x10*(V{0.037037037037037}*x12*x49 + V{2.31296463463574e-18}*x14*x44 - x22*x87 + x34*x86 + x37*x47*x88 + x37*x87 + x40*x86 - x48*x88 + x55*x85 + x80 + x84*x85 + x89) + x35*x37*x53 + V{-0.111111111111111};
auto x5 = x50*x53*(x31 + x51 + x67) + x60;
auto x6 = -x35*x47*(x28 + x43) - x82;
auto x7 = x10*(V{0.166666666666666}*cell[1] + x22*x90 + V{4.62592926927149e-18}*x44*x62 - x56 - x57 - x58 + V{0.0185185185185185}*x61*x72 - V{0.00462962962962963}*x61*x77 + V{0.0231481481481481}*x61*x78 + V{1.87928376564154e-18}*x62*x69 + V{1.01192202765314e-18}*x62*x70 + x64*x90 + x74*x91 + x75*x91 - x83 + V{-0.0555555555555556}) - x47*x50*(x15 + x28 + x39 + x68) + V{-0.0277777777777778};
auto x8 = x10*(x22*x92 - V{2.31296463463574e-18}*x44*x62 - V{0.037037037037037}*x61*x72 + x64*x92 + x74*x93 + x75*x93 + x81 + x89) - x35*x48 + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -V{1}*x12*x53, x16 + x18 };
}
};

}

}

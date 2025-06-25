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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{-1};
auto x12 = V{1} / (x11);
auto x13 = cell[0] + cell[2] + V{2}*cell[3] + V{2}*cell[4] + V{2}*cell[5] + cell[6] + V{1};
auto x14 = x12*x13;
auto x15 = V{1.38777878078145e-17}*x14;
auto x16 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x17 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x18 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x19 = V{1.5}*x18;
auto x20 = -x19;
auto x21 = V{3}*x17 + x20 + V{1};
auto x22 = x16 + x21;
auto x23 = -x16;
auto x24 = x21 + x23;
auto x25 = V{0.0740740740740741}*x12;
auto x26 = V{1.5}*x17;
auto x27 = x26 + V{-1};
auto x28 = x19 + x27;
auto x29 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x30 = V{3}*x18;
auto x31 = x27 + x29 - x30;
auto x32 = V{0.0555555555555555}*x14;
auto x33 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x34 = V{4.5}*(x33*x33);
auto x35 = x28 + x29;
auto x36 = x16 - x34 + x35;
auto x37 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x38 = V{4.5}*(x37*x37);
auto x39 = -x38;
auto x40 = x23 + x35 + x39;
auto x41 = V{0.111111111111111}*x12;
auto x42 = x13*x41;
auto x43 = V{1.66533453693773e-16}*cell[2] + V{1.33226762955019e-15}*cell[3] + V{6.66133814775094e-16}*cell[4] + V{8.88178419700125e-16}*cell[5] + V{1.66533453693773e-16}*cell[6] + V{1.11022302462516e-16};
auto x44 = V{0.444444444444444}*x12*x13*x28 + V{0.222222222222222}*x14*x31 - x22*x42 - x24*x42 + x32*x36 + x32*x40 + x43;
auto x45 = x22*x44;
auto x46 = x24*x44;
auto x47 = V{0.0277777777777778}*x12;
auto x48 = -x29;
auto x49 = x16 + x48;
auto x50 = V{0.00462962962962963}*x12;
auto x51 = x36*x44;
auto x52 = V{0.0833333333333334}*cell[2];
auto x53 = V{0.0833333333333334}*cell[6];
auto x54 = V{0.166666666666667}*cell[4];
auto x55 = V{1.73472347597681e-18}*x14;
auto x56 = V{0.00925925925925926}*x12;
auto x57 = -V{4.62592926927149e-18}*x12*x13*x31 - V{0.0185185185185185}*x12*x31*x44 + x22*x55 + x24*x55 + x45*x56 + x46*x56 + x52 + x53 + x54 + V{0.0555555555555556};
auto x58 = x10*(V{0.166666666666666}*cell[3] - V{0.833333333333333}*cell[5] + V{1.87928376564154e-18}*x12*x13*x36 + V{1.01192202765314e-18}*x12*x13*x40 + V{0.0231481481481481}*x12*x40*x44 - x50*x51 - x57) + V{-0.0277777777777778};
auto x59 = -1/x11;
auto x60 = x13*x59;
auto x61 = V{3.46944695195361e-18}*x60;
auto x62 = V{0.037037037037037}*x59;
auto x63 = V{0.111111111111111}*x60;
auto x64 = V{1} - x26;
auto x65 = x30 + x64;
auto x66 = x48 + x65;
auto x67 = V{0.0555555555555555}*x60;
auto x68 = x20 + x64;
auto x69 = x38 + x49 + x68;
auto x70 = -x36;
auto x71 = x22*x63 + x24*x63 - V{0.444444444444444}*x28*x60 + x43 + V{0.222222222222222}*x60*x66 + x67*x69 + x67*x70;
auto x72 = x22*x71;
auto x73 = x24*x71;
auto x74 = x66*x71;
auto x75 = V{0.00925925925925926}*x59;
auto x76 = x69*x71;
auto x77 = x70*x71;
auto x78 = V{1.44560289664734e-18}*x60;
auto x79 = -V{0.333333333333334}*cell[3] - V{0.333333333333333}*cell[5] + V{-0.0555555555555556};
auto x80 = x69*x78 + x70*x78 + x75*x76 + x75*x77 + x79;
auto x81 = -x10*(-V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[4] - V{0.333333333333333}*cell[6] + x22*x61 + x24*x61 - V{0.037037037037037}*x59*x74 - V{2.31296463463574e-18}*x60*x66 + x62*x72 + x62*x73 + x80) + V{0.111111111111111};
auto x82 = V{0.833333333333333}*cell[3];
auto x83 = x40*x44;
auto x84 = V{0.0185185185185185}*x12;
auto x85 = V{0.00925925925925926}*x12;
auto x86 = V{1.44560289664734e-18}*x14;
auto x87 = x31*x44;
auto x88 = V{3.85185988877447e-34}*x14;
auto x89 = V{0.166666666666667}*cell[2] - V{0.666666666666667}*cell[4] + V{0.166666666666667}*cell[6];
auto x90 = V{1.73472347597681e-18}*x60;
auto x91 = V{0.00925925925925926}*x59;
auto x92 = V{3.85185988877447e-34}*x60;
auto x93 = V{0.0185185185185185}*x59;
auto x0 = -x10*(-V{0.666666666666667}*cell[2] - V{2.66666666666667}*cell[3] - V{1.33333333333333}*cell[4] - V{2.66666666666667}*cell[5] - V{0.666666666666667}*cell[6] + V{3.70074341541719e-17}*x12*x13*x31 + V{1.15648231731787e-17}*x12*x13*x36 + V{1.15648231731787e-17}*x12*x13*x40 + V{0.148148148148148}*x12*x31*x44 + V{0.0740740740740741}*x12*x36*x44 + V{0.0740740740740741}*x12*x40*x44 - x15*x22 - x15*x24 - x25*x45 - x25*x46 + V{-0.444444444444444}) + V{0.444444444444444}*x12*x28*x44 + V{-0.444444444444444};
auto x1 = x44*x47*(x28 + x39 + x49) + x58;
auto x2 = -x41*x46 - x81;
auto x3 = x10*(V{0.166666666666667}*cell[5] + V{1.01192202765314e-18}*x12*x13*x36 + V{1.87928376564154e-18}*x12*x13*x40 + V{0.0231481481481481}*x12*x36*x44 - x50*x83 - x57 - x82) + x47*x51 + V{-0.0277777777777778};
auto x4 = x10*(V{0.0740740740740741}*x12*x87 + V{1.15648231731787e-17}*x14*x31 - x22*x88 - x24*x88 + x36*x86 + x40*x86 + x45*x84 + x46*x84 + x51*x85 + x79 + x83*x85 + x89) + x41*x87 + V{-0.111111111111111};
auto x5 = x47*x83 + x58;
auto x6 = -x41*x45 - x81;
auto x7 = x10*(V{0.166666666666667}*cell[5] + x22*x90 + x24*x90 - x52 - x53 - x54 + V{0.0185185185185185}*x59*x74 - V{0.00462962962962963}*x59*x76 + V{0.0231481481481481}*x59*x77 + V{4.62592926927149e-18}*x60*x66 + V{1.87928376564154e-18}*x60*x69 + V{1.01192202765314e-18}*x60*x70 + x72*x91 + x73*x91 - x82 + V{-0.0555555555555556}) - x44*x47*(x16 + x29 + x34 + x68) + V{-0.0277777777777778};
auto x8 = x10*(x22*x92 + x24*x92 + V{0.0740740740740741}*x59*x74 + V{1.15648231731787e-17}*x60*x66 - x72*x93 - x73*x93 + x80 + x89) - x41*x44*(x29 + x65) + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -V{1}*x12*x44, x17 + x18 };
}
};

}

}

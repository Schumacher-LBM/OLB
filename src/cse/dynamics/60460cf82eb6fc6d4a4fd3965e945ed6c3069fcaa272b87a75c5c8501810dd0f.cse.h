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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x12 = cell[0] + cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + V{1};
auto x13 = x11*x12;
auto x14 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x15 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x16 = V{3}*x15;
auto x17 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x18 = V{1.5}*x17;
auto x19 = x18 + V{-1};
auto x20 = x14 - x16 + x19;
auto x21 = V{0.0740740740740741}*x11;
auto x22 = V{0.111111111111111}*x11;
auto x23 = x12*x22;
auto x24 = V{0.444444444444444}*x11;
auto x25 = x12*x24;
auto x26 = V{1.5}*x15;
auto x27 = x19 + x26;
auto x28 = V{0.0555555555555555}*x13;
auto x29 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x30 = -x29;
auto x31 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x32 = V{4.5}*(x31*x31);
auto x33 = x14 + x27;
auto x34 = x30 - x32 + x33;
auto x35 = V{1} - x18;
auto x36 = x16 + x35;
auto x37 = x14 + x36;
auto x38 = -x26;
auto x39 = x29 + x38;
auto x40 = V{3}*x17 + V{1};
auto x41 = x39 + x40;
auto x42 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x43 = V{4.5}*(x42*x42);
auto x44 = x35 + x39;
auto x45 = x14 + x43 + x44;
auto x46 = V{1.66533453693773e-16}*cell[4] + V{1.33226762955019e-15}*cell[5] + V{6.66133814775094e-16}*cell[6] + V{8.88178419700125e-16}*cell[7] + V{1.66533453693773e-16}*cell[8] + V{0.222222222222222}*x13*x41 + x23*x37 + x28*x45 + V{1.11022302462516e-16};
auto x47 = -x20*x23 - x25*x27 - x28*x34 + x46;
auto x48 = x20*x47;
auto x49 = x34*x47;
auto x50 = V{0.0277777777777778}*x11;
auto x51 = -x14;
auto x52 = -x31;
auto x53 = V{1.01192202765314e-18}*x13;
auto x54 = x32 + x44 + x51;
auto x55 = V{1.87928376564154e-18}*x13;
auto x56 = V{0.0231481481481481}*x11;
auto x57 = x36 + x51;
auto x58 = x23*x57 - x25*x27 + x28*x54 + x46;
auto x59 = x54*x58;
auto x60 = V{0.00462962962962963}*x11;
auto x61 = x45*x58;
auto x62 = V{1.73472347597681e-18}*x13;
auto x63 = V{0.00925925925925926}*x11;
auto x64 = x37*x58;
auto x65 = x57*x58;
auto x66 = V{0.0185185185185185}*x11;
auto x67 = x41*x58;
auto x68 = -V{0.0833333333333334}*cell[4] - V{0.166666666666667}*cell[6] - V{0.0833333333333334}*cell[8] + V{4.62592926927149e-18}*x13*x41 + x37*x62 + V{-0.0555555555555556};
auto x69 = x57*x62 + x63*x64 + x63*x65 + x66*x67 + x68;
auto x70 = -x10*(-V{0.833333333333333}*cell[5] + V{0.166666666666667}*cell[7] + x45*x55 + x53*x54 + x56*x59 - x60*x61 + x69) + V{0.0277777777777778};
auto x71 = V{3.85185988877447e-34}*x13;
auto x72 = V{0.0185185185185185}*x11;
auto x73 = V{0.00925925925925926}*x11;
auto x74 = x45*x47;
auto x75 = x41*x47;
auto x76 = x37*x47;
auto x77 = x49*x73;
auto x78 = V{1.44560289664734e-18}*x13;
auto x79 = x34*x78;
auto x80 = V{0.333333333333333}*cell[7];
auto x81 = V{0.333333333333334}*cell[5];
auto x82 = x45*x78 - x80 - x81 + V{-0.0555555555555556};
auto x83 = x10*(V{0.166666666666667}*cell[4] - V{0.666666666666667}*cell[6] + V{0.166666666666667}*cell[8] + V{0.0740740740740741}*x11*x75 + V{1.15648231731787e-17}*x13*x41 - x20*x71 + x37*x71 + x48*x72 - x72*x76 + x73*x74 - x77 - x79 + x82) + V{-0.111111111111111};
auto x84 = V{0.166666666666666}*cell[5] - V{0.833333333333333}*cell[7] + x45*x53;
auto x85 = V{0.333333333333333}*cell[4];
auto x86 = V{0.333333333333333}*cell[8];
auto x87 = V{3.46944695195361e-18}*x13;
auto x88 = V{0.037037037037037}*x11;
auto x89 = V{0.037037037037037}*x11;
auto x90 = V{2.31296463463574e-18}*x13*x41;
auto x0 = -x10*(-V{0.666666666666667}*cell[4] - V{2.66666666666667}*cell[5] - V{1.33333333333333}*cell[6] - V{2.66666666666667}*cell[7] - V{0.666666666666667}*cell[8] + V{1.38777878078145e-17}*x11*x12*x37 + V{3.70074341541719e-17}*x11*x12*x41 + V{1.15648231731787e-17}*x11*x12*x45 + V{0.0740740740740741}*x11*x37*x47 + V{0.148148148148148}*x11*x41*x47 + V{0.0740740740740741}*x11*x45*x47 - V{1.38777878078145e-17}*x13*x20 - V{1.15648231731787e-17}*x13*x34 - x21*x48 - x21*x49 + V{-0.444444444444444}) - x24*x27*x47 + V{-0.444444444444444};
auto x1 = -(x47*x50*(x27 + x29 + x51 - V{4.5}*x52*x52) + x70);
auto x2 = x22*x47*(x30 + x38 + x40) + x83;
auto x3 = x10*(x54*x55 + x56*x61 - x59*x60 + x69 + x84) - x47*x50*(x29 + x33 - x43) + V{-0.0277777777777778};
auto x4 = x10*(V{0.333333333333333}*cell[6] + x37*x87 + x54*x78 + x57*x87 + x59*x73 + x61*x73 + x64*x88 + x65*x88 - x67*x89 + x82 - x85 - x86 - x90) - x22*x48 + V{-0.111111111111111};
auto x5 = -x49*x50 - x70;
auto x6 = x22*x75 + x83;
auto x7 = x10*(-x20*x62 - x34*x55 - x48*x63 + x49*x60 + x56*x74 + x63*x76 + x66*x75 + x68 + x84) + x50*x74 + V{-0.0277777777777778};
auto x8 = x10*(V{0.333333333333333}*cell[6] + V{3.46944695195361e-18}*x11*x12*x37 + V{1.44560289664734e-18}*x11*x12*x45 + V{0.037037037037037}*x11*x37*x47 + V{0.00925925925925926}*x11*x45*x47 - x20*x87 - x48*x88 - x75*x89 - x77 - x79 - x80 - x81 - x85 - x86 - x90 + V{-0.0555555555555556}) + x22*x76 + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { V{1}*x11*x47, x15 + x17 };
}
};

}

}

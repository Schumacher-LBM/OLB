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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x12 = cell[0] + V{2}*cell[1] + cell[2] + cell[6] + V{2}*cell[7] + V{2}*cell[8] + V{1};
auto x13 = x11*x12;
auto x14 = V{1.38777878078145e-17}*x13;
auto x15 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x16 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x17 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x18 = V{1.5}*x17;
auto x19 = V{1} - x18;
auto x20 = V{3}*x16 + x19;
auto x21 = x15 + x20;
auto x22 = -x15;
auto x23 = x20 + x22;
auto x24 = V{1.04083408558608e-17}*x13;
auto x25 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x26 = V{4.5}*(x25*x25);
auto x27 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x28 = V{1.5}*x16;
auto x29 = x27 - x28;
auto x30 = x19 + x29;
auto x31 = x15 + x26 + x30;
auto x32 = V{0.0740740740740741}*x11;
auto x33 = V{0.0555555555555555}*x13;
auto x34 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x35 = -x34;
auto x36 = V{4.5}*(x35*x35);
auto x37 = x28 + V{-1};
auto x38 = x18 + x37;
auto x39 = x15 + x38;
auto x40 = -x27 + x39;
auto x41 = -x36 + x40;
auto x42 = V{0.444444444444444}*x11;
auto x43 = x12*x42;
auto x44 = V{0.111111111111111}*x11;
auto x45 = x12*x44;
auto x46 = V{3}*x17;
auto x47 = x29 + x46 + V{1};
auto x48 = V{8.88178419700125e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{1.66533453693773e-16}*cell[6] + V{8.88178419700125e-16}*cell[7] + V{2.22044604925031e-16}*cell[8] + V{0.222222222222222}*x13*x47 + x21*x45 + x23*x45 + x31*x33 + V{1.11022302462516e-16};
auto x49 = -x38*x43 + x48;
auto x50 = -x33*x41 + x49;
auto x51 = x21*x50;
auto x52 = x23*x50;
auto x53 = x31*x50;
auto x54 = x47*x50;
auto x55 = x41*x50;
auto x56 = V{0.0277777777777778}*x11;
auto x57 = V{1.73472347597681e-18}*x13;
auto x58 = V{0.0231481481481481}*x11;
auto x59 = x22 + x30 + x36;
auto x60 = x33*x59 - x38*x43 + x48;
auto x61 = x59*x60;
auto x62 = V{8.67361737988404e-19}*x13;
auto x63 = V{0.00462962962962963}*x11;
auto x64 = x31*x60;
auto x65 = V{0.00925925925925926}*x11;
auto x66 = x21*x60;
auto x67 = x23*x60;
auto x68 = V{0.0185185185185185}*x11;
auto x69 = x47*x60;
auto x70 = -V{0.0833333333333334}*cell[2] - V{0.0833333333333334}*cell[6] - V{0.166666666666667}*cell[8] + V{4.62592926927149e-18}*x13*x47 + x21*x57 + x23*x57 + V{-0.0555555555555556};
auto x71 = x65*x66 + x65*x67 + x68*x69 + x70;
auto x72 = -x10*(-V{0.833333333333333}*cell[1] + V{0.166666666666667}*cell[7] + x31*x57 + x58*x61 + x59*x62 - x63*x64 + x71) + V{0.0277777777777778};
auto x73 = -V{4.5}*x34*x34;
auto x74 = -x33*(x40 + x73) + x49;
auto x75 = V{0.333333333333333}*cell[7];
auto x76 = V{0.00925925925925926}*x11;
auto x77 = x55*x76;
auto x78 = x10*(-V{0.333333333333334}*cell[1] - V{0.333333333333333}*cell[2] - V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[8] + V{3.46944695195361e-18}*x11*x12*x21 + V{3.46944695195361e-18}*x11*x12*x23 + V{1.73472347597681e-18}*x11*x12*x31 + V{0.037037037037037}*x11*x21*x50 + V{0.037037037037037}*x11*x23*x50 + V{0.00925925925925926}*x11*x31*x50 - V{0.037037037037037}*x11*x54 - V{1.73472347597681e-18}*x13*x41 - V{2.31296463463574e-18}*x13*x47 - x75 - x77 + V{-0.0555555555555556}) + V{-0.111111111111111};
auto x79 = V{0.166666666666667}*cell[1] - V{0.833333333333333}*cell[7] + x31*x62;
auto x80 = V{8.67361737988404e-19}*x13;
auto x81 = V{0.0740740740740741}*x11;
auto x82 = V{0.0185185185185185}*x11;
auto x83 = V{3.85185988877447e-34}*x13;
auto x84 = -V{0.333333333333333}*cell[1] + V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] - V{0.666666666666667}*cell[8] + V{1.15648231731787e-17}*x13*x47 + x21*x83 + x23*x83 + x31*x80 - x75 + V{-0.0555555555555556};
auto x0 = -x10*(-V{2.66666666666667}*cell[1] - V{0.666666666666667}*cell[2] - V{0.666666666666667}*cell[6] - V{2.66666666666667}*cell[7] - V{1.33333333333333}*cell[8] + V{0.148148148148148}*x11*x54 + V{3.70074341541719e-17}*x13*x47 + x14*x21 + x14*x23 + x24*x31 - x24*x41 + x32*x51 + x32*x52 + x32*x53 - x32*x55 + V{-0.444444444444445}) - x38*x42*x50 + V{-0.444444444444444};
auto x1 = -x55*x56 - x72;
auto x2 = x23*x44*x74 + x78;
auto x3 = x10*(x57*x59 + x58*x64 - x61*x63 + x71 + x79) - x50*x56*(-x26 + x27 + x39) + V{-0.0277777777777778};
auto x4 = x10*(x59*x80 + x61*x76 + x64*x76 - x66*x82 - x67*x82 + x69*x81 + x84) - x44*x50*(x27 + x37 - x46) + V{-0.111111111111111};
auto x5 = -x50*x56*(x22 + x27 + x38 + x73) - x72;
auto x6 = x21*x44*x74 + x78;
auto x7 = x10*(-x41*x57 + x51*x65 + x52*x65 + x53*x58 + x54*x68 + x55*x63 + x70 + x79) + x31*x56*x74 + V{-0.0277777777777778};
auto x8 = x10*(-x41*x80 - x51*x82 - x52*x82 + x53*x76 + x54*x81 - x77 + x84) + x44*x47*x74 + V{-0.111111111111111};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { V{1}*x11*x74, x16 + x17 };
}
};

}

}

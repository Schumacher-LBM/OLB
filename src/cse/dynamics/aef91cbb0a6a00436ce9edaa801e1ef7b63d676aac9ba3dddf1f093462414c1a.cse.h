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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<1, 1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0);
auto x22 = V{1.5}*x21;
auto x23 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2);
auto x24 = V{1.5}*x23;
auto x25 = x22 + x24;
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedDensity::RHO>(0));
auto x27 = x26*(cell[0] + cell[10] + V{2}*cell[11] + cell[12] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[18] + cell[1] + cell[3] + V{2}*cell[5] + cell[6] + cell[7] + V{1});
auto x28 = V{1} - x27;
auto x29 = -x28;
auto x30 = x29*x29;
auto x31 = V{1.5}*x30;
auto x32 = x31 + V{-1};
auto x33 = x25 + x32;
auto x34 = V{0.0555555555555556}*x19;
auto x35 = V{3}*cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0);
auto x36 = V{3}*x21;
auto x37 = x28*x28;
auto x38 = x26*(V{3}*cell[0] + V{3}*cell[10] + V{6}*cell[11] + V{3}*cell[12] + V{6}*cell[13] + V{3}*cell[15] + V{3}*cell[16] + V{6}*cell[17] + V{6}*cell[18] + V{3}*cell[1] + V{3}*cell[3] + V{6}*cell[5] + V{3}*cell[6] + V{3}*cell[7] + V{3});
auto x39 = x25 + x31 + x38 + V{-4};
auto x40 = V{3}*cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2);
auto x41 = V{3}*x23;
auto x42 = V{0.0277777777777778}*x19;
auto x43 = x27 + V{-1};
auto x44 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0) + x43;
auto x45 = -x44;
auto x46 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0) + x28;
auto x47 = -x46;
auto x48 = x25 - x38 + V{2};
auto x49 = x31 + x48;
auto x50 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2);
auto x51 = V{4.5}*(x50*x50);
auto x52 = x33 + x35;
auto x53 = -x40;
auto x54 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2);
auto x55 = -x54;
auto x56 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2) + x43;
auto x57 = -x56;
auto x58 = cell.template getFieldComponent<momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2) + x28;
auto x59 = V{1} - V{1.5}*x37;
auto x60 = -x24 + x35 + x59;
auto x61 = -x22 + x40;
auto x62 = -x35;
auto x63 = -x58;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x33 + V{1});
auto x1 = -cell[1]*x20 - x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x24 + x32 + x35 - x36) + V{1});
auto x2 = -cell[2]*x20 - x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(-V{4.5}*x37 + x39) + V{1});
auto x3 = -cell[3]*x20 - x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x22 + x32 + x40 - x41) + V{1});
auto x4 = -(cell[4]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x35 + x39 - V{4.5}*x45*x45) + V{1}));
auto x5 = -(cell[5]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x35 + x49 - V{4.5}*x47*x47) + V{1}));
auto x6 = -cell[6]*x20 - x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x40 - x51 + x52) + V{1});
auto x7 = -(cell[7]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x52 + x53 - V{4.5}*x55*x55) + V{1}));
auto x8 = -(cell[8]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x39 + x40 - V{4.5}*x57*x57) + V{1}));
auto x9 = -(cell[9]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x39 + x53 - V{4.5}*x58*x58) + V{1}));
auto x10 = -cell[10]*x20 + x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x36 + x60) + V{-1});
auto x11 = -cell[11]*x20 - x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(-V{3}*x30 + x48) + V{1});
auto x12 = -cell[12]*x20 + x34*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x41 + x59 + x61) + V{-1});
auto x13 = -(cell[13]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x49 + x62 - V{4.5}*x44*x44) + V{1}));
auto x14 = -(cell[14]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x39 + x62 - V{4.5}*x46*x46) + V{1}));
auto x15 = -cell[15]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x51 + x60 + x61) + V{-1});
auto x16 = -(cell[16]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x33 + x40 + x62 - V{4.5}*x54*x54) + V{1}));
auto x17 = -(cell[17]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x49 + x53 - V{4.5}*x56*x56) + V{1}));
auto x18 = -(cell[18]*x20 + x42*(cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*(x40 + x49 - V{4.5}*x63*x63) + V{1}));
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { cell.template getFieldComponent<momenta::FixedDensity::RHO>(0), x21 + x23 + V{1}*x37 };
}
};

}

}

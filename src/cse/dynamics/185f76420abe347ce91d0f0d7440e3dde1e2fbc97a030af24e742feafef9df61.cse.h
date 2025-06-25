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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1};
auto x22 = cell[0] + cell[10] + cell[11] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[18] + cell[1] + cell[2] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + V{2}*cell[8] + V{1};
auto x23 = -x22/x21;
auto x24 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x25 = V{1.5}*x24;
auto x26 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x27 = V{1.5}*x26;
auto x28 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x29 = V{1.5}*x28;
auto x30 = x27 + x29 + V{-1};
auto x31 = x25 + x30;
auto x32 = V{0.0555555555555556}*x19;
auto x33 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x34 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x35 = x25 + V{-1};
auto x36 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x37 = V{0.0277777777777778}*x19;
auto x38 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x39 = x38*x38;
auto x40 = x31 + x33;
auto x41 = -x34;
auto x42 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x43 = -x42;
auto x44 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x45 = x44*x44;
auto x46 = -x36;
auto x47 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x48 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x47;
auto x49 = -x48;
auto x50 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x51 = x50*x50;
auto x52 = x31 + x34;
auto x53 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x47;
auto x54 = -x53;
auto x55 = -x27;
auto x56 = V{1} - x29;
auto x57 = x55 + x56;
auto x58 = x33 + x57;
auto x59 = -x25;
auto x60 = x34 + x59;
auto x61 = x36 + x59;
auto x62 = -x33;
auto x63 = x31 + x36;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x23*x31 + V{1});
auto x1 = -cell[1]*x20 - x32*(-x23*(V{3}*x24 - x30 - x33) + V{1});
auto x2 = -cell[2]*x20 - x32*(-x23*(V{3}*x26 - x29 - x34 - x35) + V{1});
auto x3 = -cell[3]*x20 - x32*(-x23*(-x27 + V{3}*x28 - x35 - x36) + V{1});
auto x4 = -cell[4]*x20 - x37*(-x23*(-x34 + V{4.5}*x39 - x40) + V{1});
auto x5 = -(cell[5]*x20 + x37*(x23*(x40 + x41 - V{4.5}*x43*x43) + V{1}));
auto x6 = -cell[6]*x20 - x37*(-x23*(-x36 - x40 + V{4.5}*x45) + V{1});
auto x7 = -(cell[7]*x20 + x37*(x23*(x40 + x46 - V{4.5}*x49*x49) + V{1}));
auto x8 = -cell[8]*x20 - x37*(-x23*(-x36 + V{4.5}*x51 - x52) + V{1});
auto x9 = -(cell[9]*x20 + x37*(x23*(x46 + x52 - V{4.5}*x54*x54) + V{1}));
auto x10 = -cell[10]*x20 - x32*(-x23*(V{3}*x24 + x58) + V{1});
auto x11 = -cell[11]*x20 - x32*(-x23*(V{3}*x26 + x56 + x60) + V{1});
auto x12 = -cell[12]*x20 - x32*(-x23*(V{3}*x28 + x55 + x61 + V{1}) + V{1});
auto x13 = -cell[13]*x20 - x37*(-x23*(V{4.5}*x39 + x58 + x60) + V{1});
auto x14 = -(cell[14]*x20 + x37*(x23*(x52 + x62 - V{4.5}*x42*x42) + V{1}));
auto x15 = -cell[15]*x20 - x37*(-x23*(V{4.5}*x45 + x58 + x61) + V{1});
auto x16 = -(cell[16]*x20 + x37*(x23*(x62 + x63 - V{4.5}*x48*x48) + V{1}));
auto x17 = -cell[17]*x20 - x37*(-x23*(x36 + V{4.5}*x51 + x57 + x60) + V{1});
auto x18 = -(cell[18]*x20 + x37*(x23*(x41 + x63 - V{4.5}*x53*x53) + V{1}));
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
return { -V{1}*x22/x21, x24 + x26 + x28 };
}
};

}

}

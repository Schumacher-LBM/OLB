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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = (cell[0] + cell[10] + V{2}*cell[11] + cell[12] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[18] + cell[1] + cell[3] + V{2}*cell[5] + cell[6] + cell[7] + V{1})/(cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x22 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x23 = V{1.5}*x22;
auto x24 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x25 = V{1.5}*x24;
auto x26 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x27 = V{1.5}*x26;
auto x28 = x25 + x27 + V{-1};
auto x29 = x23 + x28;
auto x30 = V{0.0555555555555556}*x19;
auto x31 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x32 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x33 = x23 + V{-1};
auto x34 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x35 = V{0.0277777777777778}*x19;
auto x36 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x37 = x36*x36;
auto x38 = x29 + x31;
auto x39 = -x32;
auto x40 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x41 = -x40;
auto x42 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x43 = x42*x42;
auto x44 = -x34;
auto x45 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x46 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x45;
auto x47 = -x46;
auto x48 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x49 = x48*x48;
auto x50 = x29 + x32;
auto x51 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x45;
auto x52 = -x51;
auto x53 = -x25;
auto x54 = V{1} - x27;
auto x55 = x53 + x54;
auto x56 = x31 + x55;
auto x57 = -x23;
auto x58 = x32 + x57;
auto x59 = x34 + x57;
auto x60 = -x31;
auto x61 = x29 + x34;
auto x0 = -cell[0]*x20 - V{0.333333333333333}*x19*(x21*x29 + V{1});
auto x1 = -cell[1]*x20 - x30*(-x21*(V{3}*x22 - x28 - x31) + V{1});
auto x2 = -cell[2]*x20 - x30*(-x21*(V{3}*x24 - x27 - x32 - x33) + V{1});
auto x3 = -cell[3]*x20 - x30*(-x21*(-x25 + V{3}*x26 - x33 - x34) + V{1});
auto x4 = -cell[4]*x20 - x35*(-x21*(-x32 + V{4.5}*x37 - x38) + V{1});
auto x5 = -(cell[5]*x20 + x35*(x21*(x38 + x39 - V{4.5}*x41*x41) + V{1}));
auto x6 = -cell[6]*x20 - x35*(-x21*(-x34 - x38 + V{4.5}*x43) + V{1});
auto x7 = -(cell[7]*x20 + x35*(x21*(x38 + x44 - V{4.5}*x47*x47) + V{1}));
auto x8 = -cell[8]*x20 - x35*(-x21*(-x34 + V{4.5}*x49 - x50) + V{1});
auto x9 = -(cell[9]*x20 + x35*(x21*(x44 + x50 - V{4.5}*x52*x52) + V{1}));
auto x10 = -cell[10]*x20 - x30*(-x21*(V{3}*x22 + x56) + V{1});
auto x11 = -cell[11]*x20 - x30*(-x21*(V{3}*x24 + x54 + x58) + V{1});
auto x12 = -cell[12]*x20 - x30*(-x21*(V{3}*x26 + x53 + x59 + V{1}) + V{1});
auto x13 = -cell[13]*x20 - x35*(-x21*(V{4.5}*x37 + x56 + x58) + V{1});
auto x14 = -(cell[14]*x20 + x35*(x21*(x50 + x60 - V{4.5}*x40*x40) + V{1}));
auto x15 = -cell[15]*x20 - x35*(-x21*(V{4.5}*x43 + x56 + x59) + V{1});
auto x16 = -(cell[16]*x20 + x35*(x21*(x60 + x61 - V{4.5}*x46*x46) + V{1}));
auto x17 = -cell[17]*x20 - x35*(-x21*(x34 + V{4.5}*x49 + x55 + x58) + V{1});
auto x18 = -(cell[18]*x20 + x35*(x21*(x39 + x61 - V{4.5}*x51*x51) + V{1}));
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
return { V{1}*x21, x22 + x24 + x26 };
}
};

}

}

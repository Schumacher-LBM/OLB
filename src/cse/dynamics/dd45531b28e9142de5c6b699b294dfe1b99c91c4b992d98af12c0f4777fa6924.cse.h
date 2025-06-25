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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FreeEnergyMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::FreeEnergy, collision::FreeEnergy, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = parameters.template get<collision::FreeEnergy::GAMMA>();
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x11 = x9 + V{-1};
auto x12 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x13 = x12 + V{1};
auto x14 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x15 = V{1.5}*x14;
auto x16 = cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x17 = V{1.5}*x16;
auto x18 = x17 + V{-1};
auto x19 = x15 + x18;
auto x20 = x12 + V{1};
auto x21 = -V{3}*cell.template getFieldComponent<descriptors::CHEM_POTENTIAL>(0)*x10 + x20;
auto x22 = V{0.0277777777777778}*x9;
auto x23 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x24 = cell.template getFieldComponent<descriptors::FORCE>(0) - cell.template getFieldComponent<descriptors::FORCE>(1);
auto x25 = -x24;
auto x26 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x27 = x19 + x26;
auto x28 = x21*x22;
auto x29 = V{0.111111111111111}*x9;
auto x30 = x21*x29;
auto x31 = -x26;
auto x32 = V{1} - x15;
auto x33 = V{3}*x16 + x32;
auto x34 = cell.template getFieldComponent<descriptors::FORCE>(0) + cell.template getFieldComponent<descriptors::FORCE>(1);
auto x35 = V{4.5}*(x34*x34);
auto x36 = V{3}*x14;
auto x37 = -x17 + x23;
auto x0 = -cell[0]*x11 + V{0.555555555555556}*x21*x9 - V{0.444444444444444}*x9*(x13*x19 + V{1});
auto x1 = -(cell[1]*x11 + x22*(x13*(-x23 + x27 - V{4.5}*x25*x25) + V{1}) + x28);
auto x2 = -cell[2]*x11 - x30 + V{0.111111111111111}*x9*(x13*(x31 + x33) + V{-1});
auto x3 = -cell[3]*x11 - x22*(x13*(x23 + x27 - x35) + V{1}) - x28;
auto x4 = -cell[4]*x11 - x29*(x13*(x18 + x23 - x36) + V{1}) - x30;
auto x5 = -(cell[5]*x11 + x22*(x13*(x19 + x23 + x31 - V{4.5}*x24*x24) + V{1}) + x28);
auto x6 = -cell[6]*x11 - x30 + V{0.111111111111111}*x9*(x13*(x26 + x33) + V{-1});
auto x7 = -cell[7]*x11 - x28 + V{0.0277777777777778}*x9*(x13*(x26 + x32 + x35 + x37) + V{-1});
auto x8 = -cell[8]*x11 - x30 + V{0.111111111111111}*x9*(x13*(x36 + x37 + V{1}) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x20, x14 + x16 };
}
};

}

}

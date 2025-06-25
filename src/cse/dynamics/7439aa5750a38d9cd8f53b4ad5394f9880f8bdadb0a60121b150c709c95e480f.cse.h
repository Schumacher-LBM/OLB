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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::MCGuo<momenta::Identity> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x12 = x11 + V{1};
auto x13 = cell.template getFieldComponent<descriptors::VELOCITY>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x14 = V{1.5}*x13;
auto x15 = cell.template getFieldComponent<descriptors::VELOCITY>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x16 = V{1.5}*x15;
auto x17 = x16 + V{-1};
auto x18 = x14 + x17;
auto x19 = V{0.5}*x9 + V{-1};
auto x20 = cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x21 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x22 = x11 + V{1};
auto x23 = V{0.0277777777777778}*x9;
auto x24 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x25 = cell.template getFieldComponent<descriptors::VELOCITY>(0) - cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x26 = -x25;
auto x27 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x28 = x18 + x27;
auto x29 = V{9}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x30 = V{6}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x31 = V{9}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x32 = V{6}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x33 = x32 + V{3};
auto x34 = x19*x22;
auto x35 = V{0.0277777777777778}*x34;
auto x36 = -x27;
auto x37 = V{1} - x14;
auto x38 = V{3}*x15 + x37;
auto x39 = x30 + V{-3};
auto x40 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x41 = -V{0.333333333333333}*x21;
auto x42 = cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x43 = V{4.5}*(x42*x42);
auto x44 = x32 + V{-3};
auto x45 = V{0.111111111111111}*x9;
auto x46 = V{3}*x13;
auto x47 = V{0.333333333333333}*x20;
auto x48 = V{0.111111111111111}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x49 = x30 + V{3};
auto x50 = -x16 + x24;
auto x0 = -cell[0]*x10 + V{1.33333333333333}*x19*x22*(x20 + x21) - V{0.444444444444444}*x9*(x12*x18 + V{1});
auto x1 = -(cell[1]*x10 + x23*(x12*(-x24 + x28 - V{4.5}*x26*x26) + V{1}) + x35*(-cell.template getFieldComponent<descriptors::FORCE>(0)*(x29 - x30 + V{3}) + cell.template getFieldComponent<descriptors::FORCE>(1)*(-x31 + x33)));
auto x2 = -cell[2]*x10 - x34*(x39*x40 + x41) + V{0.111111111111111}*x9*(x12*(x36 + x38) + V{-1});
auto x3 = -cell[3]*x10 - x23*(x12*(x24 + x28 - x43) + V{1}) - x35*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x29 + x39) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x31 + x44));
auto x4 = -cell[4]*x10 + x19*x22*(-x44*x48 + x47) - x45*(x12*(x17 + x24 - x46) + V{1});
auto x5 = -(cell[5]*x10 + x23*(x12*(x18 + x24 + x36 - V{4.5}*x25*x25) + V{1}) + x35*(cell.template getFieldComponent<descriptors::FORCE>(0)*(-x29 + x49) - cell.template getFieldComponent<descriptors::FORCE>(1)*(x31 - x32 + V{3})));
auto x6 = -cell[6]*x10 - x34*(x40*x49 + x41) + V{0.111111111111111}*x9*(x12*(x27 + x38) + V{-1});
auto x7 = -cell[7]*x10 - x35*(cell.template getFieldComponent<descriptors::FORCE>(0)*(x29 + x49) + cell.template getFieldComponent<descriptors::FORCE>(1)*(x31 + x33)) + V{0.0277777777777778}*x9*(x12*(x27 + x37 + x43 + x50) + V{-1});
auto x8 = -cell[8]*x10 + x34*(-x33*x48 + x47) + x45*(x12*(x46 + x50 + V{1}) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x22, x13 + x15 };
}
};

}

}

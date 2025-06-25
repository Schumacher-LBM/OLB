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
struct CSE<AdvectionDiffusionCornerDynamics2D<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 1, -1>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x5 = parameters.template get<descriptors::OMEGA>();
auto x6 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x7 = x6 + V{-1};
auto x8 = V{0.0555555555555556}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x9 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x10 = x9 + V{-1};
auto x11 = x6 + V{1};
auto x12 = x9 + V{1};
auto x13 = x5 + V{-1};
auto x14 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x15 = V{0.0277777777777778}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x16 = V{0.166666666666667}*cell[0] + x11*x15 + x12*x15 + V{0.0555555555555556};
auto x17 = -x10*x15 - x15*x7 + x16;
auto x18 = x13*(V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x11 - V{1}*cell[3] + V{0.5}*x11*x17 - x14*x7 + V{0.5}*x17*x7 + V{-0.166666666666667});
auto x19 = x10*x17;
auto x20 = V{1}*cell[2];
auto x21 = -x10;
auto x22 = V{0.5}*x15*x21 - V{0.5}*x15*x7 + V{0.5}*x16;
auto x23 = V{0.166666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x0 = V{0.0555555555555556}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x11 + V{0.0555555555555556}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x12 + V{0.333333333333333}*cell[0] - x10*x8 - x7*x8 + V{-0.222222222222222};
auto x1 = -x17*x7 - x18 + V{-0.166666666666667};
auto x2 = x13*(x12*x14 - x12*x22 + x14*x21 - x20 + x21*x22 + V{-0.166666666666667}) - x19 + V{-0.166666666666667};
auto x3 = x11*x17 + x18 + V{-0.166666666666667};
auto x4 = x12*x17 - x13*(V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x12 - x10*x14 - V{0.5}*x12*x17 - V{0.5}*x19 - x20 + V{-0.166666666666667}) + V{-0.166666666666667};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
return { cell[0] - x10*x23 + x11*x23 + x12*x23 - x23*x7 + V{0.333333333333333}, cell.template getFieldComponent<descriptors::VELOCITY>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1) };
}
};

}

}

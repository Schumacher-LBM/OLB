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
struct CSE<AdvectionDiffusionEdgesDynamics<T, descriptors::D3Q7<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q7<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 2, -1, 1>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x7 = parameters.template get<descriptors::OMEGA>();
auto x8 = x7 + V{-1};
auto x9 = V{4}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x10 = V{0.125}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x11 = x10*(x9 + V{-1}) + V{0.125};
auto x12 = x8*(V{1}*cell[1] + x11);
auto x13 = V{4}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x14 = -x10*(x13 + V{1}) + V{0.125};
auto x15 = x8*(V{1}*cell[5] + x14);
auto x16 = V{0.5}*cell[3];
auto x17 = V{4}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x18 = x17 + V{-1};
auto x19 = V{0.0625}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x20 = x17 + V{1};
auto x21 = x19*x20;
cell[0] = V{0.25}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0) + V{-0.25};
cell[1] = -x11 - x12;
cell[2] = -x10*(x13 + V{-1}) + x15 + V{-0.125};
cell[3] = -x10*x18 + x8*(V{0.5}*cell[6] - x16 - x18*x19 - x21) + V{-0.125};
cell[4] = x10*(x9 + V{1}) + x12 + V{-0.125};
cell[5] = -x14 - x15;
cell[6] = V{0.125}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x20 - x8*(V{0.5}*cell[6] - x16 - x18*x19 - x21) + V{-0.125};
return { cell.template getFieldComponent<momenta::FixedDensity::RHO>(0), cell.template getFieldComponent<descriptors::VELOCITY>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1) + cell.template getFieldComponent<descriptors::VELOCITY>(2)*cell.template getFieldComponent<descriptors::VELOCITY>(2) };
}
};

}

}

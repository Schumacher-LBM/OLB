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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::HeatFluxBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumAD, momenta::NoStress, momenta::DefineToEq> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x5 = parameters.template get<descriptors::OMEGA>();
auto x6 = V{1} / (cell.template getFieldComponent<descriptors::VELOCITY>(1) + V{1});
auto x7 = cell.template getFieldComponent<momenta::FixedVelocityMomentumAD::VELOCITY>(1)*x6;
auto x8 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x9 = x8 + V{1};
auto x10 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumAD::VELOCITY>(1) + cell[0] + cell[1] + cell[3] + V{2}*cell[4] + V{1};
auto x11 = -x10;
auto x12 = x11*x6;
auto x13 = V{0.0555555555555556}*x12;
auto x14 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x15 = x14 + V{1};
auto x16 = x8 + V{-1};
auto x17 = x14 + V{-1};
auto x18 = V{0.0277777777777778}*x12;
auto x19 = V{0.0277777777777778}*x9;
auto x20 = x6*(V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[3] + V{0.111111111111111}*cell[4] + V{0.0555555555555556}) - V{0.0555555555555556}*x7 + V{1.85037170770859e-17};
auto x21 = -x12*x19 - x15*x18 + x16*x18 + x17*x18 + x20;
auto x22 = x16*x21;
auto x23 = x5 + V{-1};
auto x24 = -x16;
auto x25 = x10*x6;
auto x26 = V{0.0277777777777778}*x25;
auto x27 = -x17;
auto x28 = V{0.5}*x15*x26 + V{0.5}*x19*x25 + V{0.5}*x20 + V{0.5}*x24*x26 + V{0.5}*x26*x27;
auto x29 = V{0.0833333333333333}*x25;
auto x30 = x17*x21;
auto x31 = V{0.0833333333333333}*x12;
auto x32 = V{0.166666666666667}*x12;
auto x0 = V{0.0555555555555556}*x11*x16*x6 + V{0.0555555555555556}*x11*x17*x6 - x13*x15 - x13*x9 + x6*(V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[3] + V{0.222222222222222}*cell[4] + V{0.111111111111111}) - V{0.111111111111111}*x7 + V{-0.333333333333333};
auto x1 = -x22 + x23*(V{0.5}*cell.template getFieldComponent<momenta::FixedVelocityMomentumAD::VELOCITY>(0) + x24*x28 - x24*x29 - x28*x9 + x29*x9) + V{-0.166666666666667};
auto x2 = x23*(V{0.5}*cell.template getFieldComponent<momenta::FixedVelocityMomentumAD::VELOCITY>(1) - x15*x28 + x15*x29 + x27*x28 - x27*x29) - x30 + V{-0.166666666666667};
auto x3 = x21*x9 - x23*(V{0.5}*cell.template getFieldComponent<momenta::FixedVelocityMomentumAD::VELOCITY>(0) - x16*x31 - V{0.5}*x21*x9 - V{0.5}*x22 - x31*x9) + V{-0.166666666666667};
auto x4 = x15*x21 - x23*(V{0.5}*cell.template getFieldComponent<momenta::FixedVelocityMomentumAD::VELOCITY>(1) - V{0.5}*x15*x21 - x15*x31 - x17*x31 - V{0.5}*x30) + V{-0.166666666666667};
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
return { -x15*x32 + x16*x32 + x17*x32 - x32*x9 + x6*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[3] + V{0.666666666666667}*cell[4] + V{0.333333333333333}) - V{0.333333333333333}*x7 + V{1.11022302462516e-16}, cell.template getFieldComponent<descriptors::VELOCITY>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1) };
}
};

}

}

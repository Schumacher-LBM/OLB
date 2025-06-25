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
struct CSE<SourcedAdvectionDiffusionBGKdynamics<T, descriptors::D3Q19<FIELDS...> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{0.5}*x19 + V{-1};
auto x22 = cell.template getFieldComponent<descriptors::SOURCE>(0)*x21;
auto x23 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x24 = x23 + V{-1};
auto x25 = V{0.0277777777777778}*cell.template getFieldComponent<descriptors::SOURCE>(0);
auto x26 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + x25 + V{0.0555555555555556};
auto x27 = V{0.0555555555555556}*x22;
auto x28 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x29 = x28 + V{-1};
auto x30 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(2);
auto x31 = V{0.0138888888888889}*cell.template getFieldComponent<descriptors::SOURCE>(0) + V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x32 = x21*x25;
auto x33 = -x23;
auto x34 = x28 + V{1};
auto x35 = x30 + V{1};
auto x36 = -x28;
auto x37 = x23 + V{1};
auto x38 = -x30;
auto x0 = -cell[0]*x20 + x19*(V{0.166666666666667}*cell.template getFieldComponent<descriptors::SOURCE>(0) + V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9]) - V{0.333333333333333}*x22;
auto x1 = -cell[1]*x20 - x19*(x24*x26 + V{0.0555555555555556}) - x27;
auto x2 = -cell[2]*x20 - x19*(x26*x29 + V{0.0555555555555556}) - x27;
auto x3 = -cell[3]*x20 - x19*(x26*(x30 + V{-1}) + V{0.0555555555555556}) - x27;
auto x4 = -cell[4]*x20 - x19*(x31*(x24 + x28) + V{0.0277777777777778}) - x32;
auto x5 = -cell[5]*x20 + x19*(x31*(x33 + x34) + V{-0.0277777777777778}) - x32;
auto x6 = -cell[6]*x20 - x19*(x31*(x24 + x30) + V{0.0277777777777778}) - x32;
auto x7 = -cell[7]*x20 + x19*(x31*(x33 + x35) + V{-0.0277777777777778}) - x32;
auto x8 = -cell[8]*x20 - x19*(x31*(x29 + x30) + V{0.0277777777777778}) - x32;
auto x9 = -cell[9]*x20 + x19*(x31*(x35 + x36) + V{-0.0277777777777778}) - x32;
auto x10 = -cell[10]*x20 + x19*(x26*x37 + V{-0.0555555555555556}) - x27;
auto x11 = -cell[11]*x20 + x19*(x26*x34 + V{-0.0555555555555556}) - x27;
auto x12 = -cell[12]*x20 + x19*(x26*x35 + V{-0.0555555555555556}) - x27;
auto x13 = -cell[13]*x20 + x19*(x31*(x28 + x37) + V{-0.0277777777777778}) - x32;
auto x14 = -cell[14]*x20 + x19*(x31*(x36 + x37) + V{-0.0277777777777778}) - x32;
auto x15 = -cell[15]*x20 + x19*(x31*(x30 + x37) + V{-0.0277777777777778}) - x32;
auto x16 = -cell[16]*x20 + x19*(x31*(x37 + x38) + V{-0.0277777777777778}) - x32;
auto x17 = -cell[17]*x20 + x19*(x31*(x30 + x34) + V{-0.0277777777777778}) - x32;
auto x18 = -cell[18]*x20 + x19*(x31*(x34 + x38) + V{-0.0277777777777778}) - x32;
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
return { V{0.5}*cell.template getFieldComponent<descriptors::SOURCE>(0) + cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1}, cell.template getFieldComponent<descriptors::VELOCITY>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1) + cell.template getFieldComponent<descriptors::VELOCITY>(2)*cell.template getFieldComponent<descriptors::VELOCITY>(2) };
}
};

}

}

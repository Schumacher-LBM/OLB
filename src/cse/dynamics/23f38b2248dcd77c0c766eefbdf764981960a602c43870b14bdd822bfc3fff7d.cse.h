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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkPressure<momenta::ForcedMomentum<momenta::IncompressibleBulkMomentum> >, momenta::IncompressibleBulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::MPIncompressible, collision::OmegaFromCellTauEff<collision::BGK>, forcing::Liang<momenta::Forced> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = V{1} / (cell.template getFieldComponent<descriptors::TAU_EFF>(0));
auto x10 = V{1} - x9;
auto x11 = V{1}*cell[1];
auto x12 = V{1}*cell[4];
auto x13 = V{1}*cell[5];
auto x14 = V{1}*cell[6];
auto x15 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x16 = V{1} / (cell.template getFieldComponent<descriptors::RHO>(0));
auto x17 = V{1}*cell[7];
auto x18 = x11 - x13;
auto x19 = V{1}*cell[2];
auto x20 = V{1}*cell[3];
auto x21 = x19 + x20;
auto x22 = -x14 - x17 + x18 + x21;
auto x23 = x16*x22;
auto x24 = x15 - x23;
auto x25 = cell.template getFieldComponent<descriptors::NABLARHO>(0)*x24;
auto x26 = V{0.5}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x27 = V{1}*cell[8];
auto x28 = x17 + x27;
auto x29 = x16*(-x12 + x18 - x20 + x28);
auto x30 = x26 + x29;
auto x31 = cell.template getFieldComponent<descriptors::NABLARHO>(1)*x30;
auto x32 = V{0.666666666666667}*cell.template getFieldComponent<descriptors::RHO>(0);
auto x33 = x24*x24;
auto x34 = x30*x30;
auto x35 = x33 + x34;
auto x36 = x15 - x16*x22;
auto x37 = V{0.0833333333333333}*cell.template getFieldComponent<descriptors::RHO>(0);
auto x38 = -x26 - x29 + x36;
auto x39 = V{0.75}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x40 = V{1.5}*cell[3];
auto x41 = V{1.5}*cell[1] - V{1.5}*cell[5];
auto x42 = x16*(-V{1.5}*cell[4] + V{1.5}*cell[7] + V{1.5}*cell[8] - x40 + x41);
auto x43 = V{0.75}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x44 = V{1.5}*cell[2] - V{1.5}*cell[6] - V{1.5}*cell[7] + x40 + x41;
auto x45 = -x16*x44 + x43;
auto x46 = -x15 + x23;
auto x47 = x30 + x46;
auto x48 = V{0.5}*x33;
auto x49 = V{0.5}*x34;
auto x50 = -x48 - x49;
auto x51 = -V{0.0333333333333333}*cell.template getFieldComponent<descriptors::RHO>(0)*x35 + V{0.05}*cell[1] + V{0.05}*cell[2] + V{0.05}*cell[3] + V{0.05}*cell[4] + V{0.05}*cell[5] + V{0.05}*cell[6] + V{0.05}*cell[7] + V{0.05}*cell[8] + V{0.025}*x25 + V{0.025}*x31;
auto x52 = cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::RHO>(0);
auto x53 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<descriptors::RHO>(0);
auto x54 = -x47;
auto x55 = V{1} - V{0.5}*x9;
auto x56 = V{0.0833333333333333}*x55;
auto x57 = V{0.333333333333333}*cell.template getFieldComponent<descriptors::RHO>(0);
auto x58 = x48 + x49;
auto x59 = -V{0.133333333333333}*cell.template getFieldComponent<descriptors::RHO>(0)*x35 + V{0.2}*cell[1] + V{0.2}*cell[2] + V{0.2}*cell[3] + V{0.2}*cell[4] + V{0.2}*cell[5] + V{0.2}*cell[6] + V{0.2}*cell[7] + V{0.2}*cell[8] + V{0.1}*x25 + V{0.1}*x31;
auto x60 = V{0.333333333333333}*x55;
auto x61 = x39 + x42;
auto x62 = x24 + x30;
auto x63 = cell.template getFieldComponent<descriptors::NABLARHO>(0)*x62;
auto x64 = cell.template getFieldComponent<descriptors::NABLARHO>(1)*x62;
auto x65 = x52 + x53;
auto x66 = x30*x61;
auto x67 = x16*x44;
auto x68 = x43 - x67;
auto x0 = V{1}*(cell[0]*x10 - x9*(x11 + x12 + x13 + x14 + x21 + V{0.5}*x25 + x28 + V{0.5}*x31 + x32*x35 - x32*(x34 + x36*x36)));
auto x1 = x10*x11 - x56*(-cell.template getFieldComponent<descriptors::NABLARHO>(0)*x54 + cell.template getFieldComponent<descriptors::NABLARHO>(1)*x38 + x52 - x53) + x9*(-x37*(-x38*(-x39 - x42 + x45) - x47 - x50) + x51);
auto x2 = x10*x19 - x60*(-x25 + x52) + x9*(-x57*(x24 - x36*x45 + x58) + x59);
auto x3 = x10*x20 - x56*(-x63 - x64 + x65) + x9*(-x37*(x58 + x62 - (x30 + x36)*(x45 + x61)) + x51);
auto x4 = x10*x12 - x60*(-x31 + x53) + x9*(-x57*(x30 + x58 - x66) + x59);
auto x5 = x10*x13 + x56*(cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::RHO>(0) - cell.template getFieldComponent<descriptors::NABLARHO>(0)*x47 - cell.template getFieldComponent<descriptors::NABLARHO>(1)*x54 - x53) + x9*(-x37*(x47 - x54*(V{0.75}*cell.template getFieldComponent<descriptors::FORCE>(0) - x61 - x67) + x58) + x51);
auto x6 = x10*x14 + x60*(x25 + x52) + x9*(-x57*(-x24*x68 + x46 + x58) + x59);
auto x7 = x10*x17 + x56*(x63 + x64 + x65) + x9*(x37*(x50 + x62*(x61 + x68) + x62) + x51);
auto x8 = x10*x27 + x60*(x31 + x53) + x9*(x57*(x30 + x50 + x66) + x59);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { cell.template getFieldComponent<descriptors::RHO>(0), x35 };
}
};

}

}

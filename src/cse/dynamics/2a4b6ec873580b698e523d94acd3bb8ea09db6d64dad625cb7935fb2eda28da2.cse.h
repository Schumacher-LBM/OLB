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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FreeEnergyMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::FreeEnergy, collision::FreeEnergy, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = parameters.template get<collision::FreeEnergy::GAMMA>();
auto x21 = x19 + V{-1};
auto x22 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9];
auto x23 = x22 + V{1};
auto x24 = cell.template getFieldComponent<descriptors::FORCE>(0)*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x25 = V{1.5}*x24;
auto x26 = cell.template getFieldComponent<descriptors::FORCE>(1)*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x27 = V{1.5}*x26;
auto x28 = cell.template getFieldComponent<descriptors::FORCE>(2)*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x29 = V{1.5}*x28;
auto x30 = x27 + x29 + V{-1};
auto x31 = x25 + x30;
auto x32 = x22 + V{1};
auto x33 = -V{3}*cell.template getFieldComponent<descriptors::CHEM_POTENTIAL>(0)*x20 + x32;
auto x34 = V{0.0555555555555556}*x19;
auto x35 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(0);
auto x36 = V{3}*x24;
auto x37 = x33*x34;
auto x38 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(1);
auto x39 = V{3}*x26;
auto x40 = x25 + V{-1};
auto x41 = V{3}*cell.template getFieldComponent<descriptors::FORCE>(2);
auto x42 = V{3}*x28;
auto x43 = V{0.0277777777777778}*x19;
auto x44 = cell.template getFieldComponent<descriptors::FORCE>(0) + cell.template getFieldComponent<descriptors::FORCE>(1);
auto x45 = V{4.5}*(x44*x44);
auto x46 = x31 + x35;
auto x47 = x33*x43;
auto x48 = -x38;
auto x49 = cell.template getFieldComponent<descriptors::FORCE>(0) - cell.template getFieldComponent<descriptors::FORCE>(1);
auto x50 = -x49;
auto x51 = cell.template getFieldComponent<descriptors::FORCE>(0) + cell.template getFieldComponent<descriptors::FORCE>(2);
auto x52 = V{4.5}*(x51*x51);
auto x53 = -x41;
auto x54 = -cell.template getFieldComponent<descriptors::FORCE>(2);
auto x55 = cell.template getFieldComponent<descriptors::FORCE>(0) + x54;
auto x56 = -x55;
auto x57 = cell.template getFieldComponent<descriptors::FORCE>(1) + cell.template getFieldComponent<descriptors::FORCE>(2);
auto x58 = V{4.5}*(x57*x57);
auto x59 = x31 + x38;
auto x60 = cell.template getFieldComponent<descriptors::FORCE>(1) + x54;
auto x61 = -x60;
auto x62 = -x27;
auto x63 = V{1} - x29;
auto x64 = x62 + x63;
auto x65 = x35 + x64;
auto x66 = -x25;
auto x67 = x38 + x66;
auto x68 = x41 + x66;
auto x69 = -x35;
auto x70 = x31 + x41;
auto x0 = -cell[0]*x21 + V{0.666666666666667}*x19*x33 - V{0.333333333333333}*x19*(x23*x31 + V{1});
auto x1 = -cell[1]*x21 - x34*(x23*(x30 + x35 - x36) + V{1}) - x37;
auto x2 = -cell[2]*x21 - x34*(x23*(x29 + x38 - x39 + x40) + V{1}) - x37;
auto x3 = -cell[3]*x21 - x34*(x23*(x27 + x40 + x41 - x42) + V{1}) - x37;
auto x4 = -cell[4]*x21 - x43*(x23*(x38 - x45 + x46) + V{1}) - x47;
auto x5 = -(cell[5]*x21 + x43*(x23*(x46 + x48 - V{4.5}*x50*x50) + V{1}) + x47);
auto x6 = -cell[6]*x21 - x43*(x23*(x41 + x46 - x52) + V{1}) - x47;
auto x7 = -(cell[7]*x21 + x43*(x23*(x46 + x53 - V{4.5}*x56*x56) + V{1}) + x47);
auto x8 = -cell[8]*x21 - x43*(x23*(x41 - x58 + x59) + V{1}) - x47;
auto x9 = -(cell[9]*x21 + x43*(x23*(x53 + x59 - V{4.5}*x61*x61) + V{1}) + x47);
auto x10 = -cell[10]*x21 + V{0.0555555555555556}*x19*(x23*(x36 + x65) + V{-1}) - x37;
auto x11 = -cell[11]*x21 + V{0.0555555555555556}*x19*(x23*(x39 + x63 + x67) + V{-1}) - x37;
auto x12 = -cell[12]*x21 + V{0.0555555555555556}*x19*(x23*(x42 + x62 + x68 + V{1}) + V{-1}) - x37;
auto x13 = -cell[13]*x21 + V{0.0277777777777778}*x19*(x23*(x45 + x65 + x67) + V{-1}) - x47;
auto x14 = -(cell[14]*x21 + x43*(x23*(x59 + x69 - V{4.5}*x49*x49) + V{1}) + x47);
auto x15 = -cell[15]*x21 + V{0.0277777777777778}*x19*(x23*(x52 + x65 + x68) + V{-1}) - x47;
auto x16 = -(cell[16]*x21 + x43*(x23*(x69 + x70 - V{4.5}*x55*x55) + V{1}) + x47);
auto x17 = -cell[17]*x21 + V{0.0277777777777778}*x19*(x23*(x41 + x58 + x64 + x67) + V{-1}) - x47;
auto x18 = -(cell[18]*x21 + x43*(x23*(x48 + x70 - V{4.5}*x60*x60) + V{1}) + x47);
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
return { x32, x24 + x26 + x28 };
}
};

}

}

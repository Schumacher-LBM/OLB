/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Luiz Eduardo Czelusniak
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

#ifndef PHASE_FIELD_BOUNDARY_DYNAMICS_H
#define PHASE_FIELD_BOUNDARY_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

/**
* Implementation of Dirichlet boundary condition for the order parameter.
*/
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int direction, int orientation>
class PhaseFieldInletDynamics : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:
  // Use same MOMENTA in combined and nested (boundary) dynamics
  using CORRECTED_DYNAMICS = typename DYNAMICS::template exchange_momenta<MOMENTA>;

public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename CORRECTED_DYNAMICS::EquilibriumF;

  using parameters = typename CORRECTED_DYNAMICS::parameters;

  template<typename M>
  using exchange_momenta = PhaseFieldInletDynamics<T,DESCRIPTOR,DYNAMICS,M,direction,orientation>;

  std::type_index id() override {
    return typeid(PhaseFieldInletDynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<PhaseFieldInletDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    // Along all the commented parts of this code there will be an example based
    // on the situation where the wall's normal vector if (0,1) and the
    // numerotation of the velocites are done according to the D2Q9
    // lattice of the OpenLB library.

    // Find all the missing populations
    // (directions 3,4,5)
    constexpr auto missingIndices = util::subIndexOutgoing<DESCRIPTOR,direction,orientation>();

    V phi = MomentaF().computeRho(cell);
    V missingPhi = phi - 1.;
    V missingWeightSum = 0;
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {

      bool contains = false;
      for(unsigned i=0; i < missingIndices.size(); i++){
        if(missingIndices[i] == (unsigned)iPop){
          contains = true;
        }
      }

      if(contains){
        missingWeightSum += descriptors::t<V,DESCRIPTOR>(iPop);
      } else {
        missingPhi -= cell[iPop];
      }
    }

    for (unsigned iPop=0; iPop < missingIndices.size(); ++iPop) {
      cell[missingIndices[iPop]] = missingPhi * descriptors::t<V,DESCRIPTOR>(missingIndices[iPop]) / missingWeightSum;
    }

    return typename CORRECTED_DYNAMICS::CollisionO().apply(cell, parameters);
  }

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "PhaseFieldInletDynamics<" + CORRECTED_DYNAMICS().getName() + ">";
  };

};

/**
* Implementation of convective boundary condition for the order parameter.
*/
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int direction, int orientation>
class PhaseFieldConvectiveOutletDynamics : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:

public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename DYNAMICS::EquilibriumF;

  using parameters = meta::list<descriptors::MAX_VELOCITY>;

  std::type_index id() override {
    return typeid(PhaseFieldConvectiveOutletDynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<PhaseFieldConvectiveOutletDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {

    auto u_conv = parameters.template get<descriptors::MAX_VELOCITY>();
    Vector<int,DESCRIPTOR::d> normal;
    normal[direction] = -orientation;
    auto cell1 = cell.neighbor(normal);
    auto outlet_cell = cell.template getField<descriptors::CONV_POPS>();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      outlet_cell[iPop] = (outlet_cell[iPop] + u_conv * cell1[iPop]) / ((T)1 + u_conv);
      cell[iPop] = outlet_cell[iPop];
    }

    V phi = MomentaF().computeRho(cell);
    cell.template setField<descriptors::CONV_POPS>(outlet_cell);
    V u[DESCRIPTOR::d];
    MomentaF().computeU(cell, u);
    V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    return {phi, uSqr};
  }

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "PhaseFieldConvectiveOutletDynamics<" + DYNAMICS().getName() + ">";
  };

};

/**
* Implementation of convective boundary condition.
*/
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int direction, int orientation>
class ConvectiveOutletDynamics : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:

public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename DYNAMICS::EquilibriumF;

  using parameters = meta::list<descriptors::MAX_VELOCITY>;

  std::type_index id() override {
    return typeid(ConvectiveOutletDynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<ConvectiveOutletDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {

    auto u_conv = parameters.template get<descriptors::MAX_VELOCITY>();
    Vector<int,DESCRIPTOR::d> normal;
    normal[direction] = -orientation;
    auto cell1 = cell.neighbor(normal);
    auto outlet_cell = cell.template getField<descriptors::CONV_POPS>();
    auto rho = cell.template getField<descriptors::RHO>();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      outlet_cell[iPop] = (outlet_cell[iPop] + u_conv * cell1[iPop]) / ((T)1 + u_conv);
      cell[iPop] = outlet_cell[iPop];
    }

    cell.template setField<descriptors::CONV_POPS>(outlet_cell);
    V u[DESCRIPTOR::d];
    MomentaF().computeU(cell, u);
    V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    return {rho, uSqr};
  }

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "ConvectiveOutletDynamics<" + DYNAMICS().getName() + ">";
  };

};

}

#endif

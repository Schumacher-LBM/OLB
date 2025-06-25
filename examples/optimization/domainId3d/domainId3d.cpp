/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>
#include "helper.h"  // Will be removed once SuperLatticeFieldReductionO enables indicator support

#include "domainId3d.h"

using DUAL_DYNAMICS = DualPorousBGKDynamics<T,DESCRIPTOR>;     // Dynamics of dual problem
using OBJECTIVE = functors::L2DistanceF<functors::VelocityF>;  // Objective functional
using CONTROLS = descriptors::POROSITY;                        // Controlled field

void prepareReference(auto& reference) {
  prepareLattice(reference->getUnitConverter(), reference->getSuperLattice(), reference->getSuperGeometry());
  reference->getLatticeResults().template add<SuperLatticePhysPressure3D<T,DESCRIPTOR>,
                                              SuperLatticePhysVelocity3D<T,DESCRIPTOR>,
                                              SuperLatticePorosity3D<T,DESCRIPTOR>>();
}

void preparePrimal(auto& primal, const std::vector<T>& controls) {
  auto& converter = primal->getUnitConverter();
  auto& superGeometry = primal->getSuperGeometry();
  auto& sLattice = primal->getSuperLattice();
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,6});

  // Define primal physics
  sLattice.template defineDynamics<DYNAMICS>(bulkIndicator);
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  boundary::set<boundary::LocalVelocity>(sLattice, superGeometry.getMaterialIndicator({3}));
  boundary::set<boundary::LocalPressure>(sLattice, superGeometry.getMaterialIndicator({4}));
  sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  // Initialize primal problem
  AnalyticalConst3D<T,T> rhoF(1);
  Vector<T,3> velocity;
  AnalyticalConst3D<T,T> uF(velocity);
  sLattice.defineRhoU(bulkIndicator, rhoF, uF);
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);
  AnalyticalConst3D<T,T> one(1.);
  sLattice.template defineField<CONTROLS>(bulkIndicator, one);

  // Convert control vector to porosity field
  std::vector<T> porosity = opti::applyProjection<opti::projection::Sigmoid<T>>(controls);
  setFieldFromSerialized<CONTROLS>(porosity, sLattice, superGeometry.getMaterialIndicator({6}));
  sLattice.initialize();

  primal->getLatticeResults().template add<SuperLatticePhysPressure3D<T,DESCRIPTOR>,
                                           SuperLatticePhysVelocity3D<T,DESCRIPTOR>,
                                           SuperLatticePorosity3D<T,DESCRIPTOR>>();
}

T computeObjective(auto& reference, auto& primal) {
  const auto& converter = reference->getUnitConverter();
  auto& superGeometry = reference->getSuperGeometry();
  auto& refLattice = reference->getSuperLattice();
  auto& primalLattice = primal->getSuperLattice();
  auto objectiveDomain = superGeometry.getMaterialIndicator({1,2,3,4,6});

  // Evaluate functor for objective computation
  auto objectiveO = makeWriteFunctorO<OBJECTIVE,opti::J>(primalLattice);
  objectiveO->restrictTo(objectiveDomain);
  writePhysFunctorTo<functors::VelocityF,OBJECTIVE::Reference>(refLattice,
                                                               objectiveDomain,
                                                               converter.getConversionFactorVelocity());
  // Get solution from the reference simulation for the inverse problem
  copyFields<OBJECTIVE::Reference,OBJECTIVE::Reference>(refLattice, primalLattice);
  objectiveO->template setParameter<descriptors::CONVERSION>(converter.getConversionFactorVelocity());
  objectiveO->template setParameter<descriptors::NORMALIZE>(norm(refLattice, converter, objectiveDomain));
  objectiveO->execute();

  // Compute source term for the dual simulation
  auto objectiveDerivativeO = makeWriteFunctorO<functors::DerivativeF<OBJECTIVE,descriptors::POPULATION,DYNAMICS>,
                                                opti::DJDF>(primalLattice);
  objectiveDerivativeO->restrictTo(objectiveDomain);
  objectiveDerivativeO->template setParameter<descriptors::CONVERSION>(converter.getConversionFactorVelocity());
  objectiveDerivativeO->template setParameter<descriptors::NORMALIZE>(norm(refLattice, converter, objectiveDomain));
  objectiveDerivativeO->template setParameter<descriptors::DX>(converter.getPhysDeltaX());

  objectiveDerivativeO->execute();
  return integrate<opti::J>(primalLattice, objectiveDomain)[0];
}

void prepareDual(auto& dual, auto& primal, auto& reference) {
  auto& converter = dual->getUnitConverter();
  auto& superGeometry = dual->getSuperGeometry();
  auto& sLattice = dual->getSuperLattice();
  auto& primalLattice = primal->getSuperLattice();
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,6});

  // Define dual physics
  sLattice.template defineDynamics<DUAL_DYNAMICS>(bulkIndicator);
  boundary::set<boundary::BounceBack>(sLattice, superGeometry.getMaterialIndicator({2,3,4}));
  sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  // Initialize dual problem
  AnalyticalConst3D<T,T> rhoF(1);
  Vector<T,3> velocity;
  AnalyticalConst3D<T,T> uF(velocity);
  sLattice.defineRhoU(bulkIndicator, rhoF, uF);
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uF);

  // Provide fields required by the dual collision operator
  copyFields<CONTROLS,CONTROLS>(primalLattice, sLattice);
  writeFunctorTo<functors::PopulationF,opti::F>(primalLattice, bulkIndicator);
  copyFields<opti::F,opti::F>(primalLattice, sLattice);
  copyFields<opti::DJDF,opti::DJDF>(primalLattice, sLattice);
  sLattice.initialize();

  dual->getLatticeResults().template add<SuperLatticePhysPressure3D<T,DESCRIPTOR>,
                                         SuperLatticePhysVelocity3D<T,DESCRIPTOR>,
                                         SuperLatticePorosity3D<T,DESCRIPTOR>,
                                         SuperLatticeField3D<T,DESCRIPTOR,opti::DJDF>>();
}

std::vector<T> computeSensitivity(auto& dual, auto& primal, const std::vector<T>& controls) {
  const auto& converter = primal->getUnitConverter();
  auto& superGeometry = primal->getSuperGeometry();
  auto& primalLattice = primal->getSuperLattice();
  auto& dualLattice = dual->getSuperLattice();

  // Evaluate optimality condition
  auto optimalityO = makeWriteFunctorO<functors::OptimalityF<DYNAMICS,CONTROLS>,
                                       opti::SENSITIVITY<CONTROLS>>(dualLattice);
  // Derivative of applied projection
  std::vector<T> projectionD = opti::applyProjection<opti::projection::SigmoidD<T>>(controls);
  setFieldFromSerialized<opti::DPROJECTIONDALPHA<CONTROLS>>(projectionD, dualLattice, superGeometry.getMaterialIndicator({6}));
  // Compute jacobian of collision operator regarding control variable
  auto dCDalphaO = makeWriteFunctorO<functors::DerivativeF<functors::CollisionF<DYNAMICS>,CONTROLS,DYNAMICS>,
                                     opti::DCDALPHA<CONTROLS>>(primalLattice);
  dCDalphaO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  dCDalphaO->template setParameter<descriptors::DX>(1.0);
  dCDalphaO->execute();
  // Jacobian is computed on primal lattice as jacobian is evaluated for primal populations
  copyFields<opti::DCDALPHA<CONTROLS>,opti::DCDALPHA<CONTROLS>>(primalLattice, dualLattice);
  optimalityO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  optimalityO->execute();

  // Return serial vector containing total derivatives of objective regarding controls
  return getSerializedFromField<opti::SENSITIVITY<CONTROLS>>(dualLattice, superGeometry.getMaterialIndicator({6}));
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  // Create Simulation instances
  auto refSimulation = build(" REFERENCE ");
  auto primalSimulation = build(" PRIMAL ");
  auto dualSimulation = build(" DUAL ");

  // Create instances for IO for each simulation
  refSimulation->template create<LatticeResults<T,DESCRIPTOR>>("refSolution");
  primalSimulation->template create<LatticeResults<T,DESCRIPTOR>>("primalSolution");
  dualSimulation->template create<LatticeResults<T,DESCRIPTOR>>("dualSolution");

  // Simulate reference solution
  prepareReference(refSimulation);
  simulate(refSimulation);

  // Define routine for objective and derivative computation
  opti::OptiCaseAnalytical<T,std::vector<T>> optiCase;
  optiCase.setObjective([&](const std::vector<T>& controls) -> T {
    primalSimulation->resetLattice();
    preparePrimal(primalSimulation, controls);
    simulate(primalSimulation);
    return computeObjective(refSimulation, primalSimulation);
  });
  optiCase.setDerivative([&](const std::vector<T>& controls, std::vector<T>& derivatives) {
    dualSimulation->resetLattice();
    prepareDual(dualSimulation, primalSimulation, refSimulation);
    simulate(dualSimulation);
    derivatives = computeSensitivity(dualSimulation, primalSimulation, controls);
  });

  // Optimizer setup
  const std::size_t dimCtrl = getSerializedFieldSize<CONTROLS>(primalSimulation->getSuperLattice(),
                                                               primalSimulation->getSuperGeometry().getMaterialIndicator({6}));
  opti::OptimizerLBFGS<T,std::vector<T>> optimizer(
      dimCtrl, 1.e-10, 10, 1., 20, "Wolfe", 20, 1.e-4);
  optimizer.setStartValue(opti::projection::getInitialControl(1e-2,
    opti::projection::Sigmoid<T>(),
    refSimulation->getUnitConverter(),
    opti::Permeability));
  optimizer.optimize(optiCase);
}

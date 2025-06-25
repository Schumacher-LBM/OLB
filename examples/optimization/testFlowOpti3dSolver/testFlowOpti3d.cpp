/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger
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

/** \file A fluid flow simulation test case - execution script
 * Cf. testFlowOpti3d.h for implementation details
 */

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include "testFlowOpti3d.h"

using namespace olb;
using namespace olb::opti;
using namespace olb::parameters;

using S = FLOATING_POINT_TYPE;
constexpr unsigned numberOfDerivatives (3);
using T = util::ADf<S,numberOfDerivatives>;


/** Compute sensitivities with forward difference quotients
 * "control" variable scales the force, objective measures distance of velocity
 * or dissipation to the exact solution
 */
void sensitivitiesFdq() {
  OstreamManager clout(std::cout, "sensitivitiesFdq");
  XMLreader config("parameterAd.xml");

  auto testFlow = createLbSolver <TestFlowSolverDirectOpti<S>> (config);
  testFlow->parameters(Output()).verbose = false;
  OptiCaseFDQ optiCase(getCallable<S>(testFlow));

  // evaluate function
  std::vector<S> control (numberOfDerivatives, 2.5);
  clout << "objective = " << optiCase.evaluateObjective(control) << std::endl;

  // compute derivatives
  std::vector<S> derivatives (numberOfDerivatives, 12.);
  optiCase.computeDerivatives(control, derivatives);
  util::print(derivatives, "derivatives", clout, ' ');
}

/** Compute sensitivities with automatic differentiation
 * "control" variable scales the force, objective measures distance of velocity
 * or dissipation to the exact solution
 */
void sensitivitiesAd() {
  // preparation
  OstreamManager clout(std::cout, "sensitivitiesAd");
  XMLreader config("parameterAd.xml");

  auto testFlow = createLbSolver <TestFlowSolverDirectOpti<S>> (config);
  auto testFlowAD = createLbSolver <TestFlowSolverDirectOpti<T>> (config);
  testFlow->parameters(Output()).verbose = false;
  testFlowAD->parameters(Output()).verbose = false;
  testFlowAD->parameters(VisualizationVTK()).output = "off";

  OptiCaseAdForSolver optiCase(testFlow, testFlowAD);

  // evaluate function
  std::vector<S> control (numberOfDerivatives, 2.5);
  util::print(control, "control", clout, ' ');
  clout << "objective = " << optiCase.evaluateObjective(control) << std::endl;

  // compute (forward AD) derivatives
  std::vector<S> derivatives (numberOfDerivatives, 12.);
  optiCase.computeDerivatives(control, derivatives);
  util::print(derivatives, "derivatives", clout, ' ');
}

/** Variant for evaluation + forward AD sensitivity computation with arbitrary
 * number of control variables: "control" scales the force cuboid-wise
 * (independent of the cuboid decomposition which is used for parallelization).
 */
void sensitivitiesAd_variant(unsigned numberOfControls)
{
  // preparation
  OstreamManager clout(std::cout, "sensitivitiesAd");
  XMLreader config("parameterAd.xml");

  std::vector<S> control (numberOfControls, 2.5);
  util::print(control, "control", clout, ' ');

  // evaluate function
  auto testFlow = createLbSolver <TestFlowSolverDirectOpti<S>> (config);
  auto function = getCallable<S> (testFlow);
  clout << "objective = " << function(control, 0) << std::endl;

  // compute (forward AD) derivatives
  auto testFlowAD = createLbSolver <TestFlowSolverDirectOpti<T>> (config);
  auto derivative = getCallable<T> (testFlowAD);

  auto controlAD = util::copyAs<T,S,util::StdVector>(control);
  util::iniDiagonal(&controlAD[0], numberOfDerivatives);
  auto resultAD = derivative(controlAD, 0);
  util::print(resultAD.d(), "derivatives", clout, ' ');
}

/** Optimization with automatic differentiation
 * "control" variable scales the force, objective measures distance of velocity
 * or dissipation to the exact solution
 */
void optiAd() {
  OstreamManager clout(std::cout, "optiAd");
  XMLreader config("parameterAd.xml");

  auto testFlow = createLbSolver <TestFlowSolverDirectOpti<S>> (config);

  auto testFlowAD = createLbSolver <TestFlowSolverDirectOpti<T>> (config);
  testFlow->parameters(Output()).verbose = false;
  testFlowAD->parameters(Output()).verbose = false;
  testFlowAD->parameters(VisualizationVTK()).output = "off";

  bool discreteObjective = testFlow->parameters(Opti()).optiReferenceMode;
  auto referenceSolver = createLbSolver <TestFlowSolverDirectOpti<S>> (config);
  std::vector<S> referenceControl (numberOfDerivatives, 1.0);
  if (discreteObjective) {
    // compute numerical solution for velocity/ dissipation for
    // correct control variables
    referenceSolver->parameters(Opti()).computeObjective = false;
    referenceSolver->parameters(Output()).verbose = false;
    referenceSolver->parameters(VisualizationVTK()).output = "off";

    referenceSolver->parameters(Opti()).applyControl(referenceControl);
    referenceSolver->solve();

    testFlow->parameters(Opti()).referenceState
     = referenceSolver->parameters(Results()).solution;

    testFlowAD->parameters(Opti()).referenceState
     = referenceSolver->parameters(Results()).solution;
  }

  OptiCaseAdForSolver optiCase(testFlow, testFlowAD);

  auto optimizer = createOptimizerLBFGS<S>(config, numberOfDerivatives);

  std::vector<S> startValues (numberOfDerivatives, 0.0);
  optimizer->setControl(startValues);
  optimizer->setReferenceControl(referenceControl);

  // Execute optimization
  optimizer->optimize(optiCase);

  clout << "objective = " << optimizer->getObjective() << std::endl;
  util::print(optimizer->getControl(),    "control   ", clout, ' ');
  util::print(optimizer->getDerivative(), "derivative", clout, ' ');
}

/** Optimization with adjoint LBM
 * "control" variable scales the force, objective measures distance of velocity
 * or dissipation to the exact solution
 */
void optiAdjoint() {
  OstreamManager clout(std::cout, "optiAdjoint");
  XMLreader config("parameterAdjoint.xml");

  OptiCaseDual<S,TestFlowSolverOptiAdjoint,descriptors::FORCE,ForcedBGKdynamics> optiCase (config);
  // classical approach: define objective via functors
  //auto objective = std::make_shared<TestFlowObjective<S>>(config);
  // use generic objective handling
  auto objectiveHelp = std::make_shared<RelativeDifferenceVelocityObjectiveGeneric<S>>(config);
  auto objective = std::make_shared<GenericObjective<S,TestFlowSolverOptiAdjoint,
    RelativeDifferenceVelocityObjectiveGeneric<S>,
    ForcedBGKdynamics, 3>>(objectiveHelp);
  optiCase.setObjective(objective);
  auto optimizer = createOptimizerLBFGS<S>(config, optiCase._dimCtrl);

  S startValue;
  config.readOrWarn<S>("Optimization", "StartValue", "", startValue);
  startValue = projection::getInitialControl(startValue, optiCase);
  const auto referenceControl = getControl<S,TestFlowSolverOptiAdjoint,descriptors::FORCE,ForcedBGKdynamics>(
    optiCase, objectiveHelp->_referenceSolver);

  optimizer->setReferenceControl(referenceControl);
  //optimizer->setStartValue(startValue);
  optimizer->setStartValue(referenceControl, startValue);  // interpret startValue as a scaling factor of referenceControl

  optimizer->optimize(optiCase);
}


int main(int argc, char **argv)
{
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  if (argc > 1) {
    for (int i = 1; i < argc; ++i) {
      if (strcmp(argv[i], "sensitivities") == 0) {
        clout << "\nCompute sensitivities with forward difference quotients" << std::endl;
        sensitivitiesFdq();

        clout << "\nCompute sensitivities with AD" << std::endl;
        sensitivitiesAd();

        clout << "\nCompute sensitivities with AD for 5 blocks" << std::endl;
        sensitivitiesAd_variant(5);
      } else if (strcmp(argv[i], "opti-ad") == 0) {
        clout << "\nOptimize with AD" << std::endl;
        optiAd();
      } else if (strcmp(argv[i], "opti-adjoint") == 0) {
        clout << "\nOptimize with adjoint lbm" << std::endl;
        optiAdjoint();
      } else {
        clout << "Usage: after program name, state some of the options sensitivities, opti-ad or opti-adjoint. "
              << "If no option is selected, all of them are performed.\n";
      }
    }
  } else {
    clout << "\nCompute sensitivities with forward difference quotients" << std::endl;
    sensitivitiesFdq();

    clout << "\nCompute sensitivities with AD" << std::endl;
    sensitivitiesAd();

    clout << "\nCompute sensitivities with AD for 5 blocks" << std::endl;
    sensitivitiesAd_variant(5);

    clout << "\nOptimize with AD" << std::endl;
    optiAd();

    clout << "\nOptimize with adjoint lbm" << std::endl;
    optiAdjoint();
  }
}

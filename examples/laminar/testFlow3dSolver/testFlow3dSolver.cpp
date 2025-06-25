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
 * Cf. testFlow3dSolver.h for implementation details
 */


#include "testFlow3dSolver.h"

using namespace olb;
using namespace olb::parameters;

using S = FLOATING_POINT_TYPE;

XMLreader config("parameter.xml");


/** Standard execution
 * Result: paraview output
 */
void simulate() {
  auto testFlow = createLbSolver <TestFlowSolver<S>> (config);
  testFlow->solve();
}

/** Experimental order of convergence computation
 * Results: i.e., an eoc plot is created via gnuplot
 */
void eoc() {
  Gnuplot<S> gplot(
    "eoc",
    false,
    "set terminal png size 720, 720 font 'Arial,10'",
    Gnuplot<S>::LOGLOG,
    Gnuplot<S>::LINREG);
  gplot.setLabel("resolution", "errors");

  for (unsigned N : {31, 41, 51, 61})  // spatial resolution
  {
    // set simulation parameters
    utilities::TypeIndexedSharedPtrTuple<Params_TfBasic<S>> params;
    params.template get<Errors>() = std::make_shared<SimulationErrors<S>>();
    params.template get<Output>() = createParameters<OutputGeneral<S>, Output>(config);
    params.template get<VisualizationVTK>() = createParameters<OutputPlot<S>, VisualizationVTK>(config);

    params.template get<Simulation>() = createParameters<TfSimulationParams<S,Lattices>, Simulation>(config);
    using descriptor = TestFlowSolver<S>::TestFlowBase::descriptor;
    auto& converter = params.template get<Simulation>()->converter;
    converter = std::make_shared<UnitConverterFromResolutionAndLatticeVelocity<S,descriptor>>(
      N,
      1. / N,  // viscosity scaling
      converter->getCharPhysLength(),
      converter->getCharPhysVelocity(),
      converter->getPhysViscosity(),
      converter->getPhysDensity(),
      converter->getCharPhysPressure()
    );

    // simulate
    TestFlowSolver<S> testFlow(params);
    testFlow.solve();

    // postprocessing: call gnuplot routines
    const auto& errors = testFlow.parameters(Errors());
    gplot.setData (
      S(N),
      { errors.velocityAbsL1Error, errors.velocityAbsL2Error, errors.velocityAbsLinfError,
        errors.pressureAbsL1Error, errors.pressureAbsL2Error, errors.pressureAbsLinfError,
        errors.strainRateAbsL1Error, errors.strainRateAbsL2Error, errors.strainRateAbsLinfError,
        errors.dissipationAbsL1Error, errors.dissipationAbsL2Error, errors.dissipationAbsLinfError },
      { "velocity L1 error","velocity L2 error",
        "velocity Linf error",
        "pressure L1 error", "pressure L2 error",
        "pressure Linf error",
        "strain rate L1 error", "strain rate L2 error",
        "strain rate Linf error",
        "dissipation L1 error", "dissipation L2 error",
        "dissipation Linf error" },
      "top right",
      { 'p','p','p','p','p','p','p','p','p','p','p','p' } );
  }
  gplot.writePNG();
}


int main(int argc, char **argv)
{
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  clout << "\nExecute standard simulation and compute error norms" << std::endl;
  simulate();

}

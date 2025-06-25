/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2006-2023 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod, Adrian Kummerländer
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

/* cylinder2d.cpp:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 */

// On this state it depends, wether simulation parameters can be changed
// or not. If the simulation is running or paused, the parameters are locked.
let simulationState = SIMULATION_STATE.STOPPED;
// let simulationValues = SIMULATION_VALUES.VELOCITY;
// The simulation step is the number of steps the simulation has already
// performed. It is used to calculate the number of steps the simulation
// should perform in the next frame. This value must be reset when the
// simulation is reset.
let simulationStep = 0;
let stepsPerRun = STEPS_PER_RUN_DEFAULT;
//
// Parameters for the simulation setup. These parameters are used to set up the simulation in the C++ code. They can
// only be changed when the simulation is stopped.
//
let resolution = RESOLUTION_DEFAULT; // resolution of the model
let reynolds = REYNOLDS_DEFAULT; // Reynolds number

const startSimulation = () => {
  Module._startTimer();
  runSimulation();
};

const pauseSimulation = () => {
  Module._stopTimer();
};

const stopSimulation = () => {
  Module._stopTimer();
  resetSimulation();
  //
  // So that the simulation can be started again after a stop without changed parameters.
  //
  setupSimulation();
};

const setupSimulation = () => {
  //
  // Order matters. setupCylinder2d() must be called after setSimulationParameters() and before getValuesPointer() etc.
  //
  Module._setSimulationParameters(reynolds, resolution);
  Module._setupCylinder2d();
};

const setSimulationOutput = (outputType) => {
  //
  // Set output type independently of the other parameters (reynolds, resolution) that need explicit set up.
  // C++ Simulation run just uses a different method to calculate the output values (switch case). Default is velocity.
  //
  Module._setOutputType(outputType);
};

const resetSimulation = () => {
  // Reset output array
  Module._freeMemory();

  // Reset the canvas
  clearCanvas();

  // Reset the simulation step
  simulationStep = 0;

  // Reset Text Output
  outputTimeEL.value = 0;
  outputPerformanceEL.value = 0;
};

const runSimulation = () => {
  //
  // Send the offscreen canvas to the worker thread only once. Otherwise, it will throw a detached error, because the
  // canvas does not get send back, thus we have no control in the main thread. Create new worker thread for every
  // simulation run to handle resolution change correctly.
  //
  const offscreen = new OffscreenCanvas(
    Module._getNumberCellsX(),
    Module._getNumberCellsY()
  );
  const canvasWorker = new Worker("./js/worker/canvas_worker.js");
  canvasWorker.postMessage({canvas: offscreen}, [offscreen]);

  canvasWorker.onmessage = (event) => {
    //
    // Condition since the worker thread hangs asynchronously behind. Thus, it would repaint the canvas when
    // simulation has been stopped directly.
    //
    if (simulationState === SIMULATION_STATE.RUNNING) {
      ctx.transferFromImageBitmap(event.data.bitmap);
    }
  };

  const animate = () => {
    if (simulationState === SIMULATION_STATE.RUNNING) {
      // Run the simulation for SIMULATION_STEPS_PER_RUN steps
      const i_start = simulationStep * stepsPerRun;
      const i_end = (simulationStep + 1) * stepsPerRun;
      Module._runCylinder2d(i_start, i_end);

      // Write physical time to output text field
      outputTimeEL.value = Module._getPhysTime()
        .toFixed(2)
        .split(".")
        .join(":");
      // Write MLUPs to output text field. Only every k-th step to reduce jitter.
      if (simulationStep % (stepsPerRun / 2) === 0) {
        outputPerformanceEL.value = Module._getMLUPs().toFixed(2);
      }

      //
      // Get the velocity or pressure values from the heap
      // TODO: maybe float 32? Intensive trial and error still did not show any results
      //
      const valuesArray = new Float64Array( // Wenn Float32Array, dann komisch
        // Module.HEAPF64.buffer,
        Module.HEAPF32.buffer, // still works even with FLoat64Array
        Module._getValuesPointer(),
        Module._getValuesCells()
      );

      //
      // Check for divergence
      //
      if (isNaN(Module._getAvgRho())) {
        simulationState = SIMULATION_STATE.STOPPED;
        setButtonState(SIMULATION_STATE.STOPPED);
        stopSimulation();
        alert("Divergence detected. Simulation stopped.");
        return;
      }

      //
      // Draw the values to the canvas in worker thread. ctx is the visible canvas context from canvas.js
      //
      canvasWorker.postMessage({
        values: valuesArray,
        colors: getGradientData(),
        max: Math.max(Module._getCharPhysVelocity(), Module._getPhysMaxU()),
      });

      simulationStep++;
      requestAnimationFrame(animate);
    }
  };

  animate();
};

const buttonPlayEl = document.getElementById("button-play");
const buttonPauseEl = document.getElementById("button-pause");
const buttonStopEl = document.getElementById("button-stop");
const buttonResetEl = document.getElementById("button-reset");
const sliderResolutionEl = document.getElementById("slider-n");
const sliderReynoldsEl = document.getElementById("slider-re");
const radioVelocityEL = document.getElementById("radio-velocity");
const radioPressureEL = document.getElementById("radio-pressure");

// Running
buttonPlayEl.addEventListener("click", () => {
  simulationState = SIMULATION_STATE.RUNNING;
  setButtonState(SIMULATION_STATE.RUNNING);
  startSimulation();
});

// Paused
buttonPauseEl.addEventListener("click", () => {
  simulationState = SIMULATION_STATE.PAUSED;
  setButtonState(SIMULATION_STATE.PAUSED);
  pauseSimulation();
});

// Stopped
buttonStopEl.addEventListener("click", () => {
  simulationState = SIMULATION_STATE.STOPPED;
  setButtonState(SIMULATION_STATE.STOPPED);
  stopSimulation();
});

// Reset Parameters
buttonResetEl.addEventListener("click", () => {
  sliderResolutionEl.value = RESOLUTION_DEFAULT;
  sliderReynoldsEl.value = REYNOLDS_DEFAULT;
  radioVelocityEL.checked = true;
  radioPressureEL.checked = false;
  resolution = RESOLUTION_DEFAULT;
  reynolds = REYNOLDS_DEFAULT;
  parametersChanged();
  setSimulationOutput(SIMULATION_VALUES.VELOCITY);
  setupSimulation();
});

// Parameters changed, thus disable play button until new parameters are saved and now simulation parameters can be reset
sliderResolutionEl.addEventListener("change", (e) => {
  resolution = e.target.value;
  parametersChanged();
  setupSimulation();
});

sliderReynoldsEl.addEventListener("change", (e) => {
  reynolds = e.target.value;
  parametersChanged();
  setupSimulation();
});

radioVelocityEL.addEventListener("change", () => {
  setSimulationOutput(SIMULATION_VALUES.VELOCITY);
});

radioPressureEL.addEventListener("change", () => {
  setSimulationOutput(SIMULATION_VALUES.PRESSURE);
});

const setButtonState = (state) => {
  buttonStopEl.disabled = state === SIMULATION_STATE.STOPPED;
  sliderReynoldsEl.disabled = state !== SIMULATION_STATE.STOPPED;
  sliderResolutionEl.disabled = state !== SIMULATION_STATE.STOPPED;
  buttonPlayEl.style.visibility =
    state === SIMULATION_STATE.RUNNING ? "hidden" : "visible";
  buttonPauseEl.style.visibility =
    state === SIMULATION_STATE.RUNNING ? "visible" : "hidden";
};

const parametersChanged = () => {
  buttonResetEl.disabled = !(
    resolution !== RESOLUTION_DEFAULT || reynolds !== REYNOLDS_DEFAULT
  );
  changeResolutionRange(reynolds);
};

const changeResolutionRange = (reynolds) => {
  if (RESOLUTION_RANGE[reynolds]) {
    sliderResolutionEl.min = RESOLUTION_RANGE[reynolds]["min_resolution"];
    sliderResolutionEl.max = RESOLUTION_RANGE[reynolds]["max_resolution"];
  } else {
    sliderResolutionEl.min = 15;
    sliderResolutionEl.max = 50;
  }
  resolution = Math.max(
    sliderResolutionEl.min,
    Math.min(resolution, sliderResolutionEl.max)
  );

  changeStepsPerRun(resolution);
};

const changeStepsPerRun = (resolution) => {
  stepsPerRun = STEPS_PER_RUN[resolution] || STEPS_PER_RUN_DEFAULT / 2;
};

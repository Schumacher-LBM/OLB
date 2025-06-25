const GRADIENTS = [
  "./gradients/gradient_wave_1.png",
  "./gradients/gradient_wave_2.png",
  "./gradients/gradient_wave_3.png",
  "./gradients/gradient_wave_4.png",
  "./gradients/gradient_wave_5.png",
  "./gradients/gradient_wave_6.png",
  "./gradients/gradient_wave_7.png",
];

const REYNOLDS_DEFAULT = 20;
const RESOLUTION_DEFAULT = 10;
const STEPS_PER_RUN_DEFAULT = 20;

const RESOLUTION_RANGE = {
  10: {
    min_resolution: 2,
    max_resolution: 13,
  },
  20: {
    min_resolution: 3,
    max_resolution: 13,
  },
  30: {
    min_resolution: 4,
    max_resolution: 15,
  },
  40: {
    min_resolution: 5,
    max_resolution: 17,
  },
  50: {
    min_resolution: 6,
    max_resolution: 19,
  },
  60: {
    min_resolution: 8,
    max_resolution: 19,
  },
  70: {
    min_resolution: 9,
    max_resolution: 20,
  },
  80: {
    min_resolution: 10,
    max_resolution: 20,
  },
  90: {
    min_resolution: 12,
    max_resolution: 20,
  },
  100: {
    min_resolution: 13,
    max_resolution: 20,
  },
};
const STEPS_PER_RUN = {
  1: 20,
  2: 20,
  3: 20,
  4: 20,
  5: 20,
  6: 20,
  7: 20,
  8: 20,
  9: 20,
  10: 20,
  11: 20,
  12: 20,
  13: 20,
  14: 10,
  15: 10,
  16: 10,
  17: 10,
  18: 10,
  19: 10,
  20: 10,
};

const SIMULATION_STATE = {
  RUNNING: "running",
  PAUSED: "paused",
  STOPPED: "stopped",
};
const SIMULATION_VALUES = {
  VELOCITY: true,
  PRESSURE: false,
};

const normalize = (value, minValue, maxValue) => {
  if (value > 1 || isNaN(value)) {
    value = 1;
  }
  return (value - minValue) / (maxValue - minValue);
};

const getColorFromGradient = (percentage, data) => {
  //
  // Factor 4 because of rgba. The gradient is stored in a 1D array, which is four times longer than the gradient has
  // values. Floor to round the index and -1 to accord for the array index shift.
  //
  const index = Math.floor(((data.length - 1) / 4) * percentage);
  const colors = data.slice(index * 4, index * 4 + 4);
  return {
    r: colors[0],
    g: colors[1],
    b: colors[2],
    a: colors[3],
  };
};

const measureExecutionTime = (fn) => {
  // Record the start time
  const startTime = performance.now();

  // Call the function
  fn();

  // Record the end time
  const endTime = performance.now();

  // Calculate the execution time in milliseconds
  const executionTime = endTime - startTime;

  console.log(
    `Function executed in ${executionTime} milliseconds or ${
      executionTime / 1000
    } seconds`
  );
};

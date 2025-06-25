importScripts("../utils.js");
//
// Canvas values get passed to the worker in the simulation run only once before the Simulation animates to save
// unnecessary function calls. Max is updated dynamically since physMaxU is dynamic.
//
let canvas;
let ctx;
let max;

const paint = (values, colors) => {
  const imageData = ctx.createImageData(canvas.width, canvas.height);
  const buf = new Uint8ClampedArray(imageData.data.buffer);

  for (let y = 0; y < canvas.height; y++) {
    for (let x = 0; x < canvas.width; x++) {
      //
      // Painting done this way, as ctx.fillRect() shows huge performance issues. (The reason is that ctx.fillRect() is
      // not hardware accelerated. The only way to paint the canvas is to use the ImageData object and write the
      // pixels directly into the buffer. This is the only way to get hardware acceleration.) Resulting in drawing
      // lagging behind the actual simulation.
      //
      const index = (y * canvas.width + x) * Uint32Array.BYTES_PER_ELEMENT; // 4 channels (RGBA) per pixel

      const color = getColorFromGradient(
        protectedNormalize(values[x + y * canvas.width]),
        colors
      );

      buf[index] = color.r;
      buf[index + 1] = color.g;
      buf[index + 2] = color.b;
      buf[index + 3] = 255; // Alpha channel, it's always 255
    }
  }

  ctx.putImageData(imageData, 0, 0);
};

const protectedNormalize = (value) => {
  //
  // The normalization value can be more than 1, because charPhysVelocity and phsMaxU are not true maximum values.
  //
  value = normalize(value, 0, max);
  if (value > 1 || isNaN(value)) {
    value = 1;
  }
  return value;
};

onmessage = (event) => {
  if (!event.data) return;

  if (!canvas) canvas = event.data.canvas;
  if (!ctx) ctx = canvas.getContext("2d");

  max = event.data.max;
  const values = event.data.values;
  const colors = event.data.colors;
  if (!values || !colors) return;

  paint(values, colors);

  self.postMessage({
    bitmap: canvas.transferToImageBitmap(),
  });
};

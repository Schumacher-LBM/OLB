// Create offscreen canvas for legend to retrieve color values from
const legendCanvas = new OffscreenCanvas(500, 1);
const legendCtx = legendCanvas.getContext("2d", {
  willReadFrequently: true,
});

const paintLegend = (gradient) =>
  legendCtx.drawImage(gradient, 0, 0, legendCanvas.width, legendCanvas.height);

const getGradientData = () =>
  legendCtx.getImageData(0, 0, legendCanvas.width, legendCanvas.height).data;

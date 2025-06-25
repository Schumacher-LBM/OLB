// Purpose: Canvas related functions
// bitmap render context since we paint to an offscreen canvas and then transfer the bitmap image to the visible canvas
const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("bitmaprenderer");

const clearCanvas = () => {
  ctx.transferFromImageBitmap(null);
};

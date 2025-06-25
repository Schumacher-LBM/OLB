const outputTimeEL = document.getElementById("output-time");
const outputPerformanceEL = document.getElementById("output-mlups");

//
// Needed to reset the output text fields when the page is reloaded. However, this is only a problem for soft reloads
// in Firefox.
//
document.addEventListener("DOMContentLoaded", () => {
  outputTimeEL.value = 0;
  outputPerformanceEL.value = 0;
});

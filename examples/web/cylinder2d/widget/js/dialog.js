const dialogAnchorEl = document.getElementById("dialog-anchor");
const dialogEL = document.getElementById("dialog");
const dialogButtonEL = document.getElementById("dialog-button");

// Dialog state
dialogAnchorEl.addEventListener(
  "click",
  () => (dialogEL.open = !dialogEL.open)
);

dialogButtonEL.addEventListener("click", () => dialogEL.close());

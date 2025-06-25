// This example uses anchor as an ID reference
const anchorEl = document.body.querySelector("#gradient-anchor");
const menuEl = document.body.querySelector("#gradient-menu");
const legendEL = document.body.querySelector("#legend");
let selectedGradientEl = null;

// Create Menu items
GRADIENTS.forEach((gradient, id) => {
  const menuItemEl = document.createElement("md-menu-item");
  menuItemEl.id = id.toString();
  menuItemEl.type = "option";

  const imgEl = document.createElement("img");
  imgEl.alt = `gradient-${id + 1}`;
  imgEl.className = "menu-img";
  imgEl.slot = "headline";
  imgEl.src = gradient;

  menuItemEl.appendChild(imgEl);
  menuEl.appendChild(menuItemEl);
});

// Menu state
anchorEl.addEventListener("click", () => {
  menuEl.open = !menuEl.open;
});

// Menu item selection
menuEl.addEventListener("close-menu", (e) => {
  // Select the item that was clicked
  const menuItemEl = e.target;
  // Overwrite the previous selected item. Set new selected last to account for the case that the same item was clicked.
  selectedGradientEl.selected = false;
  selectedGradientEl = menuItemEl;
  menuItemEl.selected = true;
  // Update the legend
  legendEL.src = menuItemEl.querySelector("img").src; // Only has one img element child
  paintLegend(legendEL);
});

// Set default selected item.
// Has to be done after the DOM is loaded otherwise the menuEl.items array is empty.
document.addEventListener("DOMContentLoaded", () => {
  menuEl.items[0].selected = true;
  selectedGradientEl = menuEl.items[0];
  paintLegend(legendEL);

  outputTimeEL.value = 0;
  outputPerformanceEL.value = 0;
});

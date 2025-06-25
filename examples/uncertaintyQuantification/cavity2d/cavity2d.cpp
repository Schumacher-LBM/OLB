/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mathias J. Krause, Jonas Fietz,
 *  Jonas Latt, Jonas Kratzke, Mingliang Zhong, Stephan Simonis
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

/* cavity2d.cpp:
 * Demonstrates a 2D lid-driven cavity flow where the lid velocity is uncertain.
 * Monte Carlo sampling and stochastic collocation are used to quantify how this
 * boundary uncertainty affects the flow field.
 */

#include "cavity2d.h"
#include <time.h>

// Define one method of sampling (uncomment if needed):
// #define MonteCarloSampling
#define stochasticCollocationMethod

/**
 * @brief UQ (Uncertainty Quantification) parameter setup:
 *
 * 1) Monte Carlo Sampling (Uncomment #define MonteCarloSampling):
 *    - Uses random draws (nq samples) from the specified distributions.
 *    - The number of samples nq can be set as an integer. Increasing nq yields
 *      more statistically robust results at higher computational cost.
 *    - Convergence Rate:
 *      - The Monte Carlo error typically decreases on the order of 1/sqrt(nq).
 *      - Consequently, to gain one extra digit of accuracy, the sample size nq must
 *        increase by about a factor of 100.
 *    - Example code block:
 *         int nq = resolution / 4.0;
 *         unsigned int seed = 123456; // fixed seed for reproducibility
 *         uq.initializeMonteCarlo(nq, dist, seed);
 *
 * 2) Stochastic Collocation (Default, #define stochasticCollocationMethod):
 *    - Implements generalized polynomial chaos (GPC).
 *    - Configure the polynomial expansion order (orderGPC) and the number
 *      of quadrature points per dimension (nqPerDim).
 *    - The total number of samples is (nqPerDim)^(number_of_dimensions).
 *    - Convergence Rate:
 *      - For smooth problems, polynomial-based methods like GPC can achieve
 *        exponential (spectral) convergence as the polynomial order is increased.
 *      - However, increasing orderGPC or nqPerDim grows the number of simulations
 *        rapidly, especially in high-dimensional parameter spaces.
 *    - Example code block:
 *         int orderGPC = 2;  // polynomial chaos order
 *         int nqPerDim = 3;  // quadrature points per dimension
 *         uq.initializeGPC(orderGPC, nqPerDim, dist);
 *
 * Adjust these parameters according to the desired accuracy and computational
 * resources. Higher order or more samples generally capture broader variations
 * but increase runtime.
 */

/**
  * @brief Helper function to create a directory (with '-p' for recursive creation).
  *        Logs a warning if directory creation fails.
  *
  * @param directory Path of the directory to create.
  * @param clout     Output stream for logging.
  */
void createDirectory(const std::string& directory, OstreamManager& clout) {
  std::string command = "mkdir -p " + directory;
  int result = std::system(command.c_str());
  if (result != 0) {
    clout << "Warning: Failed to create directory " << directory
          << " (return code = " << result << ")" << std::endl;
  } else {
    clout << "Directory created (or already exists): " << directory << std::endl;
  }
}

const T maxPhysT   = 100.0; //Maximum physical time for simulation (seconds)
const T vtkSave    = 1.0;   // Interval for writing VTK output (in physical seconds)

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  OstreamManager clout( std::cout,"main" );

  int resolution = 128;
  T physVelocity = 1.0;
  auto dist = uniform(0.9 * physVelocity, 1.1 * physVelocity);

  // Create main folder to store results
  std::string foldPath = "uq/res_" + std::to_string(resolution) + "/";
  createDirectory(foldPath, clout);

  #ifdef MonteCarloSampling
    int nq = resolution/4.0; // number of samples
    UncertaintyQuantification<T> uq(UQMethod::MonteCarlo);
    unsigned int seed = 123456; // fixed seed for reproducibility
    uq.initializeMonteCarlo(nq, dist, seed);
  #elif defined(stochasticCollocationMethod)
    UncertaintyQuantification<T> uq(UQMethod::GPC);
    int orderGPC = 2;
    int nqPerDim = 5; // nqPerDim is recommended to be 2*orderGPC+1
    uq.initializeGPC(orderGPC, nqPerDim, dist);
  #endif

  auto samples = uq.getSamplingPoints();

  Vector<T,2> extend( 1,1 );
  Vector<T,2> origin( 0,0 );
  IndicatorCuboid2D<T> cuboid( extend, origin );

#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, 1.0 / resolution, singleton::mpi().getSize() );
#else
  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, 1.0 / resolution, 1 );
#endif

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,2> superGeometry( cuboidDecomposition, loadBalancer );

  clout << "Starting simulation over " << samples.size() << " samples." << std::endl;

  for(int n = 0; n < samples.size(); ++n) {
    // Create subfolder for this sample
    std::string subFoldPath = foldPath + std::to_string(n) + "/tmp/";
    createDirectory(subFoldPath, clout);
    // Redirect output to this sample's folder
    singleton::directories().setOutputDir(subFoldPath);

    // Run the cavity2d simulation for this particular sample
    simulateCavity2d(samples[n][0], resolution, n);
  }

  // === 6. Post-processing: compute mean, std, and write VTI data ===
  // The string "physVelocity" is a tag used inside the function, not the variable name.
  computeMeanAndStdAndWriteVTI<T, DESCRIPTOR>(
    uq, foldPath, "cavity2d", "physVelocity", cuboidDecomposition, superGeometry
  );

}

/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mingliang Zhong, Stephan Simonis
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

/* tgv2d.cpp:
 * This example simulates the 2D Taylor-Green vortex problem, where flow
 * parameters (such as the Reynolds number) are uncertain. Monte Carlo sampling
 * and the stochastic collocation method are used to analyze how these uncertainties
 * influence the vortex evolution.
 *
 * Some note about the usage method, order, nq
 */

#include "tgv2d.h"
#include <time.h>

// Define one method of sampling (uncomment if needed):
#define MonteCarloSampling
// #define stochasticCollocationMethod

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
 *         int nq = resolution / 3.0;
 *         unsigned int seed = 123456; // fixed seed for reproducibility
 *         uq.initializeMonteCarlo(nq, jointPerturb, seed);
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
 *         uq.initializeGPC(orderGPC, nqPerDim, jointPerturb);
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
const T smagoConst = 0.1;   // Smagorinsky constant for LES

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  OstreamManager clout(std::cout, "main");

  // === 2nd Lattice / Physics Parameters ===
  int    resolution    = 64;
  if (argc > 1) {
    resolution = std::stoi(argv[1]);
  }
  T numVortices   = 2.0;
  T L             = 1.0;
  T Ma            = 1.6 / resolution;
  T Re            = 40.0 * resolution;

  bool   exportResults  = true;            // Whether to write VTI/VTM output

  // Characteristic length in TGV (Taylor–Green Vortex) setup
  T charL               = numVortices * M_PI * L;
  T dx                  = charL / resolution;
  T dt                  = dx * Ma / std::sqrt(3.0);
  T physViscosity       = 1.0 / Re;

  // === 3rd Uncertainty Quantification Setup ===
  // Define uniform distributions over given parameter ranges.
  auto perturb1 = uniform(-0.025, 0.025);
  auto perturb2 = uniform(-0.025, 0.025);
  auto perturb3 = uniform(-0.025, 0.025);
  auto perturb4 = uniform(-0.025, 0.025);

  // Combine these four separate uniform distributions into a single
  // joint distribution representing a 4-dimensional parameter space.
  auto jointPerturb = joint<T>({perturb1, perturb2, perturb3, perturb4});


  // Create main folder to store results
  std::string foldPath = "uq/res_" + std::to_string(resolution) + "/";
  createDirectory(foldPath, clout);

#ifdef MonteCarloSampling
  int nq = resolution/3.0; // number of samples
  UncertaintyQuantification<T> uq(UQMethod::MonteCarlo);
  unsigned int seed = 123456; // fixed seed for reproducibility
  uq.initializeMonteCarlo(nq, jointPerturb, seed);
#elif defined(stochasticCollocationMethod)
  UncertaintyQuantification<T> uq(UQMethod::GPC);
  int orderGPC = 2;
  int nqPerDim = 3; // jointPerturb is 4D, so nq overall will be nqPerDim^4
  uq.initializeGPC(orderGPC, nqPerDim, jointPerturb);
#endif

  // Retrieve the actual sample points from the UQ object
  auto samples = uq.getSamplingPoints();

  // === 4th Geometry Setup ===
  Vector<T,2> extend(charL, charL);
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid(extend, origin);

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif

  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, dx, noOfCuboids);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer);

  clout << "Starting simulation over " << samples.size() << " samples." << std::endl;

  // === 5th Loop over all samples and run simulations ===
  for (std::size_t n = 0; n < samples.size(); ++n) {
    if (exportResults) {
    // Create subfolder for this sample
      std::string subFoldPath = foldPath + std::to_string(n) + "/tmp/";
      createDirectory(subFoldPath, clout);
      // Redirect output to this sample's folder
      singleton::directories().setOutputDir(subFoldPath);
    }

    // Run the TGV simulation for this particular sample
    simulateTGV(physViscosity, dx, dt, samples[n], exportResults, n);
  }

  // === 6. Post-processing: compute mean, std, and write VTI data ===
  // The string "physVelocity" is a tag used inside the function, not the variable name.

  if (exportResults) {
    computeMeanAndStdAndWriteVTI<T, DESCRIPTOR>(
      uq, foldPath, "tgv2d", "physVelocity", cuboidDecomposition, superGeometry
    );
  }
  return 0;
}
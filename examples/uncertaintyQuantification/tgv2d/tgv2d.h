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
 * This example simulates the 2D Taylor-Green vortex problem, where certain flow
 * parameters (such as the Reynolds number) are uncertain. Monte Carlo sampling
 * and the stochastic collocation method are used to analyze how these uncertainties
 * influence the vortex evolution.
 */
#include <olb.h>

 using namespace olb;
 using namespace olb::uq;
 using namespace olb::descriptors;
 using namespace olb::graphics;

 using T = FLOATING_POINT_TYPE;
 using DESCRIPTOR = D2Q9<>;

 // Choose DNS or SMAGORINSKY LES model (uncomment as needed)
 #define DNS
 //#define _SMAGORINSKY

 #ifdef DNS
   using BulkDynamics = ConstRhoBGKdynamics<T,DESCRIPTOR>;
 #elif defined(_SMAGORINSKY)
   using BulkDynamics = SmagorinskyBGKdynamics<T,DESCRIPTOR>;
 #endif


 // -------------------------------------------------------------------------
 // Global parameters for the TGV simulation
 // -------------------------------------------------------------------------
 extern const T smagoConst;  ///< Smagorinsky constant for LES
 extern const T vtkSave;  ///< Interval for writing VTK output (in physical seconds)
 extern const T maxPhysT; ///< Maximum physical time for simulation (seconds)


 /**
  * @brief Base TGV solution in 2D:
  *        u_x = -u0 cos(x) sin(y)
  *        u_y =  u0 sin(x) cos(y)
  *
  * Currently not used in the code, but kept for reference.
  *
  * @tparam T           Floating-point type
  * @tparam _DESCRIPTOR Lattice descriptor
  */
 template <typename T, typename _DESCRIPTOR>
 class Tgv2D : public AnalyticalF2D<T,T> {
 protected:
   T u0; ///< Characteristic velocity

 public:
   Tgv2D(UnitConverter<T,_DESCRIPTOR> const& converter, T /*unusedFrac*/)
     : AnalyticalF2D<T,T>(2),
       u0(converter.getCharLatticeVelocity())
   {}

   bool operator()(T output[], const T input[]) override {
     const T x = input[0];
     const T y = input[1];

     output[0] = -u0 * util::cos(x) * util::sin(y);
     output[1] =  u0 * util::sin(x) * util::cos(y);
     return true;
   }
 };

 /**
  * @brief TGV solution in 2D with user-defined perturbations:
  *        Perturbation = \sum_i deltas[i] * [sin/cos(2x) * sin/cos(2y)]
  *
  * @tparam T           Floating-point type
  * @tparam _DESCRIPTOR Lattice descriptor
  */
 template <typename T, typename _DESCRIPTOR>
 class Tgv2D4U : public AnalyticalF2D<T,T> {
 protected:
   T u0;                         ///< Characteristic velocity
   std::vector<T> deltas_;  ///< Perturbation parameters

 public:
   Tgv2D4U(UnitConverter<T,_DESCRIPTOR> const& converter,
           const std::vector<T>& deltas)
     : AnalyticalF2D<T,T>(2),
       u0(converter.getCharLatticeVelocity()),
       deltas_(deltas)
   {}

   bool operator()(T output[], const T input[]) override {
     const T x = input[0];
     const T y = input[1];

     T u_perturbation = 0.0;
     u_perturbation += deltas_[0] * util::sin(2*x) * util::sin(2*y);
     u_perturbation += deltas_[1] * util::sin(2*x) * util::cos(2*y);
     u_perturbation += deltas_[2] * util::cos(2*x) * util::sin(2*y);
     u_perturbation += deltas_[3] * util::cos(2*x) * util::cos(2*y);

    const T factor = 1.0 + 0.25 * u_perturbation;
    // const T factor = 0.25 * u_perturbation;
     output[0] =  u0 * factor * util::sin(x) * util::cos(y);
     output[1] = -u0 * factor * util::cos(x) * util::sin(y);
     return true;
   }
 };


 /**
  * @brief Prepare the geometry by assigning material indicators and cleaning up boundaries.
  */
  void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
    SuperGeometry<T,2>& superGeometry) {
    OstreamManager clout(std::cout,"prepareGeometry");
    clout << "Prepare Geometry ..." << std::endl;

    // Assign material 1 to all cells
    superGeometry.rename(0, 1);
    superGeometry.communicate();

    // Remove unnecessary boundary voxels
    superGeometry.clean();
    superGeometry.innerClean();

    // Check for geometry errors
    superGeometry.checkForErrors();
    //  superGeometry.getStatistics().print();

    clout << "Prepare Geometry ... OK" << std::endl;
  }

 /**
  * @brief Define the bulk and boundary dynamics on the lattice.
  */
 void prepareLattice(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry) {
   OstreamManager clout(std::cout,"prepareLattice");
   clout << "Prepare Lattice ..." << std::endl;

   // Relaxation frequency in lattice units
   const T omega = converter.getLatticeRelaxationFrequency();

   // Material=0 -> No Dynamics
   sLattice.defineDynamics<NoDynamics<T,DESCRIPTOR>>(superGeometry, 0);

   // Material=1 -> Bulk Dynamics (DNS or Smagorinsky)
   sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);

   // Set BGK relaxation parameter
   sLattice.setParameter<descriptors::OMEGA>(omega);

   #ifdef _SMAGORINSKY
   // Set Smagorinsky constant (used in SmagorinskyBGKdynamics)
   sLattice.setParameter<collision::LES::SMAGORINSKY>(smagoConst);
   #endif

   clout << "Prepare Lattice ... OK" << std::endl;
 }

 /**
  * @brief Sets the initial density and velocity fields within the domain.
  *        Uses Tgv2D4U for velocity with random perturbations.
  */
 void setBoundaryValues( const UnitConverter<T,DESCRIPTOR>& converter,
                         SuperLattice<T, DESCRIPTOR>& sLattice,
                         SuperGeometry<T,2>& superGeometry,
                         const std::vector<double>& random) {
   OstreamManager clout(std::cout,"setBoundaryValues");

   // Constant density initialization
   AnalyticalConst2D<T,T> rhoF(1.0);

   // TGV velocity solution with perturbations
   Tgv2D4U<T,DESCRIPTOR> uSol(converter, random);

   // Set initial equilibrium in the bulk (material=1)
   auto bulkIndicator = superGeometry.getMaterialIndicator({1});
   sLattice.iniEquilibrium(bulkIndicator, rhoF, uSol);
   sLattice.defineRhoU(bulkIndicator, rhoF, uSol);

   // Prepare the lattice for simulation
   sLattice.initialize();
 }

 /**
  * @brief Periodically writes simulation results (velocity, pressure) to VTK.
  *        Collects iteration indices in a vector.
  */
 void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                  const UnitConverter<T,DESCRIPTOR>& converter,
                  std::size_t iT,
                  util::Timer<double>& timer,
                  std::vector<int>& iTList) {
   OstreamManager clout(std::cout,"getResults");

   // Set up a VTM writer for 2D data
   SuperVTMwriter2D<T> vtmWriter("tgv2d", 1, false);

   // Write geometry and rank data on the first iteration
   if (iT == 0) {
     SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(sLattice);
     SuperLatticeRank2D<T, DESCRIPTOR>   rank(sLattice);
     vtmWriter.write(cuboid);
     vtmWriter.write(rank);
     vtmWriter.createMasterFile();
   }

   // Write data every 'vtkSave' physical seconds
   if (iT % converter.getLatticeTime(vtkSave) == 0) {
     sLattice.setProcessingContext(ProcessingContext::Evaluation);

     // Add velocity, pressure functors
     SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
     SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, converter);
     vtmWriter.addFunctor(velocity);
     vtmWriter.addFunctor(pressure);
     vtmWriter.write(iT);

     // Update timer and log to console
     timer.update(iT);
     timer.printStep(2);
     sLattice.getStatistics().print(iT, converter.getPhysTime(iT));

     // Only the main MPI rank logs iteration steps
     int rank = 0;
     #ifdef PARALLEL_MODE_MPI
       rank = singleton::mpi().getRank();
     #endif
     if (rank == 0) {
       iTList.push_back(static_cast<int>(iT));
     }
   }
 }


 void save_velocity_field(int nx, int ny, std::vector<std::vector<T>> u, std::vector<std::vector<T>> v, int iter) {
  std::string filename_u = singleton::directories().getLogOutDir() + "/u_" + std::to_string(iter) + ".dat";
  std::ofstream outputFile_u(filename_u);
  if (!outputFile_u) {
    std::cerr << "Error opening the file: " << filename_u << std::endl;
    return;
  }
  outputFile_u << std::fixed << std::setprecision(15);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      outputFile_u << u[i][j] << "\t";
    }
    outputFile_u << "\n";
  }
  outputFile_u.close();

  std::string filename_v = singleton::directories().getLogOutDir() + "/v_" + std::to_string(iter) + ".dat";
  std::ofstream outputFile_v(filename_v);
  if (!outputFile_v) {
    std::cerr << "Error opening the file: " << filename_v << std::endl;
    return;
  }
  outputFile_v << std::fixed << std::setprecision(15);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      outputFile_v << v[i][j] << "\t";
    }
    outputFile_v << "\n";
  }

  outputFile_v.close();

}

void saveResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                  UnitConverter<T,DESCRIPTOR> const& converter,
                  SuperGeometry<T,2>& superGeometry,
                  int resolution,
                  int idx,
                  int iT ) {

  OstreamManager clout( std::cout,"output" );
  const int UQIter = converter.getLatticeTime( maxPhysT * 0.2 );

  if (iT % UQIter == 0) {
    std::vector<std::vector<T>> u(resolution, std::vector<T>(resolution));
    std::vector<std::vector<T>> v(resolution, std::vector<T>(resolution));

    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocityField( sLattice, converter );
    AnalyticalFfromSuperF2D<T> intpolateVelocity( velocityField, true, 1 );

    T dx = converter.getPhysDeltaX();

    for (int i = 0; i < resolution; ++i) {
      for (int j = 0; j < resolution; ++j) {
        T position[2] = {i * dx, j * dx};
        T velocity[2];
        intpolateVelocity(velocity, position);
        u[i][j] = velocity[0];
        v[i][j] = velocity[1];
      }
    }

    save_velocity_field( resolution, resolution, u, v, iT / UQIter);
    clout << "saved " << iT / UQIter << "\tphysT: " << converter.getPhysTime(iT) << std::endl;
  }

}


 /**
 * @brief Main function for running the Taylor–Green Vortex (TGV) simulation
 *        with optional perturbations in the velocity field.
 *
 * @param physViscosity  Physical viscosity in m^2/s
 * @param dx             Lattice spacing
 * @param dt             Lattice time step
 * @param random         Vector of perturbation parameters
 * @param idx            Sample index for labeling/logging
 */
 void simulateTGV( T physViscosity,
                   T dx,
                   T dt,
                   const std::vector<T>& random,
                   bool exportResults,
                   int idx) {
   OstreamManager clout(std::cout,"simulateTGV");
   clout << "Start sample " << idx << std::endl;

   // Create a UnitConverter based on given parameters
   UnitConverter<T,DESCRIPTOR> converter(
           dx,
           dt,
     (T)   1.0,                // charPhysLength
     (T)   1.0,                // charPhysVelocity
     (T)   physViscosity,
     (T)   1.0                 // physDensity
   );

   // Print converter info on screen and save to file
   converter.print();
   converter.write("tgv2d");

 #ifdef PARALLEL_MODE_MPI
   const int noOfCuboids = singleton::mpi().getSize();
 #else
   const int noOfCuboids = 4;
 #endif

   // Define the 2D domain from (0,0) to (2π,2π)
   Vector<T,2> extend(2 * M_PI, 2 * M_PI);
   Vector<T,2> origin;
   IndicatorCuboid2D<T> cuboid(extend, origin);

   // Decompose the domain and set periodic boundaries in x, y
   CuboidDecomposition2D<T> cuboidDecomposition(cuboid, dx, noOfCuboids);
   cuboidDecomposition.setPeriodicity({true, true});

   // Build a super geometry with load balancing
   HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
   SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer);

   // Prepare geometry and lattice data structures
   prepareGeometry(converter, superGeometry);
   SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);
   prepareLattice(converter, sLattice, superGeometry);

   // Create timer for performance measurements
   util::Timer<T> timer(
     converter.getLatticeTime(maxPhysT),
     superGeometry.getStatistics().getNvoxel()
   );

   // Optional convergence check
   int interval = converter.getLatticeTime(1);
   T epsilon = static_cast<T>(1e-3);
   util::ValueTracer<T> converge(interval, epsilon);

   // Set initial fields (rho, velocity)
   setBoundaryValues(converter, sLattice, superGeometry, random);

   timer.start();
   std::size_t iT = 0;
   std::vector<int> iTList;

   // Main time loop
   while (iT <= converter.getLatticeTime(maxPhysT)) {
     sLattice.collideAndStream();                  // Collide & stream
     if (exportResults) {
      getResults(sLattice, converter, iT, timer, iTList); // Output results
     }
     saveResults( sLattice, converter, superGeometry, (converter.getReynoldsNumber() / 40)+1, idx, iT );
     ++iT;
   }

   // Write iteration indices to a file (only on the main processor)
   if (singleton::mpi().isMainProcessor()) {
     std::ofstream outFile(singleton::directories().getLogOutDir() + "iteration_log.txt");
     if (outFile.is_open()) {
       for (auto step : iTList) {
         outFile << step << "\n";
       }
       outFile.close();
     } else {
       clout << "Error: Unable to open iteration_log.txt for writing!" << std::endl;
     }
   }

   timer.stop();
   timer.printSummary();
 }
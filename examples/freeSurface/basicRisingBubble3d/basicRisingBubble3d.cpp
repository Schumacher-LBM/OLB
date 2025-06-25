/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024-2025 Danial Khazaeipoul
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

/*
 * This example presents dynamics of a single bubble rising in a column of stagnant liquid.
 * The bubble is initialized at the bottom-center of the domain and starts to rise due to the
 * buoyancy force. The free surface Lattice Boltzmann Method (FSLBM) is used to simulate the
 * bubble dynamics. It is important to note that this is a basic version and is only used for
 * demonstration purpose. The actual case must implement a Bubble Model to track and update
 * bubble density for a more realistic simulation.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

//- Define the floating point type: single- or double-precision
using T = FLOATING_POINT_TYPE;

//- Define the lattice descriptor and bulk dynamics for the Navier-Stokes equations
using NavierStokesDescriptor = D3Q27
<
  descriptors::FORCE,
  FreeSurface::MASS,
  FreeSurface::EPSILON,
  FreeSurface::CELL_TYPE,
  FreeSurface::CELL_FLAGS,
  FreeSurface::TEMP_MASS_EXCHANGE,
  FreeSurface::PREVIOUS_VELOCITY
>;
using NavierStokesBulkDynamics = SmagorinskyForcedBGKdynamics<T, NavierStokesDescriptor>;

//- Set the fields based on the initial volume fraction: a spherical gas bubble in a liquid
template<typename T>
class bubble final : public AnalyticalF3D<T, T>
{
protected:

  //- Lattice grid spacing in SI units [m]
  T deltaX_;

  //- Bubble radius in SI units [m]
  T radius_;

  //- Bubble center in SI units [m]
  std::vector<T> center_;

  //- Sets field values, i.e., independent of the descriptor dimension; always uses 3 components.
  std::array<T, 3> values_;

public:

  //- Constructor from the input file, Navier-Stokes unit converter, and field values
  bubble
  (
    XMLreader& configuration, UnitConverter<T, NavierStokesDescriptor> const& NavierStokesConverter,
    const std::array<T,3>& values
  ) : AnalyticalF3D<T, T>(1), deltaX_(NavierStokesConverter.getPhysDeltaX()), radius_(T(0)), values_(values)
  {
    //- Reading bubble radius [m] and its center [m] from the input file
    radius_ = configuration["Geometry"]["Bubble"]["Radius"].get<T>();
    configuration["Geometry"]["Bubble"]["Center"].read(center_);
  }

  //- Set the initial volume fraction field based on the bubble geometry, i.e., gas (0), interface (1), fluid (2)
  bool operator()(T output[], const T x[]) override {
    //- Initially, set the current node as fluid
    output[0] = values_[2];

    //- Calculate the distance from the bubble center and verify if the node is within the bubble.
    //- Additionally, check if the node lies within the specified limit to be considered as interface.
    Vector<T, 3> fromCtr = {x[0] - center_[0], x[1] - center_[1], x[2] - center_[2]};

    if (norm(fromCtr) <= radius_)
    {
      output[0] = values_[0];
    }
    else
    {
      const T cutoff = deltaX_ * T{1.1};
      bool found = false;

      util::forEachOffset<3>([&](int i, int j, int k) {
        if (!found) {
          Vector<T, 3> offset = {fromCtr[0] + i * cutoff, fromCtr[1] + j * cutoff, fromCtr[2] + k * cutoff};

          if (norm(offset) <= radius_) {
            output[0] = values_[1];
            found = true;
          }
        }
      });
    }

    return true;
  }
};

//- Initializing hydrostatic pressure within the domain
template<typename T>
class HydroStaticPressure3D final : public AnalyticalF3D<T, T>
{
protected:

  //- Navier-Stokes unit converter
  UnitConverter<T, NavierStokesDescriptor> const& converter_;

  //- Liquid column width in z-direction [m]
  T width_;

  //- The gravity vector [m/s^2]
  std::vector<T> gravity_;

public:

  //- Constructor from the Advection-Diffusion unit converter
  HydroStaticPressure3D
  (
    XMLreader& configuration, UnitConverter<T, NavierStokesDescriptor> const& NavierStokesConverter
  ) : AnalyticalF3D<T, T>(1), converter_(NavierStokesConverter)
  {
    //- Reading the liquid column width [m] and gravity [m/s^2] from the input file
    width_ = configuration["Geometry"]["Width"].get<T>();
    configuration["Application"]["PhysParameters"]["PhysGravity"].read(gravity_);
  }

  //- Density computed using hydrostatic pressure (p = p0 + rho * g * h),
  //  with the reference density established at the midpoint in the z-direction
  bool operator()(T output[], const T x[]) override {
    output[0] = converter_.getLatticeDensityFromPhysPressure
    (
      converter_.getCharPhysPressure() + converter_.getPhysDensity() * util::abs(gravity_[2]) * (T(0.5) * width_ - x[2])
    );

    return true;
  }
};

//- Calculate the center of mass of the bubble, i.e., works only for a single bubble
Vector<T, 3> bubbleCenterOfMass
(
    SuperGeometry<T, 3>& superGeometry,
    SuperLattice<T, NavierStokesDescriptor>& NavierStokesLattice
)
{
  //- Initialize the center of mass of the bubble
  Vector<T, 3> centerOfMass(T(0));
  T nCells = T(0);

  for (int iC = 0; iC < NavierStokesLattice.getLoadBalancer().size(); ++iC)
  {
    const auto& block = NavierStokesLattice.getBlock(iC);
    block.forCoreSpatialLocations([&](auto iX, auto iY, auto iZ) {
      auto cell = block.get({iX, iY, iZ});

      //- Calculate the center of mass of the bubble, i.e., is Gas Type (0)
      if (isCellType(cell, FreeSurface::Type::Gas))
      {
        centerOfMass += superGeometry.getBlockGeometry(iC).getPhysR({iX, iY, iZ});
        nCells += T(1);
      }
    });
  }

#ifdef PARALLEL_MODE_MPI
  for (size_t i = 0; i < centerOfMass.size(); ++i)
  {
    singleton::mpi().reduceAndBcast(centerOfMass[i], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(nCells, MPI_SUM);
#endif

  centerOfMass[0] /= nCells;
  centerOfMass[1] /= nCells;
  centerOfMass[2] /= nCells;

  return centerOfMass;
}

//- Prepare the geometry using the indicator method
void prepareGeometry
(
  XMLreader& configuration, SuperGeometry<T, 3>& superGeometry,
  UnitConverter<T, NavierStokesDescriptor> const& NavierStokesConverter
)
{
  OstreamManager clout(std::cout, "prepareGeometry");

  //- Setting the material ID to (1) for all voxels, i.e. domain-wide
  //- Note that the side walls in x- and y-directions are periodic
  superGeometry.rename(0, 1);

  //- Reading cuboid length, height, and width [m] from the input file
  const T length = configuration["Geometry"]["Length"].get<T>();
  const T height = configuration["Geometry"]["Height"].get<T>();
  const T width = configuration["Geometry"]["Width"].get<T>();
  const T epsilon = T(0.5) * NavierStokesConverter.getPhysDeltaX();
  std::vector<T> origin = {-epsilon, -epsilon, -epsilon};
  std::vector<T> extend = {length + T(2) * epsilon, height + T(2) * epsilon, T(2) * epsilon};

  //- Setting material ID 2 for the bottom wall boundary voxels
  IndicatorCuboid3D<T> bottomWall(extend, origin);
  superGeometry.rename(1, 2, bottomWall);

  //- Setting material ID 3 for the top wall boundary voxels
  origin[2] += width;
  IndicatorCuboid3D<T> topWall(extend, origin);
  superGeometry.rename(1, 3, topWall);

  //- Removes all not needed boundary voxels outside the surface
  //- Removes all not needed boundary voxels inside the surface
  //- Checks the geometry for possible errors
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Geometry created successfully." << std::endl;
}

//- Prepare lattice for the Navier-Stokes equations with the free-surface model
void prepareNavierStokesLattice
(
  XMLreader& configuration, SuperGeometry<T, 3>& superGeometry,
  UnitConverter<T, NavierStokesDescriptor> const& NavierStokesConverter,
  SuperLattice<T, NavierStokesDescriptor>& NavierStokesLattice
)
{
  OstreamManager clout(std::cout, "prepareNavierStokesLattice");

  //- Lattice relaxation frequency for the Navier-Stokes equations
  T NavierStokesOmega  = NavierStokesConverter.getLatticeRelaxationFrequency();

  //- rho0 is the initial density profile calculated from the hydro-static pressure
  HydroStaticPressure3D<T> rho0(configuration, NavierStokesConverter);

  //- velocity0 is the uniform initial velocity profile
  AnalyticalConst3D<T, T> velocity0(T(0), T(0), T(0));

  //- Setting NSEs bulk dynamics for the defined material IDs
  //- Bulk materials: (1) fluid, (2) bottomWall, and (3) topWall
  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3});
  NavierStokesLattice.defineDynamics<NoDynamics>(superGeometry, 0);
  NavierStokesLattice.defineDynamics<NavierStokesBulkDynamics>(bulkIndicator);

  //- Setting NSEs velocity boundary condition for the material (2), i.e. bottomWall
  boundary::set<boundary::InterpolatedVelocity>(NavierStokesLattice, superGeometry, 2);

  //- Setting NSEs velocity boundary condition for the material (3), i.e. topWall
  boundary::set<boundary::InterpolatedVelocity>(NavierStokesLattice, superGeometry, 3);

  //- Setting the lattice relaxation frequency and Smagorinsky constant
  NavierStokesLattice.setParameter<descriptors::OMEGA>(NavierStokesOmega);
  NavierStokesLattice.setParameter<collision::LES::SMAGORINSKY>(T(0.2));

  //- Define user-defined constants for field initialization
  AnalyticalConst3D<T, T> zero(T(0));
  AnalyticalConst3D<T, T> one(T(1));
  AnalyticalConst3D<T, T> two(T(2));
  AnalyticalConst3D<T, T> four(T(4));
  AnalyticalConst3D<T, T> zeros(T(0), T(0), T(0));

  //- Create a single bubble with the specified center and radius
  bubble<T> types(configuration, NavierStokesConverter, {T(0), T(1), T(2)});
  bubble<T> mass(configuration, NavierStokesConverter, {T(0), T(0.5), T(1)});

  //- Initialize the free surface fields on the lattice
  for (auto i : {0, 1, 2, 3}) {
    NavierStokesLattice.defineField<FreeSurface::MASS>(superGeometry, i, zero);
    NavierStokesLattice.defineField<FreeSurface::EPSILON>(superGeometry, i, zero);
    NavierStokesLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, i, zero);
    NavierStokesLattice.defineField<FreeSurface::CELL_FLAGS>(superGeometry, i, zero);
    NavierStokesLattice.defineField<FreeSurface::PREVIOUS_VELOCITY>(superGeometry, i, zeros);
    NavierStokesLattice.defineField<FORCE>(superGeometry, i, zeros);
  }

  //- Setting the initial types, volume fraction, and mass based on the bubble geometry for the material IDs
  //- Material ID 1 denotes the internal fluid nodes
  NavierStokesLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, 1, types);
  NavierStokesLattice.defineField<FreeSurface::MASS>(superGeometry, 1, mass);
  NavierStokesLattice.defineField<FreeSurface::EPSILON>(superGeometry, 1, mass);

  //- Materials ID 2 and 3 denote the bottom and top wall boundary voxels
  //- Note that epsilon is set to 1.0 for the wall boundary voxels, i.e., a wetting wall condition
  for (auto i : {2, 3}) {
    NavierStokesLattice.defineField<FreeSurface::CELL_TYPE>(superGeometry, i, four);
    NavierStokesLattice.defineField<FreeSurface::EPSILON>(superGeometry, i, one);
  }

  //- Setting the gravity force in the z-direction, i.e., in lattice units
  std::vector<T> gravity;
  configuration["Application"]["PhysParameters"]["PhysGravity"].read(gravity);
  T gScale = T(1.0) / NavierStokesConverter.getConversionFactorForce() * NavierStokesConverter.getConversionFactorMass();
  AnalyticalConst3D<T,T> gForce{gravity[0] * gScale, gravity[1] * gScale, gravity[2] * gScale};
  NavierStokesLattice.defineField<FORCE>(superGeometry, 1, gForce);

  //- Setting the initial conditions for the Navier-Stokes equations
  NavierStokesLattice.defineRhoU(bulkIndicator, rho0, velocity0);
  NavierStokesLattice.iniEquilibrium(bulkIndicator, rho0, velocity0);

  FreeSurface::initialize(NavierStokesLattice);
  NavierStokesLattice.initialize();

  clout << "Lattice created successfully." << std::endl;
}

//- Reading and setting the free-surface parameters from the input file
void setFreeSurfaceParameters
(
  XMLreader& configuration, UnitConverter<T, NavierStokesDescriptor> const& NavierStokesConverter,
  SuperLattice<T, NavierStokesDescriptor>& NavierStokesLattice
)
{
  //- Setting the drop isolated cells switch
  bool dropIsolatedCells = false;
  configuration["Application"]["PhysParameters"]["DropIsolatedCells"].read<bool>(dropIsolatedCells);

  //- Setting the anti-jitter and lonely cells removal threshold
  T transitionThreshold = configuration["Application"]["PhysParameters"]["TransitionThreshold"].get<T>();
  T lonelyThreshold = configuration["Application"]["PhysParameters"]["LonelyThreshold"].get<T>();

  //- Setting the surface tension coefficient value
  bool hasSurfaceTension = false;
  configuration["Application"]["PhysParameters"]["HasSurfaceTension"].read<bool>(hasSurfaceTension);
  T surfaceTensionCoeff = configuration["Application"]["PhysParameters"]["SurfaceTensionCoeff"].get<T>();
  T tScale = util::pow(NavierStokesConverter.getConversionFactorTime(), 2)
    / (NavierStokesConverter.getPhysDensity() * util::pow(NavierStokesConverter.getPhysDeltaX(), 3));

  //- Calculate the force conversion factor
  // T gScale = T(1.0) / NavierStokesConverter.getConversionFactorForce() * NavierStokesConverter.getConversionFactorMass();

  //- Setting the free-surface parameters
  NavierStokesLattice.setParameter<FreeSurface::DROP_ISOLATED_CELLS>(dropIsolatedCells);
  NavierStokesLattice.setParameter<FreeSurface::TRANSITION>(transitionThreshold);
  NavierStokesLattice.setParameter<FreeSurface::LONELY_THRESHOLD>(lonelyThreshold);
  NavierStokesLattice.setParameter<FreeSurface::HAS_SURFACE_TENSION>(hasSurfaceTension);
  NavierStokesLattice.setParameter<FreeSurface::SURFACE_TENSION_PARAMETER>(tScale * surfaceTensionCoeff);
}

//- Set the output data for post-processing and visualization
void postprocess
(
  XMLreader& configuration, const std::size_t& iT, util::Timer<T>& timer,
  SuperGeometry<T, 3>& superGeometry,
  UnitConverter<T, NavierStokesDescriptor> const& NavierStokesConverter,
  SuperLattice<T, NavierStokesDescriptor>& NavierStokesLattice,
  std::size_t& iTPrev, Vector<T, 3>& CenterOfMassPrev
)
{
  OstreamManager clout(std::cout, "postprocess");

  //- Setting the output file name and save time for visualization
  std::string VTKFileName;
  configuration["Output"]["VisualizationVTK"]["FileName"].read(VTKFileName);
  SuperVTMwriter3D<T> vtmWriter(VTKFileName);
  std::string CSVFileName;
  configuration["Output"]["SampleData"]["FileName"].read(CSVFileName);
  static CSV<T> csvWriter(CSVFileName, ',', {"Time", "X", "Y", "Z", "RiseVelocity"});
  const T vtkSaveTime = configuration["Output"]["VisualizationVTK"]["SaveTime"].get<T>();
  const T logTime = configuration["Output"]["VisualizationVTK"]["SaveTime"].get<T>();
  const T maxPhysTime = configuration["Application"]["PhysParameters"]["MaxPhysTime"].get<T>();
  const bool lastTimeStep = (iT + 1 == NavierStokesConverter.getLatticeTime(maxPhysTime));

  if (iT == 0) {
    SuperLatticeCuboid3D<T, NavierStokesDescriptor> cuboid(NavierStokesLattice);
    SuperLatticeRank3D<T, NavierStokesDescriptor> rank(NavierStokesLattice);

    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();

    //- Reading the material properties from the input file in SI units
    const T physDensity = NavierStokesConverter.getPhysDensity();
    const T physViscosity = NavierStokesConverter.getPhysViscosity();
    const T surfaceTensionCoeff = configuration["Application"]["PhysParameters"]["SurfaceTensionCoeff"].get<T>();
    std::vector<T> gravity;
    configuration["Application"]["PhysParameters"]["PhysGravity"].read(gravity);
    const T radius = configuration["Geometry"]["Bubble"]["Radius"].get<T>();

    //- Calculate the Bond number for the rising bubble
    const T Bond = util::abs(T(4) * physDensity * gravity[2] * util::pow(radius, 2) / surfaceTensionCoeff);
    clout << "Bond number = " << Bond << std::endl;

    //- Calculate the Morton number for the rising bubble
    const T Morton = util::abs(gravity[2] * util::pow(physViscosity, 4) * util::pow(physDensity / surfaceTensionCoeff, 3));
    clout << "Morton number = " << Morton << std::endl;
  }

  //- Write the output data for visualization in VTK format
  if (iT % NavierStokesConverter.getLatticeTime(vtkSaveTime) == 0 || lastTimeStep) {
    //- Triggers data transfer between the host (CPU) and the device (GPU)
    NavierStokesLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity3D<T, NavierStokesDescriptor> velocity(NavierStokesLattice, NavierStokesConverter);
    SuperLatticePhysPressure3D<T, NavierStokesDescriptor> pressure(NavierStokesLattice, NavierStokesConverter);
    SuperLatticeDensity3D<T, NavierStokesDescriptor> density(NavierStokesLattice);
    SuperLatticeExternalScalarField3D<T, NavierStokesDescriptor, FreeSurface::CELL_TYPE> types(NavierStokesLattice);
    SuperLatticeExternalScalarField3D<T, NavierStokesDescriptor, FreeSurface::EPSILON> epsilon(NavierStokesLattice);
    SuperLatticeExternalScalarField3D<T, NavierStokesDescriptor, FreeSurface::MASS> mass(NavierStokesLattice);

    vtmWriter.addFunctor(velocity, "Velocity");
    vtmWriter.addFunctor(pressure, "Pressure");
    vtmWriter.addFunctor(density, "Density");
    vtmWriter.addFunctor(types, "Types");
    vtmWriter.addFunctor(epsilon, "Epsilon");
    vtmWriter.addFunctor(mass, "Mass");
    vtmWriter.write(iT);

    //- Write bubble statistics in SI units: id, center of mass [m], and rise velocity [m/s]
    Vector<T, 3> centerOfMass = bubbleCenterOfMass(superGeometry, NavierStokesLattice);
    const T time = NavierStokesConverter.getPhysTime(iT);
    const T timePrev = NavierStokesConverter.getPhysTime(iTPrev);
    const T riseVelocity = (centerOfMass[2] - CenterOfMassPrev[2]) / (time - timePrev + T(1e-14));
    csvWriter.writeDataFile(time, {centerOfMass[0], centerOfMass[1], centerOfMass[2], riseVelocity});
    iTPrev = iT;
    CenterOfMassPrev = centerOfMass;
    const T Reynolds = riseVelocity * NavierStokesConverter.getCharPhysLength() / NavierStokesConverter.getPhysViscosity();
    clout << "Calculated Re = " << Reynolds << ", Experimental Re = 94" << std::endl;
    clout << "Bubble's center of mass = (" << centerOfMass[0] << ", " << centerOfMass[1] << ", " << centerOfMass[2] << ")" << std::endl;
  }

  //- Write the output data for logging in the console
  if (iT % NavierStokesConverter.getLatticeTime(logTime) == 0 || lastTimeStep) {
    timer.update(iT);
    timer.printStep();
    NavierStokesLattice.getStatistics().print(iT, NavierStokesConverter.getPhysTime(iT));
  }
}

int main(int argc, char* argv[])
{
  //-- 1st Step: Reading input files and solver initialization
  OstreamManager clout(std::cout, "main");
  initialize(&argc, &argv, false, false);

  std::string inputFileName("input.xml");
  XMLreader configuration(inputFileName);

  //- Reading and setting olb source directory
  std::string olbPath;
  configuration["Application"]["OlbDir"].read(olbPath);
  singleton::directories().setOlbDir(olbPath);

  //- Reading and setting the output directory
  std::string outputPath;
  configuration["Output"]["OutputDir"].read(outputPath);
  singleton::directories().setOutputDir(outputPath);

  //- Reading the maximum simulation time [s]
  const T maxPhysTime = configuration["Application"]["PhysParameters"]["MaxPhysTime"].get<T>();

  //- Reading cuboid length, height, and width [m] from the input file
  const T length = configuration["Geometry"]["Length"].get<T>();
  const T height = configuration["Geometry"]["Height"].get<T>();
  const T width = configuration["Geometry"]["Width"].get<T>();

  //- Setting up lattice unit converters using physical parameters in SI units
  //- The Navier- Stokes (NSE) unit converter read the input parameters automatically from
  //- the configuration file, while a few parameters must be read manually
  UnitConverter<T, NavierStokesDescriptor>* NavierStokesConverterPtr =
    createUnitConverter<T, NavierStokesDescriptor>(configuration);

  UnitConverter<T, NavierStokesDescriptor> const& NavierStokesConverter = *NavierStokesConverterPtr;
  NavierStokesConverter.print();
  NavierStokesConverter.write("NavierStokesConverter");

  //-- 2nd Step: Prepare the geometry using indicator method
  const std::vector<T> origin(3, T(0));
  const std::vector<T> extend = {length, height, width};
  IndicatorCuboid3D<T> cube(extend, origin);

  //- Construct a 3D cuboid geometry, with weights and load balancing
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidDecomposition3D<T> cuboidDecomposition(cube, NavierStokesConverter.getPhysDeltaX(), noOfCuboids);

  //- Setting periodic condition in the x and y-directions
  cuboidDecomposition.setPeriodicity({true, true, false});

  //- Construct a load balancer by assigning cuboids to threads or MPI processes
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  //- Construct a load-balanced super geometry from a 3d cuboidGeometry
  SuperGeometry<T, 3> superGeometry(cuboidDecomposition, loadBalancer, 2);

  //- Prepare the geometry by assigning material IDs to individual voxels
  prepareGeometry(configuration, superGeometry, NavierStokesConverter);

  //-- 3rd Step: Prepare Navier-Stokes lattice with the free-surface model
  SuperLattice<T, NavierStokesDescriptor> NavierStokesLattice(superGeometry);
  prepareNavierStokesLattice(configuration, superGeometry, NavierStokesConverter, NavierStokesLattice);
  FreeSurface3DSetup<T, NavierStokesDescriptor> freeSurfaceHandler(NavierStokesLattice);
  freeSurfaceHandler.addPostProcessor();
  setFreeSurfaceParameters(configuration, NavierStokesConverter, NavierStokesLattice);

  //-- 4th Step: Prepare the main time loop using Timer and convergence check
  util::Timer<T> timer(NavierStokesConverter.getLatticeTime(maxPhysTime), superGeometry.getStatistics().getNvoxel());
  timer.start();

  //- Enable statistics by tracking bubble's center of mass
  std::size_t iTPrev = 0;
  Vector<T, 3> CenterOfMassPrev = bubbleCenterOfMass(superGeometry, NavierStokesLattice);

  for (std::size_t iT = 0; iT < NavierStokesConverter.getLatticeTime(maxPhysTime); ++iT) {
    //-- 5th Step: write the output data for post-processing and visualization
    postprocess(configuration, iT, timer, superGeometry, NavierStokesConverter,
      NavierStokesLattice, iTPrev, CenterOfMassPrev
    );

    //-- 6th Step: evolve the Navier-Stokes equations using FSLBM approximation
    NavierStokesLattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();
}

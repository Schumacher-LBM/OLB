/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Philipp Spelten
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

/* gausspulse3d.cpp:
 * Benchmark from C. Bogey and C. Bailly, “Three-dimensional non-reﬂective boundary conditions for acoustic simulations: far ﬁeld formulation and validation test cases,” vol. 88, 2002.
 * Method from H. Xu and P. Sagaut, “Analysis of the absorbing layers for the weakly-compressible lattice Boltzmann methods,” Journal of Computational Physics, vol. 245, pp. 14–42, Jul. 2013, doi: 10.1016/j.jcp.2013.02.051.
 * Use the option --boundaryCondition to vary how the far field is modeled:
 *  0=eternal (3x domain size to capture the vanishing pulse as ideal reference)
 *  1=periodic,
 *  2=local equilibrium BC,
 *  3=damping (with periodic boundaries)
 */

#include "olb3D.h"
#include "olb3D.hh"
#include "../noiseauxiliary.h"

using namespace olb;
using namespace olb::descriptors;

using T              = FLOATING_POINT_TYPE;
using DESCRIPTOR     = D3Q19<>;
using BulkDynamics   = BGKdynamics<T, DESCRIPTOR>;
using SpongeDynamics = SpongeLayerDynamics<T, DESCRIPTOR, momenta::BulkTuple, equilibria::SecondOrder>;

const int ndim = 3; // a few things (e.g. SuperSum3D) cannot be adapted to 2D, but this should help speed it up

typedef enum { eternal, periodic, local, damping } BoundaryType;

// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter, SuperGeometry<T, ndim>& superGeometry,
                     IndicatorF3D<T>& domainFluid, BoundaryType boundarytype)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << std::endl << "Prepare Geometry ..." << std::endl;

  // all nodes to temporary type
  superGeometry.rename(0, 2);

  T dx = converter.getConversionFactorLength();

  switch (boundarytype) {
  // eternal and damping: 3 is the actual fluid; periodic: 1 is the fluid
  case eternal:
  case damping:
    superGeometry.rename(2, 3, domainFluid);
  case periodic:
    superGeometry.rename(2, 1);
    break;
  case local: {
    superGeometry.rename(2, 1, {1, 1, 1});

    // Set material number for inflow
    Vector<T, ndim> origin = superGeometry.getStatistics().getMinPhysR(2) - dx;
    Vector<T, ndim> extend =
        superGeometry.getStatistics().getMaxPhysR(2) - superGeometry.getStatistics().getMinPhysR(2) + 2 * dx;
    extend[0] = 2 * dx;
    IndicatorCuboid3D<T> inflow(extend, origin);
    superGeometry.rename(2, 4, 1, inflow);

    // Set material number for outflow
    origin[0] = superGeometry.getStatistics().getMaxPhysR(2)[0] - dx;
    IndicatorCuboid3D<T> outflow(extend, origin);
    superGeometry.rename(2, 5, 1, outflow);

    origin    = superGeometry.getStatistics().getMinPhysR(2) - dx;
    extend    = superGeometry.getStatistics().getMaxPhysR(2) - superGeometry.getStatistics().getMinPhysR(2) + 2 * dx;
    extend[1] = 2 * dx;
    IndicatorCuboid3D<T> bottom(extend, origin);
    superGeometry.rename(2, 6, 1, bottom);
    origin[1] = superGeometry.getStatistics().getMaxPhysR(2)[1] - dx;
    IndicatorCuboid3D<T> top(extend, origin);
    superGeometry.rename(2, 6, 1, top);

    origin    = superGeometry.getStatistics().getMinPhysR(2) - dx;
    extend    = superGeometry.getStatistics().getMaxPhysR(2) - superGeometry.getStatistics().getMinPhysR(2) + 2 * dx;
    extend[2] = 2 * dx;
    IndicatorCuboid3D<T> front(extend, origin);
    superGeometry.rename(2, 6, 1, front);
    origin[2] = superGeometry.getStatistics().getMaxPhysR(2)[2] - dx;
    IndicatorCuboid3D<T> back(extend, origin);
    superGeometry.rename(2, 6, 1, back);
    break;
  }
  }

  // superGeometry.clean();  // deletes some nodes although it shouldn't - this is probably to do with the assumption that 1 is fluid and only needs one layer next to it (rather than 20 for pml)
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout
      << "Materials:" << std::endl
      << "0 = no Material (default, should be immediately renamed to Check)" << std::endl
      << "1 = fluid (except damping and eternal, then it's the far field)" << std::endl
      << "2 = check (should be empty, but corners may have issues with the normal; check is therefore treated as inflow"
      << std::endl
      << "3 = far field (except damping and eternal, then it's the fluid domain)" << std::endl
      << "4 = inflow" << std::endl
      << "5 = outflow" << std::endl
      << "6 = around domain" << std::endl;

  SuperVTMwriter3D<T>              vtmWriter("gausspulse3d");
  SuperGeometryF<T, DESCRIPTOR::d> geometryF(superGeometry);
  vtmWriter.write(geometryF);

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(UnitConverter<T, DESCRIPTOR> const& converter, SuperLattice<T, DESCRIPTOR>& sLattice,
                    SuperGeometry<T, ndim>& superGeometry, T rho0, T u0, T amplitude, T alpha,
                    BoundaryType boundarytype, T dampingDepthPU, T lengthDomain, T dampingStrength)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << std::endl << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=3 --> bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator(
      {0, 1, 2, 3}); // for local bcs all around, corners remain at 2, so they are included here
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);

  switch (boundarytype) {
  case eternal:
  case periodic:
    break;
  case local:
    boundary::set<boundary::LocalVelocity>(sLattice, superGeometry, 4);
    boundary::set<boundary::LocalPressure>(sLattice, superGeometry, 5);
    boundary::set<boundary::LocalVelocity>(sLattice, superGeometry, 6);
    break;
  case damping:
    sLattice.defineDynamics<SpongeDynamics>(superGeometry.getMaterialIndicator(1));
    // === alternatively, set the boundary, which overwrites the dynamics
    // boundary::set<boundary::SpongeLayer>( sLattice, superGeometry, 1 );
    break;
  }

  // Initial velocity and density
  AnalyticalConst<ndim, T, T>    ux = AnalyticalConst<ndim, T, T>(u0);
  AnalyticalConst<ndim, T, T>    uy = AnalyticalConst<ndim, T, T>(0.);
  AnalyticalConst<ndim, T, T>    uz = AnalyticalConst<ndim, T, T>(0.);
  AnalyticalComposed<ndim, T, T> u(ux, uy, uz);
  bulkIndicator = superGeometry.getMaterialIndicator({0, 1, 2, 3, 4, 5, 6}); // define RhoU and fields everywhere
  AcousticPulse<3, T> densityProfile(rho0, amplitude, alpha);
  linePlot<ndim, T>(densityProfile, converter.getResolution(), converter.getPhysDeltaX(), "pulse_diag", "density [LU]",
                    diagonal2d);
  linePlot<ndim, T>(densityProfile, converter.getResolution(), converter.getPhysDeltaX(), "pulse_hline", "density [LU]",
                    horizontal);

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(bulkIndicator, densityProfile, u);
  sLattice.iniEquilibrium(bulkIndicator, densityProfile, u);

  sLattice.setParameter<descriptors::OMEGA>(omega);

  if (boundarytype == damping) {
    // Define far field observable fields (TODO: globally defined -> uses RAM!)
    sLattice.defineField<descriptors::UX>(bulkIndicator, ux);
    sLattice.defineField<descriptors::UY>(bulkIndicator, uy);
    sLattice.defineField<descriptors::UZ>(bulkIndicator, uz);
    AnalyticalConst<ndim, T, T> rhoField(rho0);
    sLattice.defineField<descriptors::DENSITY>(bulkIndicator, rhoField);

    // define and output damping layer parameter
    Vector<T, ndim>               domainLengths(lengthDomain);
    DampingTerm<3, T, DESCRIPTOR> sigma(dampingDepthPU, domainLengths, dampingStrength);
    // === alternatively, define sigma as a global constant
    // AnalyticalConst<3,T,T> sigma( dampingStrength );
    sLattice.defineField<descriptors::DAMPING>(bulkIndicator, sigma);

    // write the initialization to see whether the damping parameters were set correctly
    SuperVTMwriter<T, ndim> vtmWriter_init("gausspulse3d_init");
    vtmWriter_init.createMasterFile();
    SuperLatticePhysField3D<T, DESCRIPTOR, DAMPING> damping(sLattice, 1.);
    SuperLatticePhysField3D<T, DESCRIPTOR, UX>      uxField(sLattice, 1.);
    SuperLatticePhysField3D<T, DESCRIPTOR, UY>      uyField(sLattice, 1.);
    SuperLatticePhysField3D<T, DESCRIPTOR, UZ>      uzField(sLattice, 1.);
    SuperLatticePhysField3D<T, DESCRIPTOR, DENSITY> density(sLattice, 1.);
    damping.getName() = "dampingField";
    uxField.getName() = "uxField";
    uyField.getName() = "uyField";
    uzField.getName() = "uzField";
    density.getName() = "densityField";
    vtmWriter_init.addFunctor(damping);
    vtmWriter_init.addFunctor(uxField);
    vtmWriter_init.addFunctor(uyField);
    vtmWriter_init.addFunctor(uzField);
    vtmWriter_init.addFunctor(density);
    vtmWriter_init.write(0);
  }

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Sets fixed far field average values
void setFarFieldValues(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T, ndim>& superGeometry, T rho0, T u0,
                       BoundaryType boundarytype)
{
  OstreamManager                 clout(std::cout, "setFarFieldValues");
  AnalyticalConst<ndim, T, T>    ux(u0);
  AnalyticalConst<ndim, T, T>    uy(0.);
  AnalyticalConst<ndim, T, T>    uz(0.);
  AnalyticalComposed<ndim, T, T> u(ux, uy, uz);
  AnalyticalConst<ndim, T, T>    rho(rho0);

  auto farFieldIndicator = superGeometry.getMaterialIndicator({4, 5, 6, 2});
  sLattice.RhoU(farFieldIndicator, rho, u);
  sLattice.iniEquilibrium(farFieldIndicator, rho, u);
}

// write data to termimal and file system
void setPlotData(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T, ndim>& superGeometry,
                 UnitConverter<T, DESCRIPTOR> const& converter, size_t iT, Gnuplot<T>& gplot_l2_abs, T Lp0,
                 size_t fluidMaterial)
{
  // write plot data
  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  T Lpi = L2Norm<ndim, T, DESCRIPTOR>(sLattice, superGeometry, converter, fluidMaterial);
  gplot_l2_abs.setData(T(iT), Lpi / Lp0);
}

void getGraphicalResults(SuperLattice<T, DESCRIPTOR>& sLattice, UnitConverter<T, DESCRIPTOR> const& converter,
                         size_t iT, SuperGeometry<T, ndim>& superGeometry, T amplitude)
{
  const std::string name("gausspulse3d");
  if (iT == 0) {
    SuperVTMwriter3D<T> vtmWriter(name);
    // Writes geometry, cuboid no. and rank no. to file system
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR>   rank(sLattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  // vtk output
  sLattice.scheduleBackgroundOutputVTK([&, name, iT](auto task) {
    SuperVTMwriter3D<T>        vtmWriter(name);
    SuperLatticePhysVelocity3D velocityF(sLattice, converter);
    SuperLatticePhysPressure3D pressureF(sLattice, converter);
    vtmWriter.addFunctor(velocityF);
    vtmWriter.addFunctor(pressureF);
    task(vtmWriter, iT);
  });

  // output pressure image
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
  BlockReduction3D2D<T>                     pressureReduction(pressure, Vector<T, ndim>({0, 0, 1}));
  heatmap::plotParam<T>                     jpeg_ParamP;
  jpeg_ParamP.maxValue       = converter.getPhysPressure(+amplitude / 200);
  jpeg_ParamP.minValue       = converter.getPhysPressure(-amplitude / 200);
  jpeg_ParamP.colour         = "rainbow";
  jpeg_ParamP.fullScreenPlot = true;
  heatmap::write(pressureReduction, iT, jpeg_ParamP);

  std::stringstream ss;
  ss << std::setw(4) << std::setfill('0') << iT;
  T                          dist        = converter.getPhysDeltaX();
  T                          ndatapoints = converter.getResolution(); // number of data points on line
  AnalyticalFfromSuperF3D<T> pressure_interpolation(pressure, true, true);
  T                          pmin(converter.getPhysPressure(-amplitude / 50));
  T                          pmax(converter.getPhysPressure(+amplitude / 50));
  linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_hline_" + ss.str(), "pressure [PU]",
                    horizontal, false, true, pmin, pmax);
  linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_vline_" + ss.str(), "pressure [PU]", vertical,
                    false, true, pmin, pmax);
  linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_diagonal_" + ss.str(), "pressure [PU]",
                    diagonal2d, false, true, pmin, pmax);
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  CLIreader args(argc, argv);

  // === get Command Line Arguments
  std::string outdir = args.getValueOrFallback<std::string>("--outdir", "");
  // lattice definitions
  const T   charL = args.getValueOrFallback("--charL", 1.); // domain size [m]
  const int res   = args.getValueOrFallback("--res", 101);  // number of points across charL
  const T   rho0  = args.getValueOrFallback("--rho0", 1.);  // background density
  const T   Ma    = args.getValueOrFallback("--Ma", 0.1);   // characteristic lattice velocity
  const T   charV = args.getValueOrFallback("--charV", 1.); // characteristic physical velocity
  // timing and outputs
  size_t maxLatticeT = args.getValueOrFallback("--iTmax", 300); // maximum number of iterations
  size_t nplot       = args.getValueOrFallback("--nplot", 100); // minimum number of plot points
  bool   doImages    = args.contains("--doImages");             // otherwise won't do images
  size_t nout        = args.getValueOrFallback("--nout", 5);    // minimum number of vtk outputs
  size_t iTout       = args.getValueOrFallback("--iTout", 0);   // vtk outputs every iTout iterations
  T      tout        = args.getValueOrFallback("--tout", 0.);   // vtk outputs every tout seconds
  // boundary condition variation
  const int boundaryCondition =
      args.getValueOrFallback("--boundaryCondition", 3); // boundaryCondition (see beginning of file)
  const size_t dampingDepthLU = args.getValueOrFallback("--dampingDepthLU", 20); // number of points for sponge layer
  const T      dampingStrength =
      args.getValueOrFallback("--dampingStrength", 1.); // maximum damping strength of sponge layer
  const int eternalscale =
      args.getValueOrFallback("--eternalscale", 3); // how much to extend the domain for "eternal" layer
  // Gauss pulse
  const T amplitude = args.getValueOrFallback("--amplitude", 1e-3);  // physical pressure amplitude of Gauss pulse
  const T bPulse    = args.getValueOrFallback("--bPulse", 1. / 20.); // parameter for the Gauss pulse shape

  // === set output directory depending on input values
  std::stringstream outdir_mod;
  if (outdir == "") {
    outdir_mod << "./tmp";
    switch (boundaryCondition) {
    case 0:
      outdir_mod << "_eternal";
      break;
    case 1:
      outdir_mod << "_periodic";
      break;
    case 2:
      outdir_mod << "_local";
      break;
    case 3:
    default:
      outdir_mod << "_damping";
      break;
    }
    if (Ma != 0.1)
      outdir_mod << "_Ma" << Ma;
    if (amplitude != 1e-3)
      outdir_mod << "_a" << amplitude;
    if (charL != 1.)
      outdir_mod << "_l" << charL;
    if (res != 101)
      outdir_mod << "_res" << res;
    if (boundaryCondition == 0)
      outdir_mod << "_scale" << eternalscale;
    if (boundaryCondition == 3)
      outdir_mod << "_bd" << dampingDepthLU << "x" << dampingStrength;
  }
  else {
    outdir_mod << outdir;
  }

  // === set OLB output directory and create a TeeBuffer that writes to both the console and the file
  singleton::directories().setOutputDir(outdir_mod.str() + "/");
  std::ofstream   fileStream(outdir_mod.str() + "/output.txt");
  DoubleBuffer    doubleBuffer(std::cout.rdbuf(), fileStream.rdbuf());
  std::streambuf* originalCoutBuffer = std::cout.rdbuf(&doubleBuffer);
  OstreamManager  clout(std::cout, "main");
  clout << std::endl << "outdir set to " << outdir_mod.str() << std::endl;
  T alpha = log(2.) / (bPulse * bPulse);
  clout << "Gauss pulse: b=" << bPulse << "; alpha=" << alpha << std::endl;

  // === set boundarytype depending on input values
  BoundaryType boundarytype;
  switch (boundaryCondition) {
  case 0:
    boundarytype = eternal;
    clout << "Boundary condition is solved by just extending the domain to " << eternalscale << " times" << std::endl;
    break;
  case 1:
    boundarytype = periodic;
    clout << "Boundary condition type specified to periodic." << std::endl;
    break;
  case 2:
    boundarytype = local;
    clout << "Boundary condition type specified to local." << std::endl;
    break;
  case 3:
  default:
    boundarytype = damping;
    clout << "Boundary condition type specified to periodic." << std::endl;
    break;
  }

  // === initialize the unit converter
  clout << std::endl
        << "Input to unitconverter: res=" << res << ", charL=" << charL << ", charV=" << charV << ", Ma=" << Ma
        << ", rho0=" << rho0 << std::endl;
  T                                                                  charVLU = Ma / std::sqrt(3.);
  UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> const converter((size_t)res,    // resolution
                                                                               (T)charVLU,     // charLatticeVelocity
                                                                               (T)charL,       // charPhysLength
                                                                               (T)charV,       // charPhysVelocity
                                                                               (T)T(1.48e-05), // physViscosity
                                                                               (T)rho0         // physDensity
  );
  converter.print();

  // === change domain size depending the boundary condition
  T dampingDepthPU = converter.getPhysLength(dampingDepthLU);
  T lengthDomain;
  switch (boundarytype) {
  case eternal:
    lengthDomain = eternalscale * charL;
    break; // extend the domain
  case periodic:
    lengthDomain = charL;
    break; // reference domain size
  case local:
    lengthDomain = charL + 2 * converter.getConversionFactorLength();
    break; // add one layer in each direction
  case damping:
    lengthDomain = charL + 2 * dampingDepthPU;
    break; // add boundary layer
  }
  clout << "Fluid Domain = " << charL << "^3; Simulation Domain = " << lengthDomain << "^3" << std::endl;
  Vector<T, ndim>      originDomain(-lengthDomain / 2.);
  Vector<T, ndim>      extendDomain(lengthDomain);
  IndicatorCuboid3D<T> domainSim(extendDomain, originDomain);
  Vector<T, ndim>      originFluid(-T(0.5) * charL);
  Vector<T, ndim>      extendFluid(charL);
  IndicatorCuboid3D<T> domainFluid(extendFluid, originFluid);
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 8;
#endif
  CuboidDecomposition<T, ndim> CuboidDecomposition(domainSim, converter.getConversionFactorLength(), noOfCuboids);
  switch (boundarytype) {
  case eternal:
  case periodic:
  case damping:
    CuboidDecomposition.setPeriodicity({true, true, true});
    break;
  case local:
    CuboidDecomposition.setPeriodicity({false, false, false});
    break;
  }

  // === create loadBalancer and superGeometry
  HeuristicLoadBalancer<T> loadBalancer(CuboidDecomposition);
  // Instantiation of a superGeometry
  size_t overlap = 3;
  if (noOfCuboids > 1 && (boundarytype == damping)) {
    clout << "overlap=" << overlap << "; dampingDepth_LU=" << dampingDepthLU
          << "; setting overlap>=dampingDepthLU to allow the boundary to be so large." << std::endl;
    overlap = std::max(overlap, dampingDepthLU);
  }
  SuperGeometry<T, ndim> superGeometry(CuboidDecomposition, loadBalancer, overlap);

  // === set fluidMaterial to evaluate the core domain correctly
  prepareGeometry(converter, superGeometry, domainFluid, boundarytype);
  size_t fluidMaterial = 3;
  switch (boundarytype) { // for correct normals, fluid material number must be 1
  case periodic:
  case local:
    fluidMaterial = 1;
    break;
  default:
    break;
  }
  clout << "Actual fluid material number is " << fluidMaterial
        << "; number of fluid voxels: " << superGeometry.getStatistics().getNvoxel(fluidMaterial) << std::endl;

  // === prepare Lattice and set boundaryConditions
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(converter, sLattice, superGeometry, rho0, charVLU, amplitude, alpha, boundarytype, dampingDepthPU,
                 lengthDomain, dampingStrength);

  // === Initialize pressure L2 norm plot
  Gnuplot<T> gplot_l2_abs("l2_absolute");
  gplot_l2_abs.setLabel("time []", "absolute L2 norm []");
  T Lp0 = L2Norm<ndim, T, DESCRIPTOR>(sLattice, superGeometry, converter, fluidMaterial);

  // === calculate output intervals
  // nout is the minimum number of vtk outputs --> take max between nout and nout derived from iTout or tout
  size_t nout_from_iTout = 0, nout_from_tout = 0;
  if (iTout != 0) {
    nout_from_iTout = size_t(maxLatticeT / iTout);
    nout            = std::max(nout, nout_from_iTout);
  }
  if (tout != 0) {
    nout_from_tout = size_t(maxLatticeT / tout);
    nout           = std::max(nout, nout_from_tout);
  }
  size_t iTvtk  = std::max(int(maxLatticeT / nout), 1);
  size_t iTplot = std::min(std::max(int(maxLatticeT / nplot), 1), 25);

  clout << "Timing setup:" << std::endl
        << "maxLatticeT=" << maxLatticeT << "; maxPhysT=" << converter.getPhysTime(maxLatticeT)
        << "; dt=" << converter.getPhysDeltaT() << std::endl
        << "iTout=" << iTout << "; tout=" << tout << "; iTvtk=" << iTvtk << "; iTplot=" << iTplot << std::endl;

  // === 4th Step: Main Loop with Timer ===
  clout << std::endl << "starting simulation..." << std::endl;
  util::Timer<T> timer(maxLatticeT, superGeometry.getStatistics().getNvoxel());
  timer.start();

  size_t iT = 0;
  while (iT < maxLatticeT) {
    // === Definition of Initial and Boundary Conditions ===
    if (boundarytype == local) {
      setFarFieldValues(sLattice, superGeometry, rho0, charVLU, boundarytype);
    }
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);
    // === Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === Computation and Output of the Results ===
    if (iT % iTplot == 0)
      setPlotData(sLattice, superGeometry, converter, iT, gplot_l2_abs, Lp0, fluidMaterial);
    if (doImages && iT % iTvtk == 0)
      getGraphicalResults(sLattice, converter, iT, superGeometry, amplitude);
    if (iT % iTvtk == 0) {
      // lattice and timer statistics console output
      timer.update(iT);
      timer.printStep();
      sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    }
    iT++;
  }
  clout << "Simulation stopped after " << iT << "/" << maxLatticeT << " iterations." << std::endl;

  // output the development of pressure L2 norm over time
  gplot_l2_abs.setYrange(1e-5, 1);
  gplot_l2_abs.setLogScale(2);
  gplot_l2_abs.writePNG(-1, -1, "gplot_l2_abs");

  timer.stop();
  timer.printSummary();
  std::cout.rdbuf(originalCoutBuffer);
}

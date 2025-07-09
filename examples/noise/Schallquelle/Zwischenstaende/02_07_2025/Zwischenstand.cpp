
/*Notiz:
Das Programm kompiliert Fehlerfrei und kann durchgeführt werden. Die Durchführung benötigt viel Zeit und ausgegeben wird ein Würfel, der auf 4 Seiten von innen nach außen eine Druckerhöhung erfährt. Die Ergebnisse sind nicht wie erwartet, aber ein erster Schritt.
Bei dem ausführen des Programms muss mit --iTmax XXX eine Zahl angegeben werden, wieviele Iterationsdurchläufe das Programm durchläuft
Terminalbefehl: make; ./cavity2d --iTmax 30*/

/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 * Dieses Programm simuliert eine Schallquell in einem 3D Medium.
 * Dabei soll sich eine isolierte frei schwingende Sinuswelle im Raum verteilen
 *
 */



#include <olb.h>
#include "olb3D.h"
#include "olb3D.hh"
#include "../noiseauxiliary.h"
#include "SchallwelleGeschwindigkeit.h"
#include "SchallwelleDruck.h"

using namespace olb;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = descriptors::D3Q19<>;
using BulkDynamics   = BGKdynamics<T, DESCRIPTOR>;
using SpongeDynamics = SpongeLayerDynamics<T, DESCRIPTOR, momenta::BulkTuple, equilibria::SecondOrder>;

const int ndim = 3; // a few things (e.g. SuperSum3D) cannot be adapted to 2D, but this should help speed it up

const T physDeltaX        = 0.02;   // grid spacing [m]
const T physDeltaT        = 0.00078125;  // temporal spacing [s]
const T physLength        = 1.0;         // length of the squared cuboid [m]
const T physLidVelocity   = 1.0;         // velocity imposed on lid [m/s]
const T physViscosity     = 0.001;        // kinetic viscosity of fluid [m*m/s]
const T physDensity       = 1.0;         // fluid density [kg/(m*m*m)]


const T physMaxT          = 0.5;        // maximal simulation time [s]

typedef enum { periodic, local } BoundaryType;


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
  case periodic:
    superGeometry.rename(2, 1);
    break;
  case local: {
    superGeometry.rename(2, 1);
    Vector<T,3> center(0., 0.,0.); // Zentrum der Welle
    T radius = dx;
    IndicatorSphere3D<T> pointSource(center, radius);
    superGeometry.rename(1,3,pointSource);
    //superGeometry.rename(1,3,1,pointSource);
    break;
  }
  }

  superGeometry.getStatistics().print();
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
        {1,3}); // for local bcs all around, corners remain at 2, so they are included here
    sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
  
    switch (boundarytype) {
    case periodic: {
      // Fuer die freie Welle nutzen
      // Verortung der Schallquelle
      Vector<T,3> center(0., 0.,0.); // Zentrum der Welle

      //Schallquelle
      SchallwelleDru<3,T> SchallwelleDru(1e-3,0.314,0.314,0.,0, center);//SchallwelleDruck(T amplitude, T Wellenzahl, T phase, Vector<T, ndim> x0 = Vector<T, ndim>(0.))
      //SchallwelleGesch<3,T> SchallwelleGeschwind(1e-3,0.314, 0. , 1., 0.577, center); // T amplitude, T Wellenzahl, T phase,T rho0,T Geschwindigkeit,Vector<T, ndim> x0 = Vector<T,  ndim>(0.))
      // TODO: rhoF durch std::sin(iT)

      //AnalyticalConst3D<T,T> rhoF(SchallwelleDru);       // Definiert rho auf 1
      AnalyticalConst3D<T,T> uWall(0, 0, 0);   // Definiert die Geschwindigkeit überall auf 0 in X und Y Richtung
      AnalyticalConst3D<T,T> uLid(0, 0, 0);    // Definiert die Geschwindigkeit an der lid auf eine Geschwindigkeit in X Richtung und 0 in Y Richtung

      break;
    }
    case local: // Fuer die erzwungene Welle nutzen
      //boundary::set<boundary::LocalPressure>(sLattice, superGeometry, 3);
      break;
    }

    sLattice.setParameter<descriptors::OMEGA>(omega);
  
    // Make the lattice ready for simulation
    sLattice.initialize();
  
    clout << "Prepare Lattice ... OK" << std::endl;
  }

// Hier werden die Startbedingungen für eine Gausschen Puls definiert
void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
  SuperLattice<T, DESCRIPTOR>& sLattice,
  std::size_t iT, SuperGeometry<T,ndim>& superGeometry)
{
  auto domain = superGeometry.getMaterialIndicator({3}); // Die Materialien werden Fluid (1), gleitfreie Grenzen (2) und Geschiwindigkeit gesetzt (3) Quelle: Handbuch

  AnalyticalConst3D<T,T> rhoF(1.+0.1*std::sin(iT*0.025*std::numbers::pi_v<T>/4));       // assuming T=40 iterations
  AnalyticalConst3D<T,T> uInf(0, 0,0);                                                    // Definiert die Geschwindigkeit überall auf 0 in X und Y Richtung
  
  if(iT==0)
  {
    AnalyticalConst3D<T,T> rhoF(1.);       // Definiert rho auf 1
    AnalyticalConst3D<T,T> uInf(0, 0,0);   // Definiert die Geschwindigkeit überall auf 0 in X und Y Richtung
  }

    sLattice.defineRhoU(domain, rhoF, uInf);
    // Initialize populations to equilibrium state
    sLattice.iniEquilibrium(domain, rhoF, uInf);

    // ALERNATIVE: define only Rho, find u from dynamics, then iniEqui
    //sLattice.defineRho(domain, rhoF);
    //auto uInf = sLattice.getU(domain);
    //sLattice.iniEquilibrium(domain, rhoF, uInf);
}  // setBoundaryValues

void getGraphicalResults(SuperLattice<T, DESCRIPTOR>& sLattice, UnitConverter<T, DESCRIPTOR> const& converter,
  size_t iT, SuperGeometry<T, ndim>& superGeometry, T amplitude)
{
    const std::string name("Schallquelle");
    if (iT == 0) {
      SuperVTMwriter3D<T> vtmWriter(name);
      // Writes geometry, cuboid no. and rank no. to file system
      SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
      SuperLatticeRank3D<T, DESCRIPTOR>   rank(sLattice);
      vtmWriter.write(cuboid);
      vtmWriter.write(rank);
      vtmWriter.createMasterFile();
    }  // iT==0

    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // vtk output
    sLattice.scheduleBackgroundOutputVTK([&, name, iT](auto task) {
    SuperVTMwriter3D<T>        vtmWriter(name);
    SuperLatticePhysVelocity3D velocityF(sLattice, converter);
    SuperLatticePhysPressure3D pressureF(sLattice, converter);
    vtmWriter.addFunctor(velocityF);
    vtmWriter.addFunctor(pressureF);
    task(vtmWriter, iT);
    });  // scheduleBackgroundOutputVTK

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
                      horizontal, false, false, pmin, pmax);  // TODO setRange=true (before pmin, pmax)
    //linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_vline_" + ss.str(), "pressure [PU]", vertical,
    //                    false, true, pmin, pmax);
    //linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_diagonal_" + ss.str(), "pressure [PU]",
    //                    diagonal2d, false, true, pmin, pmax);
}  // getGraphicalResults






int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  CLIreader args(argc, argv);
  size_t maxLatticeT = args.getValueOrFallback("--iTmax", 0); // maximum number of iterations
  
  
  OstreamManager clout( std::cout,"main" ); // writing all output first in a userdefined Buffer of type OMBuf. On a flush it spits out at first the userdefined text in squared brackets and afterwards everything from the buffer

  
  // Provide the unit converter the characteristic entities
  const UnitConverter<T,DESCRIPTOR> converter (
    physDeltaX,        // physDeltaX: spacing between two lattice cells in [m]
    physDeltaT,        // physDeltaT: time step in [s]
    physLength,        // charPhysLength: reference length of simulation geometry in [m]
    physLidVelocity,   // charPhysVelocity: highest expected velocity during simulation in [m/s]
    physViscosity,     // physViscosity: physical kinematic viscosity in [m^2/s]
    physDensity        // physDensity: physical density [kg/m^3]
  );
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  BoundaryType boundarytype = local;
  Vector<T,ndim> originFluid(-0.5, -0.5, -0.5);
  Vector<T,ndim> extendFluid(physLength, physLength, physLength);
  IndicatorCuboid3D<T> domainFluid(extendFluid, originFluid);

  Vector<T,ndim> extend{physLength, physLength, physLength};
  Vector<T,ndim> origin{-0.5, -0.5, -0.5};
  IndicatorCuboid3D<T> cuboid(extend, origin);
  CuboidDecomposition3D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());
  cuboidDecomposition.setPeriodicity({true,true,true});
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  
  SuperGeometry<T,ndim> superGeometry(cuboidDecomposition, loadBalancer);
  prepareGeometry(converter, superGeometry,domainFluid,boundarytype);

  // === 3rd Step: Prepare Lattice ===
  // Vorläufige Loesung mit zufälligen Zahlen. Bitte die richtigen Zahlen noch hinzufuegen!
  T rho0 = 1.0;
  T u0 = converter.getCharLatticeVelocity();
  T amplitude = 1e-3;
  T alpha = 0.314;  // oder berechnet aus einem Pulsparameter
  T dampingDepthPU = 0.1;
  T lengthDomain = physLength;
  T dampingStrength = 1.0;

  SuperLattice<T,DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(converter, sLattice, superGeometry, rho0, u0, amplitude, alpha,
    boundarytype, dampingDepthPU, lengthDomain, dampingStrength);

  // === 4th Step: Main Loop with Timer ===
  std::size_t iTmax = converter.getLatticeTime(physMaxT);
  if ( maxLatticeT != 0 ) iTmax = maxLatticeT;
  std::size_t iTvtk = int(std::max(iTmax/5., 1.));
  std::size_t iTtimer = int(std::max(iTmax/20., 1.));
  
  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    if (boundarytype == local) setBoundaryValues(converter, sLattice, iT, superGeometry);
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    if ( iT%iTvtk == 0 ) getGraphicalResults(sLattice, converter, iT, superGeometry, amplitude);
    if ( iT%iTtimer == 0 ) {timer.update(iT); timer.printStep();}
  }

  timer.stop();
  timer.printSummary();
}



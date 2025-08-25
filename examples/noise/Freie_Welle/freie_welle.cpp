
/*Notiz: FREIE WELLE
Das Programm kompiliert Fehlerfrei und kann durchgeführt werden. 
Bei dem ausführen des Programms muss mit --iTmax XXX eine Zahl angegeben werden, wieviele Iterationsdurchläufe das Programm durchläuft
Terminalbefehl: make clean; make; ./Freie_Welle --iTmax 150 --outdir tmp_periodic_05  --> Zusätzlich kann die Amplitude angegeben werden
To Do: Schallgeschw. und Amplitude an die Werte von Luft anpassen und Quellen finden. Außerdem die Größe des Mediums angeben. Messwerte nehmen
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 * Dieses Programm simuliert eine Schallquell in einem 3D Medium.
 * Dabei soll sich eine isolierte frei schwingende Sinuswelle im Raum verteilen
 *
 */

#define FEATURE_REPORTER

#include <olb.h>
#include "olb3D.h"
#include "olb3D.hh"
#include "../noiseauxiliary.h"
#include "SchallwelleGeschwindigkeit.h"
#include "SchallwelleDichte.h"
#include "functors/analytical/analyticalF.hh"

using namespace olb;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = descriptors::D3Q19<>;
using BulkDynamics   = BGKdynamics<T, DESCRIPTOR>;
using SpongeDynamics = SpongeLayerDynamics<T, DESCRIPTOR, momenta::BulkTuple, equilibria::SecondOrder>;

const int ndim = 3; // a few things (e.g. SuperSum3D) cannot be adapted to 2D, but this should help speed it up

const T physDeltaX        = 0.02;   // grid spacing [m]
const T physLength        = 2.4;         // length of the cuboid [m]
const T physspan          = 0.46;
const T physwidth         = 0.46;
const T physLidVelocity   = 0.57735;         // velocity imposed on lid [m/s] Fuer die Machzahl relevant (Vorher 1.0, jetzt ein zehntel der Schallgeschwindigkeit)
const T physViscosity     = 0.001392;     // kinetic viscosity of fluid [m*m/s] Fuer die Relaxationszeit verantwortlich
const T physDensity       = 1.189;         // fluid density of air (20°C)[kg/(m*m*m)]
const T physMaxT          = 0.5;        // maximal simulation time [s]
const T physDeltaT        = physDeltaX / 343.46;//((0.68255-0.5)/3)/physViscosity*physDeltaX*physDeltaX;// 0,68255, weil Tau 0,68255 sein soll. Vorher: physDeltaX/343.46;  // temporal spacing [s] t=physDeltaX/c_s (Vorher 0.00078125, Jetzt: 5,8e-5)
typedef enum { periodic, local } BoundaryType;
 


struct PressureO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    auto particle = cells.template get<names::Points>();
    const V rho = cells.template get<names::NavierStokes>().computeRho();
    const V pressure = util::pressureFromDensity<V,DESCRIPTOR>(rho);
    particle.template setField<descriptors::SCALAR>(pressure);
  }
};




 
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
   case periodic:{
     superGeometry.rename(2, 1);
     break;
   } 
     
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
 }//prepareGeometry
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
   } //prepareLattice
 
 void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
   SuperLattice<T, DESCRIPTOR>& sLattice,
   std::size_t iT, SuperGeometry<T,ndim>& superGeometry,BoundaryType boundarytype, T amplitude, T rho0)
 {
 
   if (boundarytype == local) {
     auto domain = superGeometry.getMaterialIndicator({3});
 
     // 2 Sinusperioden in 40 Zeitschritten
     T Kreisfrequenz = 2. * std::numbers::pi_v<T> /40.0;
     // optionaler Einschwingfaktor
     T envelope =std::sin(std::min(1.0, iT / 40.0) * std::numbers::pi_v<T> / 2.0);
 
     // Sinus-Anregung
     AnalyticalConst3D<T,T> rhoF(1. + envelope * 1e-3 * std::sin(iT * Kreisfrequenz)); // Hier überprüfen ob das nicht eine Amplitude ist!
     AnalyticalConst3D<T,T> uInf(0., 0., 0.);
    
      if(iT==0)
      {
        AnalyticalConst3D<T,T> rhoF(1.);       // Definiert rho auf 1
        AnalyticalConst3D<T,T> uInf(0, 0,0);   // Definiert die Geschwindigkeit überall auf 0 in X und Y Richtung
      }

     sLattice.defineRhoU(domain, rhoF, uInf);
     sLattice.iniEquilibrium(domain, rhoF, uInf);
   }
 
   if (boundarytype == periodic && iT==0) {
     auto domain = superGeometry.getMaterialIndicator({1});
    
     T wellenzahl=2. * std::numbers::pi_v<T>/0.7 ;//k=2pi/lamda
     T kreisfrequenz= 343.45*wellenzahl;//2. * std::numbers::pi_v<T> * 2.0 / 40.0; // Theoretisch muesste die Kreisfrequenz bei 4 pi liegen (w=cs*k).  :40 wird gerechnet, weil 40 Iterationen getätigt werden sollen um auf 2 Perioden zu kommen
     T phase =0.;
     //T time = converter.getPhysTime(iT);  // physikalische Zeit aus Lattice-Zeit
     T time= iT;
     T cs=sqrt(T(1)/descriptors::invCs2<T,DESCRIPTOR>());

     olb::SchallwelleRho<3, T, DESCRIPTOR> schallquelle(rho0, amplitude, wellenzahl, kreisfrequenz, phase, time, converter);
     olb::SchallwelleGesch<3,T, DESCRIPTOR> schallquelle_geschwindigkeit(amplitude,wellenzahl,kreisfrequenz,phase,time,rho0,cs,converter);
     // AnalyticalConst3D<T,T> rhoF(schallquelle);
     AnalyticalConst3D<T,T> uInf(0., 0., 0.);
 
 
      // Hier wird die Geschwindigkeit im Terminal ausgegeben
    //   Vector<T,3> punkt = {0.01, 0.0, 0.0};  // Punkt, an dem ausgewertet wird
    //   T u[3];  // Ergebnis wird hier gespeichert
    //   schallquelle_geschwindigkeit(u, punkt.data());
    //   std::cout << "[iT=" << iT << ", t=" << time << "s] Geschwindigkeit an "<<punkt<<": "
    //   << "u = (" << u[0] << ", " << u[1] << ", " << u[2] << ")\n";

    //  // Hier wird der Druck im Terminal ausgegeben
     
    //  T p[3];  // Ergebnis wird hier gespeichert
    //  schallquelle(p, punkt.data());
    //  std::cout << "[iT=" << iT << ", t=" << time << "s] Druck an "<<punkt<<": "
    //  << "rho = (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
 
     // sLattice.defineRhoU(domain, schallquelle, uInf);
     // sLattice.iniEquilibrium(domain, schallquelle, uInf);
     sLattice.defineRhoU(domain, schallquelle, schallquelle_geschwindigkeit);
     sLattice.iniEquilibrium(domain, schallquelle, schallquelle_geschwindigkeit);
   }
 } //setBoundaryValues
 
//  void getGraphicalResults(SuperLattice<T, DESCRIPTOR>& sLattice, UnitConverter<T, DESCRIPTOR> const& converter,
//    size_t iT, SuperGeometry<T, ndim>& superGeometry, T amplitude)
//  {
//      const std::string name("Schallquelle");
//      if (iT == 0) {
//        SuperVTMwriter3D<T> vtmWriter(name);
//        // Writes geometry, cuboid no. and rank no. to file system
//        SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
//        SuperLatticeRank3D<T, DESCRIPTOR>   rank(sLattice);
//        vtmWriter.write(cuboid);
//        vtmWriter.write(rank);
//        vtmWriter.createMasterFile();
//      }  // iT==0
 
//      sLattice.setProcessingContext(ProcessingContext::Evaluation);
 
//      // vtk output
//      sLattice.scheduleBackgroundOutputVTK([&, name, iT](auto task) {
//      SuperVTMwriter3D<T>        vtmWriter(name);
//      SuperLatticePhysVelocity3D velocityF(sLattice, converter);
//      SuperLatticePhysPressure3D pressureF(sLattice, converter);
//      vtmWriter.addFunctor(velocityF);
//      vtmWriter.addFunctor(pressureF);
//      task(vtmWriter, iT);
//      });  // scheduleBackgroundOutputVTK
 
//      //output pressure image
//      SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
//      BlockReduction3D2D<T>                     pressureReduction(pressure, Vector<T, ndim>({0, 0, 1}));
//      heatmap::plotParam<T>                     jpeg_ParamP;
//      jpeg_ParamP.maxValue       = converter.getPhysPressure(+amplitude / 200);
//      jpeg_ParamP.minValue       = converter.getPhysPressure(-amplitude / 200);
//      jpeg_ParamP.colour         = "rainbow";
//      jpeg_ParamP.fullScreenPlot = true;
//      heatmap::write(pressureReduction, iT, jpeg_ParamP);
 
//      std::stringstream ss;
//      ss << std::setw(4) << std::setfill('0') << iT;
//      T                          dist        = converter.getPhysDeltaX();
//      T                          ndatapoints = converter.getResolution(); // number of data points on line
//      AnalyticalFfromSuperF3D<T> pressure_interpolation(pressure, true, true);
//      T                          pmin(converter.getPhysPressure(-amplitude / 50));
//      T                          pmax(converter.getPhysPressure(+amplitude / 50));
//      linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_hline_" + ss.str(), "pressure [PU]",
//                        horizontal, false, false, pmin, pmax);  // TODO setRange=true (before pmin, pmax)
//      //linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_vline_" + ss.str(), "pressure [PU]", vertical,
//      //                    false, true, pmin, pmax);
//      // linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_diagonal_" + ss.str(), "pressure [PU]",
//      //                     diagonal2d, false, false, pmin, pmax);
//  }  // getGraphicalResults
 
 
int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  CLIreader args(argc, argv);
  std::string outdir = args.getValueOrFallback<std::string>("--outdir", "");
  outdir += "_reporter";
  size_t maxLatticeT = args.getValueOrFallback("--iTmax", 0); // maximum number of iterations
  T amplitude = args.getValueOrFallback("--a", 1e-3); // maximum number of iterations
  if (outdir == "") outdir = "./tmp/";
  else outdir = "./" + outdir + "/";
  singleton::directories().setOutputDir(outdir);
  
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
  BoundaryType boundarytype = periodic;
  Vector<T,ndim> originFluid(-physLength/2., -physwidth/2., -physspan/2.);
  Vector<T,ndim> extendFluid(physLength, physwidth, physspan);
  IndicatorCuboid3D<T> domainFluid(extendFluid, originFluid);
  // -----------Variabeln definiere Messungen
  size_t nplot                  = args.getValueOrFallback( "--nplot",             100 );  
  size_t iTout                  = args.getValueOrFallback( "--iTout",             0   );  
    
  //-----------------------------
  Vector<T,ndim> extend{physLength, physwidth, physspan};
  Vector<T,ndim> origin{-physLength/2., -physwidth/2., -physspan/2.};
  IndicatorCuboid3D<T> cuboid(extend, origin);
  CuboidDecomposition3D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());
  cuboidDecomposition.setPeriodicity({true,true,true});
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  
  SuperGeometry<T,ndim> superGeometry(cuboidDecomposition, loadBalancer);
  prepareGeometry(converter, superGeometry, domainFluid, boundarytype);

  // === 3rd Step: Prepare Lattice ===
  T rho0 = 1.0; 
  T u0 = converter.getCharLatticeVelocity();
  T alpha = 0.314;  // oder berechnet aus einem Pulsparameter
  T dampingDepthPU = 0.1;
  T lengthDomain = physLength;
  T dampingStrength = 1.0;
 
  SuperLattice<T,DESCRIPTOR> sLattice(superGeometry);
  prepareLattice( converter, sLattice, superGeometry, rho0, u0, amplitude, alpha,
                  boundarytype, dampingDepthPU, lengthDomain, dampingStrength);
 
  // === 4th Step: Main Loop with Timer ===
  std::size_t iTmax = converter.getLatticeTime(physMaxT);
  if ( maxLatticeT != 0 ) iTmax = maxLatticeT;
 // T vtkanzahl=iTmax;///10.;//*10.;
  std::size_t iTvtk = 1;//int(std::max(iTmax/vtkanzahl, 1.));
  std::size_t iTtimer = int(std::max(iTmax/20., 1.));
  // === calculate output intervals
  // nout is the minimum number of vtk outputs --> take max between nout and nout derived from iTout or tout
  
  size_t iTplot = std::min(std::max(int(maxLatticeT / nplot), 1), 25);

  clout << "Timing setup:" << std::endl
        << "maxLatticeT=" << maxLatticeT << "; maxPhysT=" << converter.getPhysTime( maxLatticeT ) << "; dt=" << converter.getPhysDeltaT() << std::endl
        << "iTout=" << iTout << "; tout=" << converter.getPhysTime( iTout ) << "; iTvtk=" << iTvtk << "; iTplot=" << iTplot << std::endl;

   util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
   timer.start();
  // Zwischenschritt: Messwerte nehmen //-----------Vorgegeben---
  // #if defined(FEATURE_TWOD)
  // Vector<T,ndim> measurePhysR{-0.0875,0.0052};
  // Vector<int,3> measureLatticeR{};
  //#elif defined(FEATURE_THREED)
  Vector<T,ndim> measurePhysR{0.566,0.0,0.230};
  Vector<int,4> measureLatticeR{};
  //#endif
  // get nearest lattice point to measurePhysR
  if (auto latticeR = cuboidDecomposition.getLatticeR(measurePhysR)) {
    measureLatticeR = *latticeR;
  }
  LatticeR<2> measureWatchpointR;
  // #if defined(FEATURE_TWOD)
  // SuperD<T,descriptors::D2<fields::PHYS_R,descriptors::SCALAR>> watchpointsD(loadBalancer);
  // #elif defined(FEATURE_THREED)
  SuperD<T,descriptors::D3<fields::PHYS_R,descriptors::SCALAR>> watchpointsD(loadBalancer);
  if (loadBalancer.isLocal(measureLatticeR[0])) {
    auto& blockD = watchpointsD.getBlock(measureLatticeR[0]);
    // #if defined(FEATURE_TWOD)
    // blockD.resize({blockD.getNx()+1,1,1});
    // #elif defined(FEATURE_THREED)
    blockD.resize({blockD.getNx()+1,1,1});
    //#endif
    auto watchpoint = blockD.get(blockD.getNcells()-1);
    watchpoint.template setField<fields::PHYS_R>(measurePhysR);
    measureWatchpointR[0] = measureLatticeR[0];
    measureWatchpointR[1] = blockD.getNcells()-1;
  }
  watchpointsD.setProcessingContext(ProcessingContext::Simulation);

  SuperLatticePointCoupling pressureO(PressureO{},
                                      names::NavierStokes{}, sLattice,
                                      names::Points{}, watchpointsD);

  CSV<T> csvWriter("Welle", ';', {"iT", "t", "p"}, ".csv");

  if (measureWatchpointR[1] == 0 && singleton::mpi().getRank() == 0) {
    std::cout << "[WARNUNG] Messpunkt wurde NICHT lokal gefunden! Druckwerte bleiben leer." << std::endl;
  }
  if (boundarytype == periodic) {setBoundaryValues(converter, sLattice, 0, superGeometry, boundarytype, amplitude, rho0);}
  //--------------------------------------- FOR SCHLEIFE-------------------------------------------------------------------------------------------
  for (std::size_t iT=0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    if (boundarytype == local) {setBoundaryValues(converter, sLattice, iT, superGeometry, boundarytype, amplitude,rho0);}
    
    // ------------------------- Messwerte nehmen

    std::vector<T> measurements(1);
    pressureO.execute();
    watchpointsD.setProcessingContext(ProcessingContext::Evaluation);
    if (loadBalancer.isLocal(measureWatchpointR[0])) {
      auto measureCell = watchpointsD.getBlock(loadBalancer.loc(measureWatchpointR[0]))
                                    .get(measureWatchpointR[1]);
      T measurePressure = measureCell.template getField<descriptors::SCALAR>();
      measurements[0] += converter.getPhysPressure(measurePressure);
    }
    std::vector<T> globalMeasurements(1);

    #ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceVect(measurements, globalMeasurements, MPI_SUM);
    singleton::mpi().bCast(globalMeasurements.data(), globalMeasurements.size());
    #else
    globalMeasurements = measurements;
    #endif
    csvWriter.writeDataFile(iT, {converter.getPhysTime(iT), T{globalMeasurements[0]}});
    
     
    

    // std::cout << "[iT=" << iT << ", t=" << converter.getPhysTime(iT) << "s] "
    //           << "Gemessener Druck am Punkt ("
    //           << measurePhysR[0] << ", "
    //           << measurePhysR[1] << ", "
    //           << measurePhysR[2] << ") (PU): "
    //           << T{globalMeasurements[0]}<< std::endl;

        
    //if ( iT%iTvtk == 0 ) {getGraphicalResults(sLattice, converter, iT, superGeometry, amplitude);}





    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    if ( iT%iTtimer == 0 ) {timer.update(iT); timer.printStep();}
    }
 
   timer.stop();
   timer.printSummary();
 }
 
 
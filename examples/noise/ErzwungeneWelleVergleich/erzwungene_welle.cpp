
/*Notiz: Erzwungene Welle
Das Programm kompiliert Fehlerfrei und kann durchgeführt werden. 
Bei dem ausführen des Programms muss mit --iTmax XXX eine Zahl angegeben werden, wieviele Iterationsdurchläufe das Programm durchläuft
Terminalbefehl: make clean; make; ./Freie_Welle --iTmax 150 --peakN 2 (Mit welchem Wellenberg wird die Geschwindigkeit berechnet?) --outdir tmp_periodic_05  --> Zusätzlich kann die Amplitude angegeben werden a--
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
#include <array>   // NEU, für std::array
#include <cmath>   // falls nicht schon vorhanden (atan2, cos, sin, sqrt)


using namespace olb;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = descriptors::D3Q19<>;
using BulkDynamics   = BGKdynamics<T, DESCRIPTOR>;
using SpongeDynamics = SpongeLayerDynamics<T, DESCRIPTOR, momenta::BulkTuple, equilibria::SecondOrder>;

const int ndim = 3; // a few things (e.g. SuperSum3D) cannot be adapted to 2D, but this should help speed it up
const T lambda_phys = T(0.6);         // Lambda verändern!
const T physDeltaX        = 0.02;     // grid spacing [m]
const T physLength        = 1.;       // length of the cuboid [m]
const T physspan          = 0.46;
const T physwidth         = 0.46;
const T physLidVelocity   = 1.0;      // velocity imposed on lid [m/s] Fuer die Machzahl relevant (Vorher 1.0, jetzt ein zehntel der Schallgeschwindigkeit)
const T physViscosity     = 1.5e-3;   // kinetic viscosity of fluid [m*m/s] Fuer die Relaxationszeit verantwortlich    
const T physDensity       = 1.;       // fluid density of air (20°C)[kg/(m*m*m)]
const T physMaxT          = 0.5;      // maximal simulation time [s]
const T charL   = 1;      // z.B. deine Wellenlänge
const int res   = 120;               // ~30–40 Zellen pro λ
const T Ma      = 0.01;             // kleine Machzahl
const T charV   = 0.003;              // charakteristische phys. Geschwindigkeit (z.B. U' = p'/(rho*c))
const T rho0    = 1.0;
const T nu_phys = 1.5e-5;
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
  std::cout << "dx: "<<dx;
   switch (boundarytype) {
   // eternal and damping: 3 is the actual fluid; periodic: 1 is the fluid
   case periodic:{
     superGeometry.rename(2, 1);
     break;
   } 
     
   case local: {
    
    //superGeometry.rename(1,3,1,pointSource);

      // --- SPONGE REGION anlegen: äußere Hülle wird Material=4 Sponge Layer fuer den Fall local-----------------
    
    {
      // Dicke der Hülle in Zellen:
      const int spongeCells = 12;                       // z.B. 8 Zellen
      const T   dx          = converter.getPhysDeltaX();
      const T   sx          = spongeCells*dx;
      const T   sy          = spongeCells*dx;
      const T   sz          = spongeCells*dx;

      // Gesamtgröße (m): 2.4 x physwidth x physspan, Ursprung wie oben gewählt
      const T Lx = 2.4/physLength, Ly = physwidth/physLength, Lz = physspan/physLength;

      // INNERER „fluid“-Kern als Cuboid (ohne Sponge-Hülle)
      Vector<T,3> innerExtend(Lx - 2*sx, Ly - 2*sy, Lz - 2*sz);
      Vector<T,3> innerOrigin(- (Lx - 2*sx)/2., - (Ly - 2*sy)/2., - (Lz - 2*sz)/2.);
      IndicatorCuboid3D<T> inner(innerExtend/physLength, innerOrigin/physLength);

      // Alles was aktuell 1 ist -> 4 (erst mal komplette Hülle)
      superGeometry.rename(1, 4);

      // Im inneren Bereich wieder 4 -> 1 (nur der Kern bleibt "echtes" Fluid)
      superGeometry.rename(4, 1, inner);
    }
    superGeometry.rename(2, 1);
    Vector<T,3> center(-0.0, 0., 0.); // Zentrum der Welle
    T radius = dx;
    IndicatorSphere3D<T> pointSource(center, radius);
    superGeometry.rename(1,3,pointSource);
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
      // Fluid (Kern + evtl. Quell-Material 3)
      auto fluid  = superGeometry.getMaterialIndicator({1,3});
      auto sponge = superGeometry.getMaterialIndicator({4});

      sLattice.defineDynamics<BulkDynamics>(fluid);
      sLattice.defineDynamics<SpongeDynamics>(sponge);

      // Relaxation wie gehabt
      const T omega = converter.getLatticeRelaxationFrequency();
      sLattice.setParameter<descriptors::OMEGA>(omega);

      // (Optional) Sponge-Stärke gleichmäßig setzen (klein anfangen!)
      #ifdef descriptors_SPONGE_STRENGTH_EXISTS
      {
        const T sigma = 0.2; // 0..1: 0=kein Dämpfen, 1=stark; starte konservativ
        AnalyticalConst3D<T,T> sigmaF(sigma);
        sLattice.defineParameter<descriptors::SPONGE_STRENGTH>(sponge, sigmaF);
      }
      #endif

      // Falls dein OpenLB die SPONGE_STRENGTH-Parameter-Flag nicht kennt,
      // einfach den #ifdef-Block ignorieren – der SpongeDynamics dämpft
      // dann mit seinem Default-Profil.


       break;
     }
 
     sLattice.setParameter<descriptors::OMEGA>(omega);
   
     // Make the lattice ready for simulation
     sLattice.initialize();
   
     clout << "Prepare Lattice ... OK" << std::endl;
   } //prepareLattice
 
 void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
   SuperLattice<T, DESCRIPTOR>& sLattice,
   std::size_t iT, SuperGeometry<T,ndim>& superGeometry,BoundaryType boundarytype, T amplitude, T rho0, T lambda_phys)
  {  
 
    if (boundarytype == local) {
      auto domain = superGeometry.getMaterialIndicator({3});
      T wellenzahl=2. * std::numbers::pi_v<T>/lambda_phys ;//k=2pi/lamda
      T kreisfrequenz=343.46*wellenzahl;
      T phase =0.;
      //T time = converter.getPhysTime(iT);  // physikalische Zeit aus Lattice-Zeit
      T time= converter.getPhysTime(iT);
      T cs=sqrt(T(1)/descriptors::invCs2<T,DESCRIPTOR>());
      //T envelope = std::sin(std::min(1.0, iT / 40.0) * std::numbers::pi_v<T> / 2.0);
      int dumb = 2;
      olb::SchallwelleRho<3, T, DESCRIPTOR> schallquelle(rho0, amplitude, wellenzahl, kreisfrequenz, phase, time, dumb,converter);
      olb::SchallwelleGesch<3,T, DESCRIPTOR> schallquelle_geschwindigkeit(amplitude,wellenzahl,kreisfrequenz,phase,time,rho0,cs,converter);
      //olb::SchallwelleGesch<3,T, DESCRIPTOR> schallquelle_geschwindigkeit(amplitude,wellenzahl,kreisfrequenz,phase,time,rho0,cs,boundarytype, converter);
      // AnalyticalConst3D<T,T> rhoF( schallquelle);
      //AnalyticalConst3D<T,T> schallquelle_geschwindigkeit(0, 0,0);  
 


      // sLattice.defineRhoU(domain, schallquelle, uInf);
      // sLattice.iniEquilibrium(domain, schallquelle, uInf);
      if(iT==0)
       {
         AnalyticalConst3D<T,T> schallquelle(1.);       // Definiert rho auf 1
         AnalyticalConst3D<T,T> schallquelle_geschwindigkeit(0., 0.,0.);   // Definiert die Geschwindigkeit überall auf 0 in X und Y Richtung
       }
      sLattice.defineRhoU(domain, schallquelle, schallquelle_geschwindigkeit);
      sLattice.iniEquilibrium(domain, schallquelle, schallquelle_geschwindigkeit);
    }
     
 
   if (boundarytype == periodic && iT==0) {
     auto domain = superGeometry.getMaterialIndicator({1});
    
     
     T wellenzahl=2. * std::numbers::pi_v<T>/lambda_phys ;//k=2pi/lamda
     T kreisfrequenz=343.46*wellenzahl;
     T phase =0.;
     //T time = converter.getPhysTime(iT);  // physikalische Zeit aus Lattice-Zeit
     T time= converter.getPhysTime(iT);
     T cs=sqrt(T(1)/descriptors::invCs2<T,DESCRIPTOR>());
     T dumb=1.;

     olb::SchallwelleRho<3, T, DESCRIPTOR> schallquelle(rho0, amplitude, wellenzahl, kreisfrequenz, phase, time, dumb, converter);
     olb::SchallwelleGesch<3,T, DESCRIPTOR> schallquelle_geschwindigkeit(amplitude,wellenzahl,kreisfrequenz,phase,time,rho0,cs,converter);
     // AnalyticalConst3D<T,T> rhoF( schallquelle);

 
 
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
 
     //output pressure image
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
     T                          pmin(converter.getPhysPressure(-amplitude / 20));
     T                          pmax(converter.getPhysPressure(+amplitude / 20));
     linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_hline_" + ss.str(), "pressure [PU]",
                       horizontal, false, false, pmin, pmax);  // TODO setRange=true (before pmin, pmax)
     //linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_vline_" + ss.str(), "pressure [PU]", vertical,
     //                    false, true, pmin, pmax);
     // linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_diagonal_" + ss.str(), "pressure [PU]",
     //                     diagonal2d, false, false, pmin, pmax);
 }  // getGraphicalResults
 
 
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
  // Welche Spitze auswerten? 1=erster Wellenberg, 2=zweiter, ...
  int peakN = args.getValueOrFallback("--peakN", 1);

  
  OstreamManager clout( std::cout,"main" ); // writing all output first in a userdefined Buffer of type OMBuf. On a flush it spits out at first the userdefined text in squared brackets and afterwards everything from the buffer

  // Provide the unit converter the characteristic entities
  UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> converter(
    (size_t)res,             // resolution = Nx auf charL
    Ma/std::sqrt(3.),        // charLatticeVelocity
    charL,                   // charPhysLength
    charV,                   // charPhysVelocity
    nu_phys,                 // physViscosity
    rho0                     // physDensity
    );
    converter.print();


  // --- Wellenzahl in physikalischen und lattice Einheiten ---
  // Nutze dieselbe λ bzw. wellenzahl wie in setBoundaryValues (hier: λ_phys = 0.5 m)
  
            
  const T k_phys = 2.*std::numbers::pi_v<T> / lambda_phys;         // [rad/m]
  const T k_lat  = k_phys * converter.getPhysDeltaX(); // [rad per lattice cell]
  const T k2_lat = k_lat * k_lat;

  // theoretisches c_s (lattice und physisch)
  const T cs_lat  = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());          // ≈ 1/√3
  const T cs_phys = (converter.getPhysDeltaX()/converter.getPhysDeltaT()) * cs_lat;  // nur Info







  // === 2nd Step: Prepare Geometry ===
  BoundaryType boundarytype = local;
  Vector<T,ndim> originFluid(-2.4/physLength/2., -physwidth/physLength/2., -physspan/physLength/2.);
  Vector<T,ndim> extendFluid(2.4/physLength, physwidth/physLength, physspan/physLength);
  IndicatorCuboid3D<T> domainFluid(extendFluid, originFluid);
  // -----------Variabeln definiere Messungen
  size_t nplot                  = args.getValueOrFallback( "--nplot",             100 );  
  size_t iTout                  = args.getValueOrFallback( "--iTout",             0   );  
    
  //----------------------------- Geometrie aufspannen
  Vector<T,ndim> extend{2.47/physLength, physwidth/physLength, physspan/physLength};
  Vector<T,ndim> origin{-2.4/physLength/2., -physwidth/physLength/2., -physspan/physLength/2.};
  IndicatorCuboid3D<T> cuboid(extend, origin);
  CuboidDecomposition3D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());
  cuboidDecomposition.setPeriodicity({false,false,false});
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
  // ---------------------------------------Zwischenschritt: Messwerte nehmen //-----------Vorgegeben---
  // #if defined(FEATURE_TWOD)
  // Vector<T,ndim> measurePhysR{-0.0875,0.0052};
  // Vector<int,3> measureLatticeR{};
  //#elif defined(FEATURE_THREED)
  // --- Messpunkte in physikalischen Koordinaten (m)
  // --- Zwei Messpunkte in physikalischen Koordinaten (m)
  std::array<Vector<T,ndim>,2> measurePhysR = {
    Vector<T,ndim>{0.546, 0.0, 0.230},
    Vector<T,ndim>{0.586, 0.0, 0.230}
  };
  std::array<Vector<int,4>,2> measureLatticeR{};

  SuperD<T,descriptors::D3<fields::PHYS_R,descriptors::SCALAR>> watchpointsD(loadBalancer);

  for (int k=0; k<2; ++k) {
    if (auto latticeR = cuboidDecomposition.getLatticeR(measurePhysR[k])) {
      measureLatticeR[k] = *latticeR;

      if (loadBalancer.isLocal(measureLatticeR[k][0])) {
        auto& blockD = watchpointsD.getBlock(measureLatticeR[k][0]);
        // Platz schaffen: 1 Zelle je Watchpoint am Blockende
        blockD.resize({blockD.getNx()+1,1,1});
        auto watchpoint = blockD.get(blockD.getNcells()-1);
        watchpoint.template setField<fields::PHYS_R>(measurePhysR[k]);

        // Lokalen Zellenindex merken (wir nutzen measureLatticeR[k][1] dafür)
        measureLatticeR[k][1] = blockD.getNcells()-1;
      }
    } else if (singleton::mpi().getRank()==0) {
      std::cout << "[WARNUNG] Messpunkt " << k << " wurde NICHT lokal gefunden!\n";
    }
  }

  watchpointsD.setProcessingContext(ProcessingContext::Simulation);

  // Kopplung: schreibt den berechneten Druck in die Watchpoints
  SuperLatticePointCoupling pressureO(PressureO{},
                                      names::NavierStokes{}, sLattice,
                                      names::Points{}, watchpointsD);

  // Zeitreihen und Konstanten für spätere Auswertung
  std::vector<T> p1, p2; 
  p1.reserve(iTmax); 
  p2.reserve(iTmax);

  const T dtPhys = converter.getPhysDeltaT();
  // oben sicherstellen: #include <cmath>
  const T dx = std::sqrt(
    (measurePhysR[1][0] - measurePhysR[0][0]) * (measurePhysR[1][0] - measurePhysR[0][0]) +
    (measurePhysR[1][1] - measurePhysR[0][1]) * (measurePhysR[1][1] - measurePhysR[0][1]) +
    (measurePhysR[1][2] - measurePhysR[0][2]) * (measurePhysR[1][2] - measurePhysR[0][2])
  );


  // Für die optionale Phasenmethode: physikalische Kreisfrequenz der Anregung bestimmen
  T omegaPerStep = T(0);
  if (boundarytype == local) {
    // in setBoundaryValues(local) wurde sin(iT * 2π/40) verwendet
    omegaPerStep =  2. * std::numbers::pi_v<T> /40.0;
  } else {
    // periodic-Zweig bei dir nutzt 2π * 2 / 40 (zwei Perioden in 40 Schritten)
    omegaPerStep =  2. * std::numbers::pi_v<T>*2. /40.0;
  }
  const T omegaPhys = omegaPerStep / dtPhys;

  CSV<T> csvWriter("Welle", ';', {"iT", "t", "p1", "p2"}, ".csv");
  CSV<T> csvSummary("cp_vs_k", ';',
    {"k_lat", "k2_lat", "cp_lat_xcorr", "cp_lat_phase", "cs_lat"},
    ".csv");

  // Messung Plot Amplitudenverlauf
  // --- Druck-Functor einmal definieren ---
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressureF(sLattice, converter);
  AnalyticalFfromSuperF3D<T> pressureInterp(pressureF, true, true);

  // --- Abtastlinie anlegen (x von links nach rechts, y=z=0) ---
  const int NxLine = 200; // Auflösung für den Plot
  std::vector<T> x_phys(NxLine), p_max(NxLine, T(0)), p_rss(NxLine, T(0));
  std::vector<T> p_snap(NxLine, T(0)); // Snapshot-Werte
  std::size_t n_accum = 0;

  const T x0 = -physLength*2.4/2.;       // wie in deiner Geometrie
  const T x1 =  physLength*2.4/2.;
  const T y0 =  T(0), z0 = T(0);
  for (int i=0; i<NxLine; ++i) {
    x_phys[i] = x0 + (x1 - x0) * ( (i + T(0.5)) / T(NxLine) );
  }

  //=== Uebergabe der BoundaryConditions
  if (boundarytype == periodic) {setBoundaryValues(converter, sLattice, 0, superGeometry, boundarytype, amplitude, rho0,lambda_phys);}
  //--------------------------------------- FOR SCHLEIFE-------------------------------------------------------------------------------------------
  for (std::size_t iT=0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    if (boundarytype == local) {setBoundaryValues(converter, sLattice, iT, superGeometry, boundarytype, amplitude,rho0,lambda_phys);}
    
    // ------------------------- Messwerte nehmen

    // ------------------------- Messwerte nehmen (2 Sensoren)
    pressureO.execute();
    watchpointsD.setProcessingContext(ProcessingContext::Evaluation);

    // lokale Messwerte sammeln
    std::vector<T> localP(2, T(0)), globalP(2, T(0));
    for (int k=0; k<2; ++k) {
      if (loadBalancer.isLocal(measureLatticeR[k][0])) {
        auto& blk = watchpointsD.getBlock(loadBalancer.loc(measureLatticeR[k][0]));
        auto cell = blk.get(measureLatticeR[k][1]);
        const T pu = cell.template getField<descriptors::SCALAR>(); // Lattice-Druck (PU)
        localP[k] += converter.getPhysPressure(pu);                 // in phys. Druck [Pa]
      }
    }

    #ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceVect(localP, globalP, MPI_SUM);
    singleton::mpi().bCast(globalP.data(), globalP.size());
    #else
    globalP = localP;
    #endif

    // Zeitreihen füllen und CSV schreiben
    p1.push_back(globalP[0]);
    p2.push_back(globalP[1]);
    csvWriter.writeDataFile(iT, {converter.getPhysTime(iT), globalP[0], globalP[1]});


    // std::cout << "[iT=" << iT << ", t=" << converter.getPhysTime(iT) << "s] "
    //           << "Gemessener Druck am Punkt ("
    //           << measurePhysR[0] << ", "
    //           << measurePhysR[1] << ", "
    //           << measurePhysR[2] << ") (PU): "
    //           << T{globalMeasurements[0]}<< std::endl;

        
    if ( iT%iTvtk == 0 ) {getGraphicalResults(sLattice, converter, iT, superGeometry, amplitude);}


    //===Zwischenschritt Amplitudenverlauf=====
    // --- p(x,t) auf der Linie auslesen ---
      for (int i=0; i<NxLine; ++i) {
        T out;
        T pos[3] = { x_phys[i], y0, z0 };
        pressureInterp(&out, pos); // out in [Pa], weil *Phys*Pressure
        const T a = std::abs(out);
        if (a > p_max[i]) p_max[i] = a;      // Peak-Hüllkurve
        p_rss[i] += out*out;                 // für RMS
      }
      // Optional: Snapshot zu einem gewünschten Zeitpunkt sichern
      // Beispiel: Snapshot beim Maximum am ersten Sensor (wenn du t* kennst):
      // if (iT == iT_snapshot) { p_snap = momentane Werte; }
      ++n_accum;


    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    if ( iT%iTtimer == 0 ) {timer.update(iT); timer.printStep();}
    }

    #include <limits>
#include <algorithm>

// --- Hilfsfunktion: alle Peak-Zeiten (lokale Maxima) finden, mit Parabel-Refinement
auto findPeakTimes = [&](const std::vector<T>& p, T dt, int guardSamples, T minAmp){
  std::vector<T> peaks;
  const int N = (int)p.size();
  const int iStart = std::min(std::max(guardSamples, 1), N-3);

  for (int i = iStart+1; i < N-1; ++i) {
    if (p[i] > minAmp && p[i] > p[i-1] && p[i] > p[i+1]) {
      // Parabolische Subsample-Interpolation um i
      const T a = p[i-1] - 2*p[i] + p[i+1];
      const T b = p[i-1] - p[i+1];
      T delta = T(0);
      if (std::abs(a) > T(1e-30)) {
        delta = b / (2*a);                 // ideal in [-0.5, 0.5]
        delta = std::clamp(delta, T(-0.5), T(0.5));
      }
      peaks.push_back( (i + delta) * dt ); // physische Zeit [s]
    }
  }
  return peaks;
};

// --- Parameter für die Peak-Suche
const int guard = 5; // ein paar erste Samples ignorieren
T estAmp = T(0);
for (int i = guard; i < (int)std::min<std::size_t>(p1.size(), guard+50); ++i) {
  estAmp = std::max(estAmp, std::abs(p1[i]));
}
const T minAmp = estAmp * T(0.2); // 20% der frühen Spitze als Schwellwert

// --- Alle frühen Peaks beider Sensoren
auto peaks1 = findPeakTimes(p1, dtPhys, guard, minAmp);
auto peaks2 = findPeakTimes(p2, dtPhys, guard, minAmp);

// --- gewünschten Peak wählen (1=erster, 2=zweiter, ...)
//     (CLI: --peakN 2 für zweiten Wellenberg)
const int n = std::max(1, peakN);
if ((int)peaks1.size() >= n && (int)peaks2.size() >= n) {
  const T t1 = peaks1[n-1];
  const T t2 = peaks2[n-1];

  if (t2 != t1 && std::isfinite(t1) && std::isfinite(t2)) {
    const T cp_peak_phys = dx / std::abs(t2 - t1); // [m/s]
    const T cp_peak_lat  = cp_peak_phys * converter.getPhysDeltaT() / converter.getPhysDeltaX();
    const T cs_lat_here  = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());
    const T ratio        = cp_peak_lat / cs_lat_here;

    if (singleton::mpi().getRank()==0) {
      std::cout << "[cp|PEAK#" << n << "] t1="<<t1<<" s, t2="<<t2<<" s"
                << " -> c_p="<<cp_peak_phys<<" m/s"
                << " | c_p_lat="<<cp_peak_lat
                << " | c_p/c_s="<<ratio << "\n";
    }

    // CSV schreiben (eigene Datei oder an deine bestehende anhängen)
    
    CSV<T> csvPeak("cp_peak_n", ';',
      {"k_lat","k2_lat","peakN","cp_phys","cp_lat","cs_lat","cp_over_cs"}, ".csv");
    csvPeak.writeDataFile(0,{k_lat, k2_lat, peakN, cp_peak_phys, cp_peak_lat, cs_lat_here, ratio});
    
  } else {
    if (singleton::mpi().getRank()==0) std::cout << "[cp|PEAK#" << n << "] ungültige Peak-Zeiten.\n";
  }
} else {
  if (singleton::mpi().getRank()==0) {
    std::cout << "[cp|PEAK#" << n << "] Nicht genug Peaks gefunden: "
              << "sensor1="<<peaks1.size()<<", sensor2="<<peaks2.size()<<"\n";
  }
}







    // ------------------------- Auswertung: Cross-Correlation & Phasenmethode

    // Optional: Einschwingtransienten verwerfen (z.B. erste 1-2 Perioden)
    auto drop_front = [&](std::vector<T>& v, std::size_t n){
      if (v.size()>n) v.erase(v.begin(), v.begin()+n);
    };
    {
      // 1 Periode ≈ (2π / omegaPhys) Sekunden
      const T Tper = 2.*std::numbers::pi_v<T> / std::max(omegaPhys, T(1e-12));
      const std::size_t Ndrop = (std::size_t)std::ceil(1.0 * Tper / dtPhys); // 1 Periode
      drop_front(p1, Ndrop);
      drop_front(p2, Ndrop);
    }

    // Cross-Correlation (einfach, normalisiert, Lag um 0 herum suchen)
    auto xcorrLag = [&](const std::vector<T>& a, const std::vector<T>& b)->int {
      const int N = (int)std::min(a.size(), b.size());
      if (N<=3) return 0;
      // maximaler Lag heuristisch begrenzen
      const int maxLag = std::min( (int)std::round(0.5 * dx / std::max(dtPhys, T(1e-12))), N-1 );

      // Mittelwerte entfernen
      T ma=0, mb=0; 
      for(int i=0;i<N;++i){ ma+=a[i]; mb+=b[i]; } 
      ma/=N; mb/=N;

      T bestC = -1e300; 
      int bestLag = 0;
      for (int lag=-maxLag; lag<=maxLag; ++lag) {
        T num=0, da=0, db=0;
        for (int i=0;i<N;++i) {
          const int j = i+lag;
          if (j<0 || j>=N) continue;
          const T aa = a[i]-ma;
          const T bb = b[j]-mb;
          num += aa*bb; da += aa*aa; db += bb*bb;
        }
        if (da>0 && db>0) {
          const T c = num / std::sqrt(da*db);
          if (c>bestC) { bestC=c; bestLag=lag; }
        }
      }
      return bestLag;
    };

    int lag = xcorrLag(p1,p2);
    T dt = lag * dtPhys;
    T cp_xcorr = dx / std::max(std::abs(dt), T(1e-12));

    if (singleton::mpi().getRank()==0) {
      std::cout << "[cp|XCORR] dx="<<dx<<" m, lag="<<lag<<" Samples, dt="<<dt
                <<" s -> c_p="<<cp_xcorr<<" m/s\n";
    }

    // -------- Optional: Phasenmethode (benötigt korrekte omegaPhys) ----------------------------------------------------------------------------------------------------------------------------
    auto complexProj = [&](const std::vector<T>& p)->std::pair<T,T>{
      T A=0, B=0; // Re=A, Im=B (mit -sin für Im)
      const std::size_t N = p.size();
      for (std::size_t n=0;n<N;++n){
        const T t = n*dtPhys;
        A += p[n]*std::cos(omegaPhys*t);
        B += p[n]*std::sin(omegaPhys*t);
      }
      const T phase = std::atan2(-B, A); // Phase ∈ (-π, π]
      return {A, phase};
    };

    if (p1.size()>=8 && p2.size()>=8) {
      auto [A1,phi1] = complexProj(p1);
      auto [A2,phi2] = complexProj(p2);

      T dphi = phi2 - phi1;
      while (dphi >  M_PI) dphi -= 2*M_PI;
      while (dphi < -M_PI) dphi += 2*M_PI;

      const T cp_phase = std::abs(omegaPhys * dx / std::max(std::abs(dphi), T(1e-12)));

      if (singleton::mpi().getRank()==0) {
        std::cout << "[cp|PHASE] dphi="<<dphi<<" rad, omega="<<omegaPhys
                  <<" rad/s -> c_p="<<cp_phase<<" m/s\n";
      }

        // Umrechnung nach lattice-Einheiten:
            const T cp_lat_xcorr = cp_xcorr * converter.getPhysDeltaT() / converter.getPhysDeltaX();
            const T cp_lat_phase = cp_phase * converter.getPhysDeltaT() / converter.getPhysDeltaX();
            csvSummary.writeDataFile(0, {k_lat, k2_lat, cp_lat_xcorr, cp_lat_phase, cs_lat});

    }
    
    //===Amplitudenverlauf eintragen=======
    // RMS aus Summe der Quadrate
    std::vector<T> p_rms(NxLine);
    for (int i=0; i<NxLine; ++i) {
      p_rms[i] = std::sqrt(p_rss[i] / std::max<std::size_t>(n_accum,1));
    }

    // --- CSVs schreiben ---
    // 1) Peak-Hüllkurve
    CSV<T> csvAmpPeak("amplitude_peak_vs_x", ';', {"x_phys_m", "p_peak_Pa"}, ".csv");
    for (int i=0; i<NxLine; ++i) csvAmpPeak.writeDataFile(i, {x_phys[i], p_max[i]});

    // 2) RMS-Hüllkurve
    CSV<T> csvAmpRms("amplitude_rms_vs_x", ';', {"x_phys_m", "p_rms_Pa"}, ".csv");
    for (int i=0; i<NxLine; ++i) csvAmpRms.writeDataFile(i, {x_phys[i], p_rms[i]});

    // 3) (optional) Snapshot
    CSV<T> csvSnap("amplitude_snapshot_vs_x", ';', {"x_phys_m", "p_snapshot_Pa"}, ".csv");
    for (int i=0; i<NxLine; ++i) csvSnap.writeDataFile(i, {x_phys[i], p_snap[i]});



   timer.stop();
   timer.printSummary();
 }
 
 
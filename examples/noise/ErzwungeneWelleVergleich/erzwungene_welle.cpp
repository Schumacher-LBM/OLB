
/* Notiz: Erzwungene Welle
 * Das Programm kompiliert Fehlerfrei und kann durchgeführt werden. 
 * Bei dem ausführen des Programms muss mit --iTmax XXX eine Zahl angegeben werden, wieviele Iterationsdurchläufe das Programm durchläuft
 * Terminalbefehl: make clean; make; ./erzwungenewelle&>log.txt --iTmax 150 --peakN 2 (Mit welchem Wellenberg wird die Geschwindigkeit berechnet?) --outdir tmp_periodic_05  --> Zusätzlich kann die Amplitude angegeben werden a--
 * To Do: Schallgeschw. und Amplitude an die Werte von Luft anpassen und Quellen finden. Außerdem die Größe des Mediums angeben. Messwerte nehmen
 */
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
 #include <limits>
 #include <algorithm>
 
 using namespace olb;
 
 using T = FLOATING_POINT_TYPE;
 using DESCRIPTOR = descriptors::D3Q19<>;
 using BulkDynamics   = BGKdynamics<T, DESCRIPTOR>;
 using SpongeDynamics = SpongeLayerDynamics<T, DESCRIPTOR, momenta::BulkTuple, equilibria::SecondOrder>;
 
 const int ndim = 3; // a few things (e.g. SuperSum3D) cannot be adapted to 2D, but this should help speed it up
 const T lambda_phys       = T(0.38);         // Lambda verändern!
 const int lambda_lat      = 10.;            //Wieviele Zellen sollen für eine Wellenlänge genutzt werden
 const int nWaves          = 6.;
 // const T physDeltaX        = 0.02;     // grid spacing [m]
 const T physLength        = 1.;       // length of the cuboid [m]
 const T Domainelength     = 10.*lambda_phys;
 const T physspan          = 10.*lambda_phys;
 const T physwidth         = 10*lambda_phys;
 const T physLidVelocity   = 1.0;      // velocity imposed on lid [m/s] Fuer die Machzahl relevant (Vorher 1.0, jetzt ein zehntel der Schallgeschwindigkeit)
 const T physViscosity     = 1.5e-2;   // kinetic viscosity of fluid [m*m/s] Fuer die Relaxationszeit verantwortlich    
 const T physDensity       = 1.;       // fluid density of air (20°C)[kg/(m*m*m)]
 const T physMaxT          = 0.5;        // maximal simulation time [s]
 const T physDeltaT        = 0.00078125;// Messung 1: 0.00078125;//((0.68255-0.5)/3)/physViscosity*physDeltaX*physDeltaX;// 0,68255, weil Tau 0,68255 sein soll. Vorher: physDeltaX/343.46;  // temporal spacing [s] t=physDeltaX/c_s (Vorher 0.00078125, Jetzt: 5,8e-5)
 // Alte Werte
 // const T charL   = 1;      // z.B. deine Wellenlänge
 // const int res   = 90;               // ~30–40 Zellen pro λ
 // const T Ma      = 0.01;             // kleine Machzahl
 // const T charV   = 0.003;              // charakteristische phys. Geschwindigkeit (z.B. U' = p'/(rho*c))
 // const T rho0    = 1.0;
 // const T nu_phys = 1.5e-5;
 
 
 
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
                       IndicatorF3D<T>& domainFluid)
  {
    OstreamManager clout(std::cout, "prepareGeometry");
    clout << std::endl << "Prepare Geometry ..." << std::endl;
  
    // all nodes to temporary type
    superGeometry.rename(0, 2);
  
      
     
     //superGeometry.rename(1,3,1,pointSource);
 
       // --- SPONGE REGION anlegen: äußere Hülle wird Material=4 Sponge Layer fuer den Fall local-----------------
     
       // Dicke der Hülle in Zellen:
       const int spongeCells = 2;                       // z.B. 8 Zellen
       const T   dx          = converter.getPhysDeltaX();
       std::cout << "dx: "<<dx;
       const T   sx          = spongeCells*dx;
       const T   sy          = spongeCells*dx;
       const T   sz          = spongeCells*dx;
 
       // Gesamtgröße (m): 2.4 x physwidth x physspan, Ursprung wie oben gewählt
       const T Lx = 2.4/physLength, Ly = physwidth/physLength, Lz = physspan/physLength;
 
       // INNERER „fluid“-Kern als Cuboid (ohne Sponge-Hülle)
       Vector<T,3> innerExtend(Lx - 2*sx, Ly - 2*sy, Lz - 2*sz);
       Vector<T,3> innerOrigin(- (Lx - 2*sx)/2., - (Ly - 2*sy)/2., - (Lz - 2*sz)/2.);
       clout << "innerExtend=[" << innerExtend[0] << '\n';
       IndicatorCuboid3D<T> inner(innerExtend/physLength, innerOrigin/physLength);
 
 
       // Im inneren Bereich wieder 4 -> 1 (nur der Kern bleibt "echtes" Fluid)
       superGeometry.rename(2, 1, inner);
       // Alles was nicht 1 war (also noch 2 ist) -> 4 (sollte Hülle sein)
       superGeometry.rename(2, 4);
       
     Vector<T,3> center(-0.0, 0., 0.); // Zentrum der Welle
     T radius = 2.*dx;
     IndicatorSphere3D<T> pointSource(center, radius);
     superGeometry.rename(1,3,pointSource);
  
    superGeometry.getStatistics().print();
 
  }//prepareGeometry
 
 
    // Set up the geometry of the simulation
    void prepareLattice(UnitConverter<T, DESCRIPTOR> const& converter, SuperLattice<T, DESCRIPTOR>& sLattice,
                        SuperGeometry<T, ndim>& superGeometry, T rho0, T u0, T amplitude, T alpha,
                        T dampingDepthPU, T lengthDomain, T dampingStrength)
    {
      OstreamManager clout(std::cout, "prepareLattice");
      clout << std::endl << "Prepare Lattice ..." << std::endl;
    
      const T omega = converter.getLatticeRelaxationFrequency();
 
 
 // In setBoundaryValues(): kreisfrequenz = omegaPhys;
 // Bei complexProj/Phasenmethode: ebenfalls omegaPhys benutzen.
 
      // Material=3 --> bulk dynamics
      auto bulkIndicator = superGeometry.getMaterialIndicator(
          {1,3}); // for local bcs all around, corners remain at 2, so they are included here
      sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
    
        //boundary::set<boundary::LocalPressure>(sLattice, superGeometry, 3);
       // Fluid (Kern + evtl. Quell-Material 3)
       auto fluid  = superGeometry.getMaterialIndicator({1,3});
       auto sponge = superGeometry.getMaterialIndicator({4});
 
       sLattice.defineDynamics<BulkDynamics>(fluid);
       sLattice.defineDynamics<SpongeDynamics>(sponge);
 
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
 
 
  
       // Relaxation wie gehabt
      sLattice.setParameter<descriptors::OMEGA>(omega);
    
      // Make the lattice ready for simulation
      sLattice.initialize();
    
      clout << "Prepare Lattice ... OK" << std::endl;
    } //prepareLattice
  
  void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
    SuperLattice<T, DESCRIPTOR>& sLattice,
    std::size_t iT, SuperGeometry<T,ndim>& superGeometry, T amplitude, T rho0, T lambda_phys,T omegaPhys)
   {  
  
       auto sourceIndicator = superGeometry.getMaterialIndicator({3});
 
       const T cs_lat = std::sqrt(T(1)/descriptors::invCs2<T,DESCRIPTOR>());
       const T c_phys = (converter.getPhysDeltaX()/converter.getPhysDeltaT()) * cs_lat;
       
 
       T wellenzahl=2. * std::numbers::pi_v<T>/lambda_phys ;//k=2pi/lamda
       //const T kreisfrequenz = c_phys * wellenzahl;
       T phase =0.;
       const T kreisfrequenz = omegaPhys;
 
 
       T time= converter.getPhysTime(iT);  // physikalische Zeit aus Lattice-Zeit
 
       // === ALT: schallquelle variabel über den Raum und angepasste Geschwindigkeit aus freier Welle
       // olb::SchallwelleRho<3, T, DESCRIPTOR> schallquelle(rho0, amplitude, wellenzahl, kreisfrequenz, phase, time, dumb,converter);
       // olb::SchallwelleGesch<3,T, DESCRIPTOR> schallquelle_geschwindigkeit(amplitude,wellenzahl,kreisfrequenz,phase,time,rho0,cs,converter);
       // === ALT
 
       if (iT==0) {
         // einmalig initialisieren
         auto fluid = superGeometry.getMaterialIndicator({1,3,4}); // oder {1,3}
         AnalyticalConst3D<T,T> u0(0,0,0);
         AnalyticalConst3D<T,T> rho1230(1.);
         sLattice.defineRhoU(fluid, rho1230, u0);
         sLattice.iniEquilibrium(fluid, rho1230, u0);
       } else {
 
         T phi =  - kreisfrequenz * time + phase; //+waveNumber * x
         
         // z.B. 1 Periodendauer rampen
         T Tramp = (2*std::numbers::pi_v<T>)/kreisfrequenz;  // Rampdauer in s
         T envelope = std::sin(std::min(T(1), time/Tramp) * std::numbers::pi_v<T>/T(2));
 
         T p_phys = envelope*amplitude * std::sin(phi);  // [Pa]
         T rhoT   = converter.getLatticeDensityFromPhysPressure(p_phys);
 
         // danach in jedem Schritt KEIN iniEquilibrium mehr!
         AnalyticalConst3D<T,T> rhoF( rhoT );
         sLattice.defineRho( sourceIndicator, rhoF );
       }
       
       
       // sLattice.defineRhoU(domain, schallquelle, schallquelle_geschwindigkeit);
       // sLattice.iniEquilibrium(domain, schallquelle, schallquelle_geschwindigkeit);
 
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
      T                          ndatapoints =500.; //converter.getResolution(); // number of data points on line
      AnalyticalFfromSuperF3D<T> pressure_interpolation(pressure, true, true);
      T                          pmin(converter.getPhysPressure(-amplitude / 20));
      T                          pmax(converter.getPhysPressure(+amplitude / 20));
      linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_hline_" + ss.str(), "pressure [PU]",
                        horizontal, false, false, pmin, pmax);  // TODO setRange=true (before pmin, pmax)
     //  linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_vline_" + ss.str(), "pressure [PU]", vertical,
     //                     false, true, pmin, pmax);
      // linePlot<ndim, T>(pressure_interpolation, ndatapoints, dist, "pressure_diagonal_" + ss.str(), "pressure [PU]",
      //                     diagonal2d, false, false, pmin, pmax);
  }  // getGraphicalResults
  
  
 int main(int argc, char* argv[])
 {
   // === 1st Step: Initialization ===
   initialize(&argc, &argv);
   CLIreader args(argc, argv);
   std::string outdir = args.getValueOrFallback<std::string>("--outdir", "");
   outdir += "tmp_reporter/";
   singleton::directories().setOutputDir(outdir);
   size_t iTmax = args.getValueOrFallback("--iTmax", 100); // maximum number of iterations
   size_t iTvtk = args.getValueOrFallback("--iTvtk", 1); // maximum number of iterations
   T amplitude = args.getValueOrFallback("--a", 2e-3); // maximum number of iterations
   
   const int ndim = 3; // a few things (e.g. SuperSum3D) cannot be adapted to 2D, but this should help speed it up
   const T physDeltaX          = lambda_phys / lambda_lat;  
   const int Nx                = nWaves*lambda_lat; 
   const T physLength         = 1.;       // length of the cuboid [m]
   const T domainlenth        =Nx*physDeltaX*2.; 
    
 
   // Welche Spitze auswerten? 1=erster Wellenberg, 2=zweiter, ...
   int peakN = args.getValueOrFallback("--peakN", 1);
 
   
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
 
   // --- Wellenzahl in physikalischen und lattice Einheiten ---
   // Nutze dieselbe λ bzw. wellenzahl wie in setBoundaryValues (hier: λ_phys = 0.5 m)
   
             
   const T k_phys = 2.*std::numbers::pi_v<T> / lambda_phys;         // [rad/m]
   const T k_lat  = k_phys * converter.getPhysDeltaX(); // [rad per lattice cell]
   const T k2_lat = k_lat * k_lat;
   const T cs_lat  = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());
   const T c_phys  = (converter.getPhysDeltaX()/converter.getPhysDeltaT()) * cs_lat;
   const T omegaPhys = c_phys * k_phys;   // EINDEUTIGE physikalische Frequenz
   
 
   // === 2nd Step: Prepare Geometry ===
   Vector<T,ndim> originFluid(-Domainelength/physLength/2., -physwidth/physLength/2., -physspan/physLength/2.);
   Vector<T,ndim> extendFluid(Domainelength/physLength, physwidth/physLength, physspan/physLength);
   IndicatorCuboid3D<T> domainFluid(extendFluid, originFluid);
   // -----------Variabeln definiere Messungen
   size_t nplot                  = args.getValueOrFallback( "--nplot",             100 );  
   size_t iTout                  = args.getValueOrFallback( "--iTout",             0   );  
     
   //----------------------------- Geometrie aufspannen
   Vector<T,ndim> extend{Domainelength/physLength, physwidth/physLength, physspan/physLength};
   Vector<T,ndim> origin{-Domainelength/physLength/2., -physwidth/physLength/2., -physspan/physLength/2.};
   IndicatorCuboid3D<T> cuboid(extend, origin);
   CuboidDecomposition3D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());
   cuboidDecomposition.setPeriodicity({false,false,false});
   HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
   
   SuperGeometry<T,ndim> superGeometry(cuboidDecomposition, loadBalancer);
   prepareGeometry(converter, superGeometry, domainFluid);
 
   // === 3rd Step: Prepare Lattice ===
   T rho0 = 1.0; 
   T u0 = converter.getCharLatticeVelocity();
   T alpha = 0.314;  // oder berechnet aus einem Pulsparameter
   T dampingDepthPU = 0.1;
   T lengthDomain = physLength;
   T dampingStrength = 1.0;
  
   SuperLattice<T,DESCRIPTOR> sLattice(superGeometry);
   prepareLattice( converter, sLattice, superGeometry, rho0, u0, amplitude, alpha,
                   dampingDepthPU, lengthDomain, dampingStrength);
  
   // === 4th Step: Main Loop with Timer ===
   std::size_t iTtimer = int(std::max(iTmax/20., 1.));
   // === calculate output intervals
   // nout is the minimum number of vtk outputs --> take max between nout and nout derived from iTout or tout
   
   size_t iTplot = std::min(std::max(int(iTmax / nplot), 1), 25);
 
   clout << "Timing setup:" << std::endl
         << "iTmax=" << iTmax << "; maxPhysT=" << converter.getPhysTime( iTmax ) << "; dt=" << converter.getPhysDeltaT() << std::endl
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
 //#ifdef FEATURE_WATCHPOINTS
 //Messpunkt X-Achse nehmen
   std::array<Vector<T,ndim>,2> measurePhysR = {
     Vector<T,ndim>{1.0, 0.0, 0.0},
     Vector<T,ndim>{1.0+lambda_phys/4., 0.0, 0.0}
   };
   std::array<Vector<int,4>,2> measureLatticeR{};
 
   // Messpunkt diagonal nehmen:
   std::array<Vector<T,ndim>,2> measurePhysDiagonal = {
     Vector<T,ndim>{physwidth/4, physwidth/4.,   physspan/4.},
     Vector<T,ndim>{physwidth/4+lambda_phys/4., physwidth/4.+lambda_phys/4.,   physspan/4.+lambda_phys/4.}
   };
   std::array<Vector<int,4>,2> measureLatticeDiagonal{};
 
   SuperD<T,descriptors::D3<fields::PHYS_R,descriptors::SCALAR>> watchpointsD(loadBalancer);
   SuperD<T,descriptors::D3<fields::PHYS_R,descriptors::SCALAR>> watchpointsDiagonal(loadBalancer);
 // Messpunkt 1
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
       std::cout << "[WARNUNG] X-Achse Messpunkt " << k << " wurde NICHT lokal gefunden!\n";
     }
   }
 
 // Messpunkt 2
 for (int k=0; k<2; ++k) {
   if (auto latticeDiagonal = cuboidDecomposition.getLatticeR(measurePhysDiagonal[k])) {
     measureLatticeDiagonal[k] = *latticeDiagonal;
 
     if (loadBalancer.isLocal(measureLatticeDiagonal[k][0])) {
       auto& blockDiagonal = watchpointsDiagonal.getBlock(measureLatticeDiagonal[k][0]);
       // Platz schaffen: 1 Zelle je Watchpoint am Blockende
       blockDiagonal.resize({blockDiagonal.getNx()+1,1,1});
       auto watchpointDiagonal = blockDiagonal.get(blockDiagonal.getNcells()-1);
       watchpointDiagonal.template setField<fields::PHYS_R>(measurePhysDiagonal[k]);
 
       // Lokalen Zellenindex merken (wir nutzen measureLatticeR[k][1] dafür)
       measureLatticeDiagonal[k][1] = blockDiagonal.getNcells()-1;
     }
   } else if (singleton::mpi().getRank()==0) {
     std::cout << "[WARNUNG] DiagonalMesspunkt " << k << " wurde NICHT lokal gefunden!\n";
   }
 }
 
   watchpointsD.setProcessingContext(ProcessingContext::Simulation);
   watchpointsDiagonal.setProcessingContext(ProcessingContext::Simulation);
  
 
   // Kopplung: schreibt den berechneten Druck in die Watchpoints
   SuperLatticePointCoupling pressureO(PressureO{},
                                       names::NavierStokes{}, sLattice,
                                       names::Points{}, watchpointsD);
                                       
   SuperLatticePointCoupling pressureO_diag(PressureO{},
                                        names::NavierStokes{}, sLattice,
                                        names::Points{}, watchpointsDiagonal);
   // Zeitreihen und Konstanten für spätere Auswertung
   std::vector<T> p1, p2; 
   p1.reserve(iTmax); 
   p2.reserve(iTmax);
    // Messpunkt 2
   std::vector<T> p3, p4; 
   p3.reserve(iTmax); 
   p4.reserve(iTmax);

   const T dtPhys = converter.getPhysDeltaT();
   
   const T dx = std::sqrt(
     (measurePhysR[1][0] - measurePhysR[0][0]) * (measurePhysR[1][0] - measurePhysR[0][0]) +
     (measurePhysR[1][1] - measurePhysR[0][1]) * (measurePhysR[1][1] - measurePhysR[0][1]) +
     (measurePhysR[1][2] - measurePhysR[0][2]) * (measurePhysR[1][2] - measurePhysR[0][2])
   );
   const T dxDiagonal = std::sqrt(
    (measurePhysDiagonal[1][0] - measurePhysDiagonal[0][0]) * (measurePhysDiagonal[1][0] - measurePhysDiagonal[0][0]) +
    (measurePhysDiagonal[1][1] - measurePhysDiagonal[0][1]) * (measurePhysDiagonal[1][1] - measurePhysDiagonal[0][1]) +
    (measurePhysDiagonal[1][2] - measurePhysDiagonal[0][2]) * (measurePhysDiagonal[1][2] - measurePhysDiagonal[0][2])
  );
 
   CSV<T> csvWriter("Welle", ';', {"iT", "t", "p1(xAchse)", "p2(XAchse)","p3(diagonal)","p4(diagonal)"}, ".csv");
   /*CSV<T> csvSummary("cp_vs_k", ';',
     {"k_lat", "k2_lat", "cp_lat_xcorr", "cp_lat_phase", "cs_lat"},
     ".csv");*/
 
   // --- Druck-Functor einmal definieren ---
   SuperLatticePhysPressure3D<T, DESCRIPTOR> pressureF(sLattice, converter);
   AnalyticalFfromSuperF3D<T> pressureInterp(pressureF, true, true);
 
   // --- Abtastlinie anlegen (x von links nach rechts, y=z=0) ---
   const int NxLine = 200; // Auflösung für den Plot
   std::vector<T> x_phys(NxLine), p_max(NxLine, T(0)), p_rss(NxLine, T(0));
   std::vector<T> p_snap(NxLine, T(0)); // Snapshot-Werte
   std::size_t n_accum = 0;
 
   const T x0 = -physLength*Domainelength/2.;       // wie in deiner Geometrie
   const T x1 =  physLength*Domainelength/2.;
   const T y0 =  T(0), z0 = T(0);
   for (int i=0; i<NxLine; ++i) {
     x_phys[i] = x0 + (x1 - x0) * ( (i + T(0.5)) / T(NxLine) );
   }
  // #endif
 
  const T cellsPerLambdaConverter = lambda_phys / converter.getPhysDeltaX();
  clout << "Gitterpunkte pro Wellenlänge (aus Converter): " 
      << cellsPerLambdaConverter << std::endl;
  clout << "physViscosity: "<<physViscosity<< std::endl;
 
 // --- Abtastlinien anlegen ---

 // Linie 1: entlang x-Achse (wie bisher)
 std::vector<T> x_phys_x(NxLine), p_max_x(NxLine, T(0)), p_rss_x(NxLine, T(0));
 std::vector<T> p_snap_x(NxLine, T(0));
 std::size_t n_accum_x = 0;

 // Linie 2: Diagonale (Beispiel: von (x0,y0,z0) nach (x1,y1,z1))
 std::vector<T> x_phys_diag(NxLine), y_phys_diag(NxLine), z_phys_diag(NxLine);
 std::vector<T> p_max_diag(NxLine, T(0)), p_rss_diag(NxLine, T(0));
 std::vector<T> p_snap_diag(NxLine, T(0));
 std::size_t n_accum_diag = 0;

 // Diagonale: z.B. in y-Richtung von -physwidth/2 nach +physwidth/2
 const T y_start = -physwidth/2.;
 const T y_end   =  physwidth/2.;
 const T z_diag  = 0.; // konstant, wenn du nur in x-y-Ebene diagonal willst

 for (int i=0; i<NxLine; ++i) {
     const T s = (i + T(0.5)) / T(NxLine); // Parameter [0,1]

     // Linie entlang x (wie bisher)
     x_phys_x[i] = x0 + (x1 - x0) * s;

     // Diagonale
     x_phys_diag[i] = x0 + (x1 - x0) * s;
     y_phys_diag[i] = y_start + (y_end - y_start) * s;
     z_phys_diag[i] = z_diag;
 }
   //--------------------------------------- FOR SCHLEIFE-------------------------------------------------------------------------------------------
   for (std::size_t iT=0; iT < iTmax; ++iT) {
     // === 5th Step: Definition of Initial and Boundary Conditions ===
     setBoundaryValues(converter, sLattice, iT, superGeometry, amplitude,rho0,lambda_phys,omegaPhys);
     
    // #ifdef FEATURE_WATCHPOINTS
     // ------------------------- Messwerte nehmen
     // ------------------------- Messwerte nehmen (2 Sensoren)
     pressureO.execute();
     pressureO_diag.execute();
     watchpointsD.setProcessingContext(ProcessingContext::Evaluation);
     watchpointsDiagonal.setProcessingContext(ProcessingContext::Evaluation);
 
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

     std::vector<T> localPD(2, T(0)), globalPD(2, T(0));
     for (int k=0; k<2; ++k) {
       if (loadBalancer.isLocal(measureLatticeDiagonal[k][0])) {
        auto& blkD  = watchpointsDiagonal.getBlock(loadBalancer.loc(measureLatticeDiagonal[k][0]));
        auto  cellD = blkD.get(measureLatticeDiagonal[k][1]);
        const T puD = cellD.template getField<descriptors::SCALAR>();
        localPD[k] += converter.getPhysPressure(puD);           // in phys. Druck [Pa]
       }
     }
     #ifdef PARALLEL_MODE_MPI
     singleton::mpi().reduceVect(localPD, globalPD, MPI_SUM);
     singleton::mpi().bCast(globalPD.data(), globalPD.size());
     #else
     globalPD = localPD;
     #endif


 
     // Zeitreihen füllen und CSV schreiben
     p1.push_back(globalP[0]);
     p2.push_back(globalP[1]);
     //csvWriter.writeDataFile(iT, {converter.getPhysTime(iT), globalP[0], globalP[1]});
     p3.push_back(globalPD[0]);
     p4.push_back(globalPD[1]);
     csvWriter.writeDataFile(iT, {converter.getPhysTime(iT), globalP[0], globalP[1], globalPD[0], globalPD[1]});
 
 
     if ( iT%iTvtk == 0 ) {getGraphicalResults(sLattice, converter, iT, superGeometry, amplitude);}
 
    // #ifdef FEATURE_WATCHPOINTS
     //===Zwischenschritt Amplitudenverlauf=====
     // --- p(x,t) auf der Linie auslesen ---
       //===Zwischenschritt Amplitudenverlauf=====

    // Linie 1: entlang x-Achse
    for (int i=0; i<NxLine; ++i) {
      T out;
      T pos_x[3] = { x_phys_x[i], y0, z0 };
      pressureInterp(&out, pos_x); // [Pa]
      const T a = std::abs(out);
      if (a > p_max_x[i]) p_max_x[i] = a;
      p_rss_x[i] += out*out;
    }
    // Optional Snapshot auf der x-Achse
    // if (iT == iT_snapshot) { p_snap_x[i] = out; }
    ++n_accum_x;

    // Linie 2: Diagonale
    for (int i=0; i<NxLine; ++i) {
      T out;
      T pos_d[3] = { x_phys_diag[i], y_phys_diag[i], z_phys_diag[i] };
      pressureInterp(&out, pos_d); // [Pa]
      const T a = std::abs(out);
      if (a > p_max_diag[i]) p_max_diag[i] = a;
      p_rss_diag[i] += out*out;
    }
    // Optional Snapshot auf der Diagonalen
    // if (iT == iT_snapshot) { p_snap_diag[i] = out; }
    ++n_accum_diag;

       //#endif
       
     // === 6th Step: Collide and Stream Execution ===
     sLattice.collideAndStream();
     // === 7th Step: Computation and Output of the Results ===
     if ( iT%iTtimer == 0 ) {timer.update(iT); timer.printStep();}
     }
 
     
 //#ifdef FEATURE_WATCHPOINTS
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
 auto peaks3 = findPeakTimes(p3, dtPhys, guard, minAmp);
 auto peaks4 = findPeakTimes(p4, dtPhys, guard, minAmp);
 // --- gewünschten Peak wählen (1=erster, 2=zweiter, ...)
 //     (CLI: --peakN 2 für zweiten Wellenberg)
 const int n = std::max(1, peakN);

// Prüfen, ob ALLE vier Sensorpaare genug Peaks haben
if ((int)peaks1.size() >= n && (int)peaks2.size() >= n &&
    (int)peaks3.size() >= n && (int)peaks4.size() >= n) {

  const T t1 = peaks1[n-1];
  const T t2 = peaks2[n-1];
  const T t3 = peaks3[n-1];
  const T t4 = peaks4[n-1];

  // ---------------- Messpunkt 1: X-Achse ----------------
  T cp_peak_phys      = std::numeric_limits<T>::quiet_NaN();
  T cp_peak_lat       = std::numeric_limits<T>::quiet_NaN();
  T cs_lat_here       = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());
  T ratio             = std::numeric_limits<T>::quiet_NaN();

  if (t2 != t1 && std::isfinite(t1) && std::isfinite(t2)) {
    cp_peak_phys = dx / std::abs(t2 - t1); // [m/s]
    cp_peak_lat  = cp_peak_phys * converter.getPhysDeltaT() / converter.getPhysDeltaX();
    ratio        = cp_peak_lat / cs_lat_here;

    if (singleton::mpi().getRank() == 0) {
      std::cout << "[cp|PEAK#" << n << " X-Achse] t1=" << t1 << " s, t2=" << t2 << " s"
                << " -> c_p=" << cp_peak_phys << " m/s"
                << " | c_p_lat=" << cp_peak_lat
                << " | c_p/c_s=" << ratio << "\n";
    }
  } else {
    if (singleton::mpi().getRank() == 0) {
      std::cout << "[cp|PEAK#" << n << " X-Achse] ungültige Peak-Zeiten.\n";
    }
  }

  // ---------------- Messpunkt 2: Diagonale ----------------
  T cp_peak_phys_Diagonal     = std::numeric_limits<T>::quiet_NaN();
  T cp_peak_lat_Diagonal      = std::numeric_limits<T>::quiet_NaN();
  T cs_lat_here_Diagonal      = cs_lat_here;  // gleiches cs
  T ratio_Diagonal            = std::numeric_limits<T>::quiet_NaN();

  if (t4 != t3 && std::isfinite(t3) && std::isfinite(t4)) {
    cp_peak_phys_Diagonal = dxDiagonal / std::abs(t4 - t3); // [m/s]
    cp_peak_lat_Diagonal  = cp_peak_phys_Diagonal * converter.getPhysDeltaT() / converter.getPhysDeltaX();
    ratio_Diagonal        = cp_peak_lat_Diagonal / cs_lat_here_Diagonal;

    if (singleton::mpi().getRank() == 0) {
      std::cout << "[cp|PEAK#" << n << " Diagonale] t3=" << t3 << " s, t4=" << t4 << " s"
                << " -> c_p=" << cp_peak_phys_Diagonal << " m/s"
                << " | c_p_lat=" << cp_peak_lat_Diagonal
                << " | c_p/c_s=" << ratio_Diagonal << "\n";
    }
  } else {
    if (singleton::mpi().getRank() == 0) {
      std::cout << "[cp|PEAK#" << n << " Diagonale] ungültige Peak-Zeiten.\n";
    }
  }

  // ---------------- CSV schreiben (beide Messmethoden zusammen) ----------------
  if (singleton::mpi().getRank() == 0) {
    CSV<T> csvPeak("cp_peak_n", ';',
      {"k_lat","k2_lat","peakN",
       "cp_phys","cp_lat","cs_lat","cp_over_cs",
       "cp_phys_Diagonal","cp_lat_Diagonal","cs_lat_Diagonal","cp_over_cs_Diagonal"},
      ".csv");

    csvPeak.writeDataFile(0, {
      k_lat, k2_lat, T(peakN),
      cp_peak_phys, cp_peak_lat, cs_lat_here,        ratio,
      cp_peak_phys_Diagonal, cp_peak_lat_Diagonal, cs_lat_here_Diagonal, ratio_Diagonal
    });
  }

} else {
  // Dieser else gehört zum großen if(...) über die Peak-Anzahlen
  if (singleton::mpi().getRank() == 0) {
    std::cout << "[cp|PEAK#" << n << "] Nicht genug Peaks gefunden: "
              << "sensor1=" << peaks1.size()
              << ", sensor2=" << peaks2.size()
              << ", sensor3=" << peaks3.size()
              << ", sensor4=" << peaks4.size()
              << "\n";
  }
}

 //#endif
 //#ifdef FEATURE_WATCHPOINTS
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
     
       const int maxLag = std::min(N/4, 200);  // z.B.: bis 1/4 der Datenlänge, max 200
 
 
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
            // csvSummary.writeDataFile(0, {k_lat, k2_lat, cp_lat_xcorr, cp_lat_phase, cs_lat});
 
     }
     
     //===Amplitudenverlauf eintragen=======
     // RMS X-Achse
std::vector<T> p_rms_x(NxLine);
for (int i=0; i<NxLine; ++i) {
    p_rms_x[i] = std::sqrt(p_rss_x[i] / std::max<std::size_t>(n_accum_x,1));
}

// RMS Diagonale
std::vector<T> p_rms_diag(NxLine);
for (int i=0; i<NxLine; ++i) {
    p_rms_diag[i] = std::sqrt(p_rss_diag[i] / std::max<std::size_t>(n_accum_diag,1));
}

// --- CSVs schreiben ---

    // 1) Peak-Hüllkurve X
    CSV<T> csvAmpPeakX("amplitude_peak_vs_xaxis", ';', {"x_phys_m", "p_peak_Pa"}, ".csv");
    for (int i=0; i<NxLine; ++i)
        csvAmpPeakX.writeDataFile(i, {x_phys_x[i], p_max_x[i]});

    // 2) RMS-Hüllkurve X
    CSV<T> csvAmpRmsX("amplitude_rms_vs_xaxis", ';', {"x_phys_m", "p_rms_Pa"}, ".csv");
    for (int i=0; i<NxLine; ++i)
        csvAmpRmsX.writeDataFile(i, {x_phys_x[i], p_rms_x[i]});

    // 3) Peak-Hüllkurve Diagonale
    CSV<T> csvAmpPeakDiag("amplitude_peak_vs_diag", ';', {"x_phys_m", "y_phys_m", "z_phys_m", "p_peak_Pa"}, ".csv");
    for (int i=0; i<NxLine; ++i)
        csvAmpPeakDiag.writeDataFile(i, {x_phys_diag[i], y_phys_diag[i], z_phys_diag[i], p_max_diag[i]});

    // 4) RMS-Hüllkurve Diagonale
    CSV<T> csvAmpRmsDiag("amplitude_rms_vs_diag", ';', {"x_phys_m", "y_phys_m", "z_phys_m", "p_rms_Pa"}, ".csv");
    for (int i=0; i<NxLine; ++i)
        csvAmpRmsDiag.writeDataFile(i, {x_phys_diag[i], y_phys_diag[i], z_phys_diag[i], p_rms_diag[i]});

 
 
    timer.stop();
    timer.printSummary();
  }
  
  
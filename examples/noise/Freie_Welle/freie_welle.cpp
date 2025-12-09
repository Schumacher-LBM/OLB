  /*Notiz: FREIE WELLE
Das Programm kompiliert Fehlerfrei und kann durchgeführt werden. 
Bei dem ausführen des Programms muss mit --iTmax XXX eine Zahl angegeben werden, wieviele Iterationsdurchläufe das Programm durchläuft
Terminalbefehl: make clean; make; ./Freie_Welle&>log.txt --iTmax 150 --peakN 2 (Mit welchem Wellenberg wird die Geschwindigkeit berechnet?) --outdir tmp_periodic_05  --> Zusätzlich kann die Amplitude angegeben werden a--
make; lam=1; RNAME="test_lam${lam}"; ./Freie_Welle&>${RNAME}.log --iTmax 150 --lambda $lam --peakN 2 --outdir $RNAME^

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
 #include <complex>
 
 using namespace olb;
 
 using T = FLOATING_POINT_TYPE;
 using DESCRIPTOR = descriptors::D3Q19<>;
 using BulkDynamics   = BGKdynamics<T, DESCRIPTOR>;
 using SpongeDynamics = SpongeLayerDynamics<T, DESCRIPTOR, momenta::BulkTuple, equilibria::SecondOrder>;
 const int ndim = 3; // a few things (e.g. SuperSum3D) cannot be adapted to 2D, but this should help speed it up
    // T lambda_phys              = T(0.6);  // gewünschte Wellenlänge in m
    const int lambda_LU        = 20.;     //Die Wellenlänge soll auf 10 LU abgebildet werden
    const int nWaves            = 3.5;               // In der Domäne sollen 6 Wellen abgebildet werden
    const T physLength          = 1.;       // length of the cuboid [m]
    const T physspan           = 0.46;
    const T physwidth         = 0.46;
    const T physLidVelocity   = 1.0;         // velocity imposed on lid [m/s] Fuer die Machzahl relevant (Vorher 1.0, jetzt ein zehntel der Schallgeschwindigkeit)
    const T physViscosity     = 1.5e-2;     // kinetic viscosity of fluid [m*m/s] Fuer die Relaxationszeit verantwortlich    
    const T physDensity       = 1.;         // fluid density of air (20°C)[kg/(m*m*m)]
    const T physMaxT          = 0.5;        // maximal simulation time [s]
    const T physDeltaT        = 0.00078125;// Messung 1: 0.00078125;//((0.68255-0.5)/3)/physViscosity*physDeltaX*physDeltaX;// 0,68255, weil Tau 0,68255 sein soll. Vorher: physDeltaX/343.46;  // temporal spacing [s] t=physDeltaX/c_s (Vorher 0.00078125, Jetzt: 5,8e-5)

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
     Vector<T,3> center(-physLength*2.4/2., 0., 0.); // Zentrum der Welle
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
    std::size_t iT, SuperGeometry<T,ndim>& superGeometry,BoundaryType boundarytype, T amplitude, T rho0, T lambda_phys)
  {
  
    if (boundarytype == local) {
      auto domain = superGeometry.getMaterialIndicator({3});
  
      // 2 Sinusperioden in 40 Zeitschritten
      T Kreisfrequenz = 2. * std::numbers::pi_v<T> /40.0;  // WENN DU DAS AENDERST, MUSST DU ES AUCH IN DER MAIN AENDERN!!!
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
     
      
      T wellenzahl=2. * std::numbers::pi_v<T>/lambda_phys ;//k=2pi/lamda

      // Wellenzahl in lattice-Einheiten
      const T k_LU  = wellenzahl* converter.getPhysDeltaX();

      // Schallgeschwindigkeit in lattice-Einheiten
      const T cs_LU = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());

      // Kreisfrequenz in lattice- und physikalischen Einheiten
      const T omega_LU   = cs_LU * k_LU;                     // ω0 = cs * k
      const T kreisfrequenz = omega_LU / converter.getPhysDeltaT(); // [1/s]

      T phase =0.;
      //T time = converter.getPhysTime(iT);  // physikalische Zeit aus Lattice-Zeit
      T time= converter.getLatticeTime(iT);
      T cs=sqrt(T(1)/descriptors::invCs2<T,DESCRIPTOR>());
 
      olb::SchallwelleRho<3, T, DESCRIPTOR> schallquelle(rho0, amplitude, wellenzahl, kreisfrequenz, phase, time, converter);
      olb::SchallwelleGesch<3,T, DESCRIPTOR> schallquelle_geschwindigkeit(amplitude,wellenzahl,kreisfrequenz,phase,time,rho0,cs,converter);
      // AnalyticalConst3D<T,T> rhoF( schallquelle);
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
      //T                          ndatapoints = converter.getResolution(); // number of data points on line
      AnalyticalFfromSuperF3D<T> pressure_interpolation(pressure, true, true);
      T                          pmin(converter.getPhysPressure(-amplitude / 50));
      T                          pmax(converter.getPhysPressure(+amplitude / 50));
      
      const T ndatapoints=300;
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
   int peakN = args.getValueOrFallback("--peakN", 2);
   T lambda_phys=args.getValueOrFallback("--lambda",0.6);
   T nPer = args.getValueOrFallback("--nPer",40.);
    const int ndim = 3; // a few things (e.g. SuperSum3D) cannot be adapted to 2D, but this should help speed it up
    // T lambda_phys              = T(0.6);  // gewünschte Wellenlänge in m
    
    const int Nx                = 80.; //nWaves*lambda_LU; 
    const T physLength         = 1.;       // length of the cuboid [m]
    const T physDeltaX          =physLength/Nx; 
    const T lambda_LU = lambda_phys/physDeltaX;
    const T domainlenth        =physLength;

 
   
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
   const T k_LU  = k_phys * converter.getPhysDeltaX(); // [rad per lattice cell]
   const T k2_LU = k_LU * k_LU;
 
   // theoretisches c_s (lattice und physisch)
   const T cs_LU  = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());          // ≈ 1/√3
   const T cs_phys = (converter.getPhysDeltaX()/converter.getPhysDeltaT()) * cs_LU;  // nur Info
 
 
 
 
 
 
 
   // === 2nd Step: Prepare Geometry ===
   BoundaryType boundarytype = periodic;
   Vector<T,ndim> originFluid(0., 0., 0.);
   Vector<T,ndim> extendFluid(domainlenth, physwidth, physspan);
   IndicatorCuboid3D<T> domainFluid(extendFluid, originFluid);
   // -----------Variabeln definiere Messungen
   size_t nplot                  = args.getValueOrFallback( "--nplot",             100 );  
   size_t iTout                  = args.getValueOrFallback( "--iTout",             0   );  
     
   //----------------------------- Geometrie aufspannen
   Vector<T,ndim> extend{domainlenth, physwidth, physspan};
   Vector<T,ndim> origin{0., 0., 0.};
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
   std::size_t iTvtk = 10;//int(std::max(iTmax/vtkanzahl, 1.));
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
     Vector<T,ndim>{domainlenth-0.4, physwidth/2., physspan/2.},
     Vector<T,ndim>{domainlenth-0.3, physwidth/2., physspan/2.}
   };
   std::array<Vector<int,4>,2> measureLatticeR{};

 
   SuperD<T,descriptors::D3<fields::PHYS_R,descriptors::SCALAR>> watchpointsD(loadBalancer);
    // --- Gitterpunkte pro Wellenlänge ausgeben ---
    const T cellsPerLambda = lambda_phys / physDeltaX;
    clout << "Gitterpunkte pro Wellenlänge: " << cellsPerLambda << std::endl;
    const T cellsPerLambdaConverter = lambda_phys / converter.getPhysDeltaX();
    clout << "Gitterpunkte pro Wellenlänge (aus Converter): " 
        << cellsPerLambdaConverter << std::endl;
    clout << "physViscosity: "<<physViscosity<< std::endl;
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
    T omegaPerStep =  cs_LU*k_LU;
   
   const T omegaPhys = omegaPerStep / dtPhys;
 
   CSV<T> csvWriter("Welle", ';', {"iT", "t", "p1", "p2"}, ".csv");
   CSV<T> csvSummary("cp_vs_k", ';',
     {"dumb","k_LU", "k2_LU", "cp_LU_xcorr", "cp_LU_phase", "cs_LU","cp_over_cs_corr","cs_ana", "cp_over_cs_ana"},
     ".csv");
 
   // Messung Plot Amplitudenverlauf
   // --- Druck-Functor einmal definieren ---
   SuperLatticePhysPressure3D<T, DESCRIPTOR> pressureF(sLattice, converter);
   AnalyticalFfromSuperF3D<T> pressureInterp(pressureF, true, true);
 
   // --- Abtastlinie anlegen (x von links nach rechts, y=z=0) ---
   const int NxLine = 300; // Auflösung für den Plot
   std::vector<T> x_phys(NxLine), p_max(NxLine, T(0)), p_rss(NxLine, T(0));
   std::vector<T> p_snap(NxLine, T(0)); // Snapshot-Werte
   std::size_t n_accum = 0;
 
   const T x0 = -physLength*2.4/2.;       // wie in deiner Geometrie
   const T x1 =  physLength*2.4/2.;
   const T y0 =  T(0), z0 = T(0);
   for (int i=0; i<NxLine; ++i) {
     x_phys[i] = x0 + (x1 - x0) * ( (i + T(0.5)) / T(NxLine) );
   }
   // --- analytische Dispersion nach Krüger mit den Initialwerten berechnet ---
    //=== Uebergabe der BoundaryConditions
  if (boundarytype == periodic) {setBoundaryValues(converter, sLattice, 0, superGeometry, boundarytype, amplitude, rho0,lambda_phys);}
    
   const T omega0= cs_LU*T(2)*std::numbers::pi_v<T> /lambda_LU;
   const T t_vi= T(2)*(converter.getLatticeRelaxationTime()-T(1)/T(2) );
   const T cp_LU_over_cs_LU_analytic= 1.-(T(1)/T(36))*(omega0/cs_LU)*(omega0/cs_LU) + (T(1)/8.) *(omega0*t_vi)*(omega0*t_vi); // Formel 3.17. Heinrichs
   const T cp_LU_analytic = cp_LU_over_cs_LU_analytic*cs_LU;
   clout<< "CP/CS_LU ANALYTISCH: "<<cp_LU_over_cs_LU_analytic << " ; omega0 analytisch: "<<omega0<<" ; tvi analytisch: " <<t_vi<< " ; Cp_LU_analytisch: "<< cp_LU_analytic<< std::endl;

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


     // Kreuzkorrelation völlig neu implementiert:
     CSV<T> csvXCorr("cp_xcorr", ';',
      {"k_LU", "k2_LU",
       "cp_phys", "cp_LU", "cp_over_cs",
       "cp_LU_analytic", "cp_over_cs_analytic"},
      ".csv");
  
     // Hier werden die ersten Perioden verworfen, damit die Welle Zeit zum einschwingen hat
    //  auto drop_front = [&](std::vector<T>& v, std::size_t n){
    //   if (v.size()>n) v.erase(v.begin(), v.begin()+n);
    // };
    
    // const T Tper = 2.*std::numbers::pi_v<T> / omegaPhys;  // omegaPhys korrekt setzen!
    // const std::size_t Ndrop = (std::size_t)std::ceil(1.0 * Tper / dtPhys); // z.B. 1 Perioden
    // drop_front(p1, Ndrop);
    // drop_front(p2, Ndrop);

    // //
    // auto subtract_mean = [](std::vector<T>& v){
    //   T m = 0;
    //   for (auto x : v) m += x;
    //   m /= (T)v.size();
    //   for (auto& x : v) x -= m;
    // };
    // subtract_mean(p1);
    // subtract_mean(p2);
    

    // // Subsample Korrelation und Berechnung
    // auto estimateLagFromXcorr = [&](const std::vector<T>& a,
    //   const std::vector<T>& b) -> T
    // {
    // const int N = (int)std::min(a.size(), b.size());
    // if (N < 10) return T(0);  // zu wenig Daten

    // // Wir nehmen nur den letzten Teil (z.B. 5 Perioden), um sicher eingeschwungen zu sein:
    // const int Nuse = N; // oder z.B. min(N, (int)std::round(5*Tper/dtPhys));
    // const int start = N - Nuse;

    // // max. Lag: z.B. +/- Nuse/4 oder +/- (dx / cp)*1.5
    // const int maxLag = std::min(Nuse/4, Nuse-1);

    // std::vector<T> R(2*maxLag+1, T(0)); // R[lag + maxLag]

    // // Kreuzkorrelation R(lag) = sum a[n]*b[n+lag]
    //   for (int lag = -maxLag; lag <= maxLag; ++lag) {
    //   T num = 0;
    //   for (int n = 0; n < Nuse; ++n) {
    //   int i = start + n;
    //   int j = i + lag;
    //   if (j < 0 || j >= N) continue;
    //   num += a[i] * b[j];
    //   }
    //   R[lag + maxLag] = num;
    //   }

    //   // Maximum finden
    //   int kMax = 0;
    //   T Rmax = std::numeric_limits<T>::lowest();
    //   for (int k = 0; k < (int)R.size(); ++k) {
    //   if (R[k] > Rmax) { Rmax = R[k]; kMax = k; }
    //   }

    //   // Index → Lag in Samples
    //   int lag0 = kMax - maxLag;

    //   // Subsample-Parabel um das Maximum (wie bei deinen Peaks)
    //   T lagFine = (T)lag0;
    //   if (kMax > 0 && kMax < (int)R.size()-1) {
    //   const T Rm = R[kMax-1];
    //   const T Rc = R[kMax];
    //   const T Rp = R[kMax+1];
    //   const T a2 = Rm - 2*Rc + Rp;
    //   const T b2 = Rm - Rp;
    //   if (std::abs(a2) > T(1e-30)) {
    //   T delta = b2 / (2*a2);           // in [-0.5, 0.5] erwartet
    //   delta = std::clamp(delta, T(-0.5), T(0.5));
    //   lagFine += delta;
    //   }
    //   }

    //   return lagFine; // in Samples (kann nicht-ganzzahlig sein)
    //   };

    //   T lagFine = estimateLagFromXcorr(p1, p2);   // in Samples
    //   T dt = lagFine * dtPhys;                    // Zeitverschiebung
    //   T cp_xcorr = dx / std::abs(dt);             // [m/s]

    //   // cp_xcorr ist in phys [m/s]
    //   T cp_xcorr_LU = cp_xcorr * converter.getPhysDeltaT() / converter.getPhysDeltaX();
    //   T cp_xcorr_over_cs = cp_xcorr_LU / cs_LU;

    //   csvXCorr.writeDataFile(0, {
    //       k_LU,
    //       k2_LU,
    //       cp_xcorr,
    //       cp_xcorr_LU,
    //       cp_xcorr_over_cs,
    //       cp_LU_analytic,
    //       cp_LU_over_cs_LU_analytic
    //   });

    // Phasenmethode implementieren
    
    // const std::size_t N = std::min(p1.size(), p2.size());
    // std::complex<T> A1(0, 0), A2(0, 0);

    // for (std::size_t n = 0; n < N; ++n) {
    //   T t = n * dtPhys;
    //   T c = std::cos(omegaPhys * t);
    //   T s = std::sin(omegaPhys * t);
    //   // e^{-i ω t} = cos(ω t) - i sin(ω t)
    //   std::complex<T> e_minus_iwt(c, -s);

    //   A1 += p1[n] * e_minus_iwt;
    //   A2 += p2[n] * e_minus_iwt;
    // }

    // // Normierung (optional, ändert nur Amplituden, nicht die Phase):
    // A1 *= (T(2) / (T)N);
    // A2 *= (T(2) / (T)N);

    // T phi1 = std::arg(A1);  // in rad, Bereich (-π, π]
    // T phi2 = std::arg(A2);
    
    // T dphi = phi2 - phi1;
    
    // // In Bereich (-π, π] „wrappen“
    // while (dphi >  M_PI) dphi -= T(2)*M_PI;
    // while (dphi < -M_PI) dphi += T(2)*M_PI;
    
    // T cp_phys = std::abs(omegaPhys * dx / std::max(std::abs(dphi), T(1e-12)));  // [m/s]
    // T cp_LU   = cp_phys * converter.getPhysDeltaT() / converter.getPhysDeltaX();
    // // T cs_LU   = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());
    // T cp_over_cs = cp_LU / cs_LU;

    // if (singleton::mpi().getRank()==0) {
    //   CSV<T> csvPhase("cp_phase_method", ';',
    //     {"k_LU", "k2_LU",
    //      "cp_phys", "cp_LU", "cp_over_cs",
    //      "cp_LU_analytic", "cp_over_cs_analytic"},
    //     ".csv");
    
    //   csvPhase.writeDataFile(0, {
    //     k_LU,
    //     k2_LU,
    //     cp_phys,
    //     cp_LU,
    //     cp_over_cs,
    //     cp_LU_analytic,
    //     cp_LU_over_cs_LU_analytic
    //   });
    // }
    
    // Messmethode an einem Punkt
    auto drop_front = [&](std::vector<T>& v, std::size_t n){
      if (v.size()>n) v.erase(v.begin(), v.begin()+n);
    };
    
    const T omega_theo = cs_LU * (T(2)*std::numbers::pi_v<T> / lambda_LU) / converter.getPhysDeltaT(); 
    // oder: omega_theo = 2*pi*cs_phys / lambda_phys;
    
    const T Tper_theo = T(2)*std::numbers::pi_v<T> / omega_theo;
    const std::size_t Ndrop = (std::size_t)std::ceil(2.0 * Tper_theo / dtPhys); // z.B. 2 Perioden
    
    drop_front(p1, Ndrop);
    auto findPeakTimes = [&](const std::vector<T>& p, T dt, int guardSamples, T minAmp){
      std::vector<T> peaks;
      const int N = (int)p.size();
      const int iStart = std::min(std::max(guardSamples, 1), N-3);
   
      for (int i = iStart/*+1*/; i < N-1; ++i) {
        if (p[i] > minAmp && p[i] > p[i-1] && p[i] > p[i+1]) {
          const T a = p[i-1] - 2*p[i] + p[i+1];
          const T b = p[i-1] - p[i+1];
          T delta = T(0);
          if (std::abs(a) > T(1e-30)) {
            delta = b / (2*a);
            delta = std::clamp(delta, T(-0.5), T(0.5));
          }
          peaks.push_back( (i + delta) * dt ); // physische Zeit [s]
        }
      }
      return peaks;
   };
   
   const int guard = 5;
    T estAmp = T(0);
    for (int i = guard; i < (int)std::min<std::size_t>(p1.size(), guard+50); ++i) {
      estAmp = std::max(estAmp, std::abs(p1[i]));
    }
    const T minAmp = estAmp * T(0.2);

    auto peaks1 = findPeakTimes(p1, dtPhys, guard, minAmp);
    if (peaks1.size() >= 3) {
      std::vector<T> periods;
      for (std::size_t i = 0; i+1 < peaks1.size(); ++i) {
          periods.push_back(peaks1[i+1] - peaks1[i]); // Δt zwischen Peaks
      }
  
      // Mittelwert der Perioden
      T Tper_num = T(0);
      for (auto T_i : periods) Tper_num += T_i;
      Tper_num /= (T)periods.size();
  
      const T omega_num = T(2)*std::numbers::pi_v<T> / Tper_num; // [1/s]
      // schon in main berechnet:
      const T k_phys = 2.*std::numbers::pi_v<T> / lambda_phys;
      const T k_LU   = k_phys * converter.getPhysDeltaX();
     // const T cs_LU  = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());

          // cp in physikalischen Einheiten:
    const T cp_phys = omega_num / k_phys; // [m/s]

    // nach LU:
    const T cp_LU   = cp_phys * converter.getPhysDeltaT() / converter.getPhysDeltaX();
    const T cp_over_cs = cp_LU / cs_LU;

    if (singleton::mpi().getRank()==0) {
      CSV<T> csvFreq("cp_from_frequency", ';',
          {"dumb","k_LU","k2_LU","cp_phys","cp_LU","cp_over_cs",
           "cp_LU_analytic","cp_over_cs_analytic"},
          ".csv");

      csvFreq.writeDataFile(0, {
          k_LU,
          k_LU*k_LU,
          cp_phys,
          cp_LU,
          cp_over_cs,
          cp_LU_analytic,
          cp_LU_over_cs_LU_analytic
      });

      std::cout << "[cp|FREQ] Tper_num="<<Tper_num<<" s, omega_num="<<omega_num
                <<" rad/s -> c_p="<<cp_phys<<" m/s, c_p/c_s="<<cp_over_cs<<"\n";
  }
      } else {
        if (singleton::mpi().getRank()==0) {
            std::cout << "[cp|FREQ] Zu wenige Peaks gefunden: "<<peaks1.size()<<"\n";
        }
}



    timer.stop();
    timer.printSummary();
  }
  
  
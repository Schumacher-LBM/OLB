  /*Notiz: FREIE WELLE
Das Programm kompiliert Fehlerfrei und kann durchgeführt werden. 
Bei dem ausführen des Programms muss mit --iTmax XXX eine Zahl angegeben werden, wieviele Iterationsdurchläufe das Programm durchläuft
Terminalbefehl: make clean; make; ./Freie_Welle&>log.txt --iTmax 150 --Nper (Anzahl Zeitschritte pro Periode) --peakN 2 (Mit welchem Wellenberg wird die Geschwindigkeit berechnet?) --outdir tmp_periodic_05  --> Zusätzlich kann die Amplitude angegeben werden a--
make; nw=1; RNAME="test_nw${nw}"; ./Freie_Welle&>${RNAME}.log --iTmax 150 --nWaves $nw --peakN 2 --outdir $RNAME^
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
 // Eigenschaften Akustik  
    //T lambda_phys             = T(0.4);  // gewünschte Wellenlänge in m
    //const int nWaves          = 3.5;               // In der Domäne sollen 6 Wellen abgebildet werden
   // const T physLength        = nWaves*lambda_phys;       // length of the cuboid [m]
    //const T physspan          = 0.46;
    //const T physwidth         = 0.46;
    const T physLength        = 4.;   // charL: feste Domänenlänge in m (NICHT mehr ändern!)
    const int Nx              = 80;     // Auflösung (fix)

    //Eigenschaften Fluid
    const T cs_phys           = 1.;
    const T nu_phys           = 3e-2;// kinem. Viskosität [m²/s]
    const T physMaxT          = 0.5; // max. physikalische Simulationszeit [s]
    const T physDensity       = 1.;
    
    // Numerik
  
    const T Ma                = 0.02;
  
 typedef enum { periodic, local } BoundaryType;
 // ------- Alles hier drueber konstant halten

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
  
    T dx_MeasureR = converter.getConversionFactorLength();
  
    switch (boundarytype) {
    // eternal and damping: 3 is the actual fluid; periodic: 1 is the fluid
    case periodic:{
      superGeometry.rename(2, 1);
      break;
    } }
  
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
      sLattice.setParameter<descriptors::OMEGA>(omega);
    
      // Make the lattice ready for simulation
      sLattice.initialize();
    
      clout << "Prepare Lattice ... OK" << std::endl;
    } //prepareLattice
  
  void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
    SuperLattice<T, DESCRIPTOR>& sLattice,
    std::size_t iT, SuperGeometry<T,ndim>& superGeometry,BoundaryType boundarytype, T amplitude, T rho0, T lambda_phys, int Nper, T physDeltaT)
  {
    if (boundarytype == periodic && iT==0) {
      auto domain = superGeometry.getMaterialIndicator({1});
 
      T wellenzahl=2. * std::numbers::pi_v<T>/lambda_phys ;//k=2pi/lamda
      //T kreisfrequenz_LU=2. * std::numbers::pi_v<T> / (T)Nper;  // 2π / Nper Schritte
      T cs=sqrt(T(1)/T(3.));
      T kreisfrequenz_LU= cs*wellenzahl;
      T kreisfrequenz_PU = kreisfrequenz_LU/physDeltaT;
      T phase =0.;
      //T time = converter.getPhysTime(iT);  // physikalische Zeit aus Lattice-Zeit
      T time= converter.getLatticeTime(iT);
      
 
      olb::SchallwelleRho<3, T, DESCRIPTOR> schallquelle(rho0, amplitude, wellenzahl, kreisfrequenz_PU, phase, time, converter);
      olb::SchallwelleGesch<3,T, DESCRIPTOR> schallquelle_geschwindigkeit(amplitude,wellenzahl,kreisfrequenz_PU,phase,time,rho0,cs,converter);
      AnalyticalConst3D<T,T> uInf(0., 0., 0.);

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

   T       nWaves            = args.getValueOrFallback<T>("--nWaves", 4);;    // kannst du von Run zu Run ändern
   T       lambda_phys       = physLength / nWaves; // daraus ergibt sich die Wellenlänge
   //outdir += "_reporter";
   size_t maxLatticeT = args.getValueOrFallback("--iTmax", 0); // maximum number of iterations
   T amplitude = args.getValueOrFallback<T>("--a", 1e-3); // maximum number of iterations
   if (outdir == "") outdir = "./tmp/";
   else outdir = "./tmp_" + outdir + "/";
   singleton::directories().setOutputDir(outdir);
   // Welche Spitze auswerten? 1=erster Wellenberg, 2=zweiter, ...
   int peakN = args.getValueOrFallback<T>("--peakN", 1);

   
   OstreamManager clout( std::cout,"main" ); // writing all output first in a userdefined Buffer of type OMBuf. On a flush it spits out at first the userdefined text in squared brackets and afterwards everything from the buffer


   // Provide the unit converter the characteristic entities
   const std::size_t res  = Nx;          // Auflösung auf charL
    const T charL         = physLength;  // charPhysLength = Domänenlänge in x
    const T charV         = cs_phys;     // charPhysVelocity = Schallgeschwindigkeit

    UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> converter(
        res,                  // Auflösung: Zellen auf charL
        Ma/std::sqrt(3.),     // charLatticeVelocity (u_char_lat)
        charL,                // charPhysLength
        charV,                // charPhysVelocity
        nu_phys,              // physViscosity
        physDensity           // physDensity
    );
    converter.print();

  // Für später praktisch:
  const T physDeltaX = converter.getPhysDeltaX();
  const T physDeltaT = converter.getPhysDeltaT();

  // Wellenlänge in Lattice-Zellen:
  const int lambda_lat = (int)std::round(lambda_phys / physDeltaX);

  // Vollständige Domänenlänge in x (nur Info, sollte == physLength sein):
  const T domainlength = physLength;

  // --- viskose Zeitskala t_v nach Krüger ---
  const T t_v = lambda_phys * lambda_phys / (4.*M_PI*M_PI * nu_phys);
  const T t_v_lat  = t_v / physDeltaT;        // [steps]
  clout << "lambda_phys = " << lambda_phys << " m\n";
  clout << "t_v (viskose Zeitskala) = " << t_v << " s\n";

   // --- Wellenzahl in physikalischen und lattice Einheiten ---
   // Nutze dieselbe λ bzw. wellenzahl wie in setBoundaryValues (hier: λ_phys = 0.5 m)
                   
   const T k_phys = 2.*std::numbers::pi_v<T> / lambda_phys;         // [rad/m]
   const T k_lat  = k_phys * converter.getPhysDeltaX(); // [rad per lattice cell]
   const T k2_lat = k_lat * k_lat;

   // theoretisches c_s (lattice und physisch)
   const T cs_lat  = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());          // ≈ 1/√3
   const T cs_phys = (converter.getPhysDeltaX()/converter.getPhysDeltaT()) * cs_lat;  // nur Info
  // -----------Variabeln definiere Messungen
      size_t nplot                  = args.getValueOrFallback( "--nplot",             100 );  
      size_t iTout                  = args.getValueOrFallback( "--iTout",             0   );  
      int Nper = args.getValueOrFallback("--Nper", 30);  // Anzahl Zeitschritte pro Periode
   // === 2nd Step: Prepare Geometry ===
   BoundaryType boundarytype = periodic;
   Vector<T,ndim> originFluid(0., 0., 0.);
   T physwidth = 3*converter.getPhysDeltaX();
   Vector<T,ndim> extendFluid(domainlength, physwidth, physwidth);
   IndicatorCuboid3D<T> domainFluid(extendFluid, originFluid);
   

   //----------------------------- Geometrie aufspannen
   Vector<T,ndim> extend{domainlength, physwidth, physwidth};
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
   // ------Zwischenschritt: Messwerte nehmen //-------
    
   std::array<Vector<T,ndim>,2> measurePhysR = {
     Vector<T,ndim>{domainlength-0.12, physwidth/2., physwidth/2.},
     Vector<T,ndim>{domainlength-0.02, physwidth/2., physwidth/2.}
   };
   std::array<Vector<int,4>,2> measureLatticeR{};

 
   SuperD<T,descriptors::D3<fields::PHYS_R,descriptors::SCALAR>> watchpointsD(loadBalancer);
    // --- Gitterpunkte pro Wellenlänge ausgeben ---
    const T cellsPerLambda = lambda_phys / physDeltaX;
    clout << "Gitterpunkte pro Wellenlänge: " << cellsPerLambda << std::endl;
    const T cellsPerLambdaConverter = lambda_phys / converter.getPhysDeltaX();
    clout << "Gitterpunkte pro Wellenlänge (aus Converter): " 
        << cellsPerLambdaConverter << std::endl;
    clout << "physViscosity: "<<nu_phys<< std::endl;
    clout << "lambda: "<<lambda_phys<<std::endl;
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
   // Abstand zwischen den Punkten berechnen
   const T dx_MeasureR = std::sqrt(
     (measurePhysR[1][0] - measurePhysR[0][0]) * (measurePhysR[1][0] - measurePhysR[0][0]) +
     (measurePhysR[1][1] - measurePhysR[0][1]) * (measurePhysR[1][1] - measurePhysR[0][1]) +
     (measurePhysR[1][2] - measurePhysR[0][2]) * (measurePhysR[1][2] - measurePhysR[0][2])
   );
 
 
   // Für die optionale Phasenmethode: physikalische Kreisfrequenz der Anregung bestimmen
   T omegaPerStep = T(0);
   omegaPerStep =  sqrt(T(1.)/T(3.))*2. * std::numbers::pi_v<T>/lambda_phys ; //2. * std::numbers::pi_v<T>*2. /Nper;
   const T omegaPhys = omegaPerStep / dtPhys;
 
   CSV<T> csvWriter("Welle", ';', {"Dumb","iT", "t", "p1", "p2"}, ".csv");
   CSV<T> csvSummary("cp_vs_k", ';',
     {"Dumb","k_lat", "k2_lat", "cp_lat_xcorr", "cp_lat_phase", "cs_lat"},
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
 
   const T x0 = domainlength/2.;       // wie in deiner Geometrie
   const T x1 =  -domainlength/2.;
   const T y0 =  T(0), z0 = T(0);
   for (int i=0; i<NxLine; ++i) {
     x_phys[i] = x0 + (x1 - x0) * ( (i + T(0.5)) / T(NxLine) );
   }
 
   //=== Uebergabe der BoundaryConditions
   if (boundarytype == periodic) {setBoundaryValues(converter, sLattice, 0, superGeometry, boundarytype, amplitude, rho0,lambda_phys, Nper,physDeltaT);}

    // --- analytische Dispersion nach Krüger mit den Initialwerten berechnet ---

   const T omega0= cs_lat*T(2)*std::numbers::pi_v<T> /lambda_lat;
   const T t_vi= T(2)*(converter.getLatticeRelaxationTime()-T(1)/T(2) );
   const T cp_LU_over_cs_LU_analytic= 1. - (T(1)/8.) *(omega0*t_vi)*(omega0*t_vi);
   
   clout<< "[CP/CS_LU ANALYTISCH: ]"<<cp_LU_over_cs_LU_analytic << "omega0 analytisch: "<<omega0<<"tvi analytisch: " <<t_vi<<std::endl;


   //--------------------------------------- FOR SCHLEIFE-------------------------------------------------------------------------------------------
   for (std::size_t iT=0; iT < iTmax; ++iT) {
     // === 5th Step: Definition of Initial and Boundary Conditions ===
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
      globalP = localP;
     
   
 
     // Zeitreihen füllen und CSV schreiben
     p1.push_back(globalP[0]);
     p2.push_back(globalP[1]);
     csvWriter.writeDataFile(iT, {converter.getPhysTime(iT), globalP[0], globalP[1]});
         
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
       ++n_accum;
 
     // === 6th Step: Collide and Stream Execution ===
     sLattice.collideAndStream();
     // === 7th Step: Computation and Output of the Results ===
     if ( iT%iTtimer == 0 ) {timer.update(iT); timer.printStep();}
     }
 
    
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
 const T minAmp = estAmp * T(0.2); // 20% der frühen Spitze als Schwellwert // Ist das zu viel?
 
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
     const T cp_peak_phys = dx_MeasureR / std::abs(t2 - t1); // [m/s]
     const T cp_peak_lat  = cp_peak_phys * converter.getPhysDeltaT() / converter.getPhysDeltaX();
     const T cs_lat_here  = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());
     const T ratio        = cp_peak_lat / cs_lat_here;

     //const T dx_MeasureR_phys    = converter.getPhysDeltaX();
     //const T dt_phys    = converter.getPhysDeltaT();
    
     const T nu_lat     = nu_phys * physDeltaT / (physDeltaX* physDeltaX);
     const T tvi_lat    = T(2)*(converter.getLatticeRelaxationTime()-T(1)/T(2));
     const T omega_lat  = omegaPerStep;           // deine Gitter-Kreisfrequenz pro Schritt
     const T omega_tvi  = omega_lat * tvi_lat;
     const T omega_tvi2 = omega_tvi * omega_tvi;

 
     if (singleton::mpi().getRank()==0) {
       std::cout << "[cp|PEAK#" << n << "] t1="<<t1<<" s, t2="<<t2<<" s"
                 << " -> c_p="<<cp_peak_phys<<" m/s"
                 << " | c_p_lat="<<cp_peak_lat
                 << " | c_p/c_s="<<ratio << "\n";
     }
 
     if (singleton::mpi().getRank()==0) {
      const T omega0_lat  = 2. * M_PI / (T)Nper;
      const T omega0_phys = omega0_lat / physDeltaT;      // mit physDeltaT von converter
  
      const T omega0_t_v  = omega0_phys * t_v;
      const T omega0_t_v2 = omega0_t_v * omega0_t_v;
  
      std::cout << "[DISP] omega0_phys = " << omega0_phys << " 1/s\n";
      std::cout << "[DISP] (omega0 * t_v)^2 = " << omega0_t_v2 << "\n";
  
      CSV<T> csvDisp("dispersion_omega_tau", ';',
          {"dumb","k_lat", "k2_lat", "omega0_phys", "t_v", "omega0_t_v2", "cp_over_cs_Simulation","cp_over_cs_analytic"},
          ".csv");
      csvDisp.writeDataFile(0, {k_lat, k2_lat, omega0_phys, t_v, omega0_t_v2, ratio, cp_LU_over_cs_LU_analytic});
  }
  


     // CSV schreiben (eigene Datei oder an deine bestehende anhängen)
     CSV<T> csvPeak("cp_peak_n", ';',
      {"dumb","k_lat","k2_lat","peakN",
       "cp_phys","cp_lat","cs_lat","cp_over_cs",
       "omega_lat","tvi_lat","omega_tvi_sq","cp_over_cs_analytisch"},
      ".csv");

    csvPeak.writeDataFile(0,{
      k_lat, k2_lat, peakN,
      cp_peak_phys, cp_peak_lat, cs_lat_here, ratio,
      omega_lat, tvi_lat, T(0), cp_LU_over_cs_LU_analytic
    });

   } else {
     if (singleton::mpi().getRank()==0) std::cout << "[cp|PEAK#" << n << "] ungültige Peak-Zeiten.\n";
   }
 } else {
   if (singleton::mpi().getRank()==0) {
     std::cout << "[cp|PEAK#" << n << "] Nicht genug Peaks gefunden: "
               << "sensor1="<<peaks1.size()<<", sensor2="<<peaks2.size()<<"\n";
   }
 }
 
 
    //  // ------------------------- Auswertung: Cross-Correlation & Phasenmethode
 
    //  // Optional: Einschwingtransienten verwerfen (z.B. erste 1-2 Perioden)
    //  auto drop_front = [&](std::vector<T>& v, std::size_t n){
    //    if (v.size()>n) v.erase(v.begin(), v.begin()+n);
    //  };
    //  {
    //    // 1 Periode ≈ (2π / omegaPhys) Sekunden
    //    const T Tper = 2.*std::numbers::pi_v<T> / std::max(omegaPhys, T(1e-12));
    //    const std::size_t Ndrop = (std::size_t)std::ceil(1.0 * Tper / dtPhys); // 1 Periode
    //    drop_front(p1, Ndrop);
    //    drop_front(p2, Ndrop);
    //  }
 
    //  // Cross-Correlation (einfach, normalisiert, Lag um 0 herum suchen)
    //  auto xcorrLag = [&](const std::vector<T>& a, const std::vector<T>& b)->int {
    //    const int N = (int)std::min(a.size(), b.size());
    //    if (N<=3) return 0;
    //    // maximaler Lag heuristisch begrenzen
    //    const int maxLag = std::min( (int)std::round(0.5 * dx_MeasureR / std::max(dtPhys, T(1e-12))), N-1 );
 
    //    // Mittelwerte entfernen
    //    T ma=0, mb=0; 
    //    for(int i=0;i<N;++i){ ma+=a[i]; mb+=b[i]; } 
    //    ma/=N; mb/=N;
 
    //    T bestC = -1e300; 
    //    int bestLag = 0;
    //    for (int lag=-maxLag; lag<=maxLag; ++lag) {
    //      T num=0, da=0, db=0;
    //      for (int i=0;i<N;++i) {
    //        const int j = i+lag;
    //        if (j<0 || j>=N) continue;
    //        const T aa = a[i]-ma;
    //        const T bb = b[j]-mb;
    //        num += aa*bb; da += aa*aa; db += bb*bb;
    //      }
    //      if (da>0 && db>0) {
    //        const T c = num / std::sqrt(da*db);
    //        if (c>bestC) { bestC=c; bestLag=lag; }
    //      }
    //    }
    //    return bestLag;
    //  };
 
    //  int lag = xcorrLag(p1,p2);
    //  T dt = lag * dtPhys;
    //  T cp_xcorr = dx_MeasureR / std::max(std::abs(dt), T(1e-12));
 
    //  if (singleton::mpi().getRank()==0) {
    //    std::cout << "[cp|XCORR] dx_MeasureR="<<dx_MeasureR<<" m, lag="<<lag<<" Samples, dt="<<dt
    //              <<" s -> c_p="<<cp_xcorr<<" m/s\n";
    //  }
 
    //  // -------- Optional: Phasenmethode (benötigt korrekte omegaPhys) ----------------------------------------------------------------------------------------------------------------------------
    //  auto complexProj = [&](const std::vector<T>& p)->std::pair<T,T>{
    //    T A=0, B=0; // Re=A, Im=B (mit -sin für Im)
    //    const std::size_t N = p.size();
    //    for (std::size_t n=0;n<N;++n){
    //      const T t = n*dtPhys;
    //      A += p[n]*std::cos(omegaPhys*t);
    //      B += p[n]*std::sin(omegaPhys*t);
    //    }
    //    const T phase = std::atan2(-B, A); // Phase ∈ (-π, π]
    //    return {A, phase};
    //  };
 
    //  if (p1.size()>=8 && p2.size()>=8) {
    //    auto [A1,phi1] = complexProj(p1);
    //    auto [A2,phi2] = complexProj(p2);
 
    //    T dphi = phi2 - phi1;
    //    while (dphi >  M_PI) dphi -= 2*M_PI;
    //    while (dphi < -M_PI) dphi += 2*M_PI;
 
    //    const T cp_phase = std::abs(omegaPhys * dx_MeasureR / std::max(std::abs(dphi), T(1e-12)));
 
    //    if (singleton::mpi().getRank()==0) {
    //      std::cout << "[cp|PHASE] dphi="<<dphi<<" rad, omega="<<omegaPhys
    //                <<" rad/s -> c_p="<<cp_phase<<" m/s\n";
    //    }
 
    //      // Umrechnung nach lattice-Einheiten:
    //          const T cp_lat_xcorr = cp_xcorr * converter.getPhysDeltaT() / converter.getPhysDeltaX();
    //          const T cp_lat_phase = cp_phase * converter.getPhysDeltaT() / converter.getPhysDeltaX();
    //          csvSummary.writeDataFile(0, {k_lat, k2_lat, cp_lat_xcorr, cp_lat_phase, cs_lat});
 
    //  }
     
    //  //===Amplitudenverlauf eintragen=======
    //  // RMS aus Summe der Quadrate
    //  std::vector<T> p_rms(NxLine);
    //  for (int i=0; i<NxLine; ++i) {
    //    p_rms[i] = std::sqrt(p_rss[i] / std::max<std::size_t>(n_accum,1));
    //  }
 
    //  // --- CSVs schreiben ---
    //  // 1) Peak-Hüllkurve
    //  CSV<T> csvAmpPeak("amplitude_peak_vs_x", ';', {"x_phys_m", "p_peak_Pa"}, ".csv");
    //  for (int i=0; i<NxLine; ++i) csvAmpPeak.writeDataFile(i, {x_phys[i], p_max[i]});
 
    //  // 2) RMS-Hüllkurve
    //  CSV<T> csvAmpRms("amplitude_rms_vs_x", ';', {"x_phys_m", "p_rms_Pa"}, ".csv");
    //  for (int i=0; i<NxLine; ++i) csvAmpRms.writeDataFile(i, {x_phys[i], p_rms[i]});
 
    //  // 3) (optional) Snapshot
    //  CSV<T> csvSnap("amplitude_snapshot_vs_x", ';', {"x_phys_m", "p_snapshot_Pa"}, ".csv");
    //  for (int i=0; i<NxLine; ++i) csvSnap.writeDataFile(i, {x_phys[i], p_snap[i]});
 
 
 
    timer.stop();
    timer.printSummary();
  }
  
  
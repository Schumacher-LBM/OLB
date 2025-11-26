/*Notiz: Dieses Programm simuliert eine freie Schallwelle über die gesamte Domänemake clean; make; ./Freie_Welle&>log.txt --iTmax 150 --Nper (Anzahl Zeitschritte pro Periode) --peakN 2 (Mit welchem Wellenberg wird die Geschwindigkeit berechnet?) --outdir tmp_periodic_05  --> Zusätzlich kann die Amplitude angegeben werden a--
Ausführen: make clean; make; ./Freie_Welle&>log.txt --iTmax 150 --Nper (Anzahl Zeitschritte pro Periode) --peakN 2 (Mit welchem Wellenberg wird die Geschwindigkeit berechnet?) --outdir tmp_periodic_05  --> Zusätzlich kann die Amplitude angegeben werden a--

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

 
 // Alle wichtigen Variablen:
    const int ndim = 3; // a few things (e.g. SuperSum3D) cannot be adapted to 2D, but this should help speed it up
    T resolution= 80; // Auflösung. Konstant halten
    T Ma=           0.02; //Machzahl

  // Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter, SuperGeometry<T, ndim>& superGeometry,
    IndicatorF3D<T>& domainFluid)
        {
        OstreamManager clout(std::cout, "prepareGeometry");
        clout << std::endl << "Prepare Geometry ..." << std::endl;

        // all nodes to temporary type
        superGeometry.rename(0, 2);
        T dx = converter.getConversionFactorLength();
        superGeometry.rename(2, 1);
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

        // Material=3 --> bulk dynamics
        auto bulkIndicator = superGeometry.getMaterialIndicator(
        {1,3}); // for local bcs all around, corners remain at 2, so they are included here
        sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
        sLattice.setParameter<descriptors::OMEGA>(omega);

        // Make the lattice ready for simulation
        sLattice.initialize();

        clout << "Prepare Lattice ... OK" << std::endl;
    } //prepareLattice

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
    // == 1st Step: Initialization
   initialize(&argc, &argv);
   CLIreader args(argc, argv);
   // Terminalbefehle zur Aenderung der Simulation
   std::string outdir = args.getValueOrFallback<std::string>("--outdir", "");
   outdir += "_reporter";
   size_t maxLatticeT = args.getValueOrFallback("--iTmax", 0); // maximum number of iterations
   T amplitude = args.getValueOrFallback("--a", 1e-3); // maximum number of iterations
   if (outdir == "") outdir = "./tmp/";
   else outdir = "./" + outdir + "/";
   singleton::directories().setOutputDir(outdir);
   int peakN = args.getValueOrFallback("--peakN", 1);
   OstreamManager clout( std::cout,"main" ); // writing all output first in a userdefined Buffer of type OMBuf. On a flush it spits out at first the userdefined text in squared brackets and afterwards everything from the buffer
   
   // Unit Conevrter einstellen:
   UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> converter(
    resolution,                  // Auflösung: Zellen auf charL
    Ma/std::sqrt(3.),     // charLatticeVelocity (u_char_lat)
    charL,                // charPhysLength
    charV,                // charPhysVelocity
    nu_phys,              // physViscosity
    physDensity           // physDensity
);
converter.print();
}
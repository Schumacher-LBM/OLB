
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006 - 2025 Mathias J. Krause, Jonas Fietz,
 *                            Jonas Latt, Jonas Kratzke, Shota Ito
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

/* cavity2d.cpp:
 * This example illustrates a minimal working example for a
 * fluid simulation with OpenLB; a flow in a cuboid, the lid-
 * driven cavity.
 */

 #include <olb.h>
 #include </home/clara/OpenLB/release-1.8.0/examples/noise/noiseauxiliary/acousticPulse.h>
 
 // Notwendige OpenLB Includes und Standardbibliotheken
#include <olb.h> // Generischer OpenLB-Header, der viele nützliche Dinge enthält
#include <iostream> // Für Konsolenausgaben (z.B. Fortschritt)
#include <cmath>    // Für std::sin, std::sqrt und M_PI (für Pi)

// === Globale Typdefinitionen und Parameter ===
using namespace olb; // Vereinfacht die Nutzung von OpenLB-Klassen
using T = FLOATING_POINT_TYPE; // Typ für Gleitkommazahlen (z.B. double)
using DESCRIPTOR = descriptors::D2Q9<>; // LBM-Modell: 2 Dimensionen, 9 Geschwindigkeitsvektoren
using BulkDynamics = ConstRhoBGKdynamics<T,DESCRIPTOR>;

// Physikalische Parameter (aus Ihrem Originalcode übernommen)
const T physDeltaX        = 0.0078125;   // Gitterzellenabstand [m]
const T physDeltaT        = 0.00078125;  // Zeitschritt [s]
const T physLength        = 1.0;         // Länge des Simulationsgebiets [m]
const T physLidVelocity   = 0.0;         // Keine Deckelgeschwindigkeit für freie Welle
const T physViscosity     = 0.001;       // Kinematische Viskosität [m*m/s]
const T physDensity       = 1.0;         // Fluiddichte [kg/(m*m*m)]
const T physMaxT          = 0.05;        // Maximale Simulationszeit [s] (reduziert für schnellere Ergebnisse)

// === Gitterdimensionen ===
// Berechnet aus physikalischen Parametern, um konsistent zu sein
const int N_X = static_cast<int>(physLength / physDeltaX); // 1.0 / 0.0078125 = 128
const int N_Y = N_X / 2; // Für eine sichtbare 2D-Welle, z.B. 64

// === Parameter für die Sinuswelle ===
const T BACKGROUND_DENSITY = 1.0;     // Hintergrunddichte des Fluids
const T WAVE_AMPLITUDE = 0.005;       // Amplitude der Dichtewelle (z.B. 0.5% Abweichung)
const T WAVE_LENGTH = 20.0;           // Wellenlänge in Gittereinheiten (z.B. 20 Gitterpunkte pro Welle)
const T WAVE_INITIAL_PHASE_X = 0.0;   // Startposition der Welle (Phasenverschiebung in x-Richtung)
const T WAVE_DIRECTION_X = 1.0;       // Ausbreitungsrichtung: 1.0 für rechts, -1.0 für links

// === Hilfsklasse für die sinusförmige Dichteverteilung ===
// Berechnet die Dichte an jedem Gitterpunkt als Sinusfunktion.
struct SineDensityField : public AnalyticalF<2, T, T> {
    T rho0;
    T amplitude;
    T k; // Wellenzahl (2*PI / Wellenlänge)
    T initialPhaseX;
    T directionX;

    SineDensityField(T rho0, T amplitude, T wavelength, T initialPhaseX, T directionX)
        : AnalyticalF<2, T, T>(1), // Gibt einen Skalarwert (Dichte) zurück
          rho0(rho0), amplitude(amplitude), initialPhaseX(initialPhaseX), directionX(directionX) {
        k = 2.0 * M_PI / wavelength; // Wellenzahl berechnen
    }

    bool operator()(T output[], const T input[]) override {
        // input[0] ist die x-Koordinate, input[1] ist die y-Koordinate (in Gittereinheiten)
        T x_coord = input[0];
        output[0] = rho0 + amplitude * std::sin(k * (x_coord - initialPhaseX) * directionX);
        return true;
    }
};

// === Hilfsklasse für das sinusförmige Geschwindigkeitsfeld ===
// Berechnet den Geschwindigkeitsvektor an jedem Gitterpunkt, konsistent mit der Dichtewelle.
struct SineVelocityField : public AnalyticalF<2, T, Vector<T, 2>> {
    T rho0;
    T amplitude;
    T k; // Wellenzahl
    T initialPhaseX;
    T directionX;
    T cs; // Schallgeschwindigkeit des Mediums

    SineVelocityField(T rho0, T amplitude, T wavelength, T initialPhaseX, T directionX, T cs)
        : AnalyticalF<2, T, Vector<T, 2>>(2), // Gibt einen Vektor (Geschwindigkeit) zurück (2 Komponenten für 2D)
          rho0(rho0), amplitude(amplitude), initialPhaseX(initialPhaseX), directionX(directionX), cs(cs) {
        k = 2.0 * M_PI / wavelength; // Wellenzahl berechnen
    }

    bool operator()(Vector<T, 2> output[], const T input[]) override {
        T x_coord = input[0];
        // Geschwindigkeitsamplitude für akustische Welle (vereinfacht)
        T velocity_magnitude = (amplitude / (rho0 * cs)) * std::sin(k * (x_coord - initialPhaseX) * directionX);

        // Die Geschwindigkeit zeigt in die Ausbreitungsrichtung (nur x-Komponente)
        output[0][0] = velocity_magnitude * directionX; // x-Komponente
        output[0][1] = 0.0;                             // y-Komponente ist 0
        return true;
    }
};

// === Funktionen zur Vorbereitung der Simulation ===

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,2>& superGeometry)
{
  // Setzt alle Materialien auf 1 (Fluid), was für eine freie Welle ideal ist.
  superGeometry.rename(0, 1);
  superGeometry.getStatistics().print();
}

void prepareLattice(const UnitConverter<T,DESCRIPTOR>& converter,
                    SuperLattice<T, DESCRIPTOR>& sLattice,
                    SuperGeometry<T,2>& superGeometry)
{
  // Definiert die Dynamik (BGK-Modell) für Material 1 (Fluid)
  sLattice.defineDynamics<ConstRhoBGKdynamics<T,DESCRIPTOR>>(superGeometry, 1);
  // Setzt den Relaxationsparameter (omega) basierend auf der physikalischen Viskosität
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

// Hier werden die Startbedingungen für die Sinuswelle definiert
void setBoundaryValues(const UnitConverter<T,DESCRIPTOR>& converter,
  SuperLattice<T, DESCRIPTOR>& sLattice,
  std::size_t iT, SuperGeometry<T,2>& superGeometry)
{
    if (iT == 0) {   // Die Welle wird nur zum Zeitpunkt t=0 initialisiert

        // Schallgeschwindigkeit im Gitter (cs = sqrt(1/3) für D2Q9)
        T cs = std::sqrt(DESCRIPTOR::cs2());

        // Erstellen der Dichte- und Geschwindigkeitsfelder für die Sinuswelle
        SineDensityField sineDensity(BACKGROUND_DENSITY, WAVE_AMPLITUDE, WAVE_LENGTH, WAVE_INITIAL_PHASE_X, WAVE_DIRECTION_X);
        SineVelocityField sineVelocity(BACKGROUND_DENSITY, WAVE_AMPLITUDE, WAVE_LENGTH, WAVE_INITIAL_PHASE_X, WAVE_DIRECTION_X, cs);

        // Den gesamten Fluidbereich (Material 1) mit der Sinuswelle initialisieren
        auto fluidDomain = superGeometry.getMaterialIndicator(1);
        sLattice.defineRhoU(fluidDomain, sineDensity, sineVelocity);
        sLattice.iniEquilibrium(fluidDomain, sineDensity, sineVelocity);

        // Finalisiert die Initialisierung des Gitters
        sLattice.initialize();
    }
}

// Funktion zum Speichern der Ergebnisse
void getResults(const UnitConverter<T,DESCRIPTOR>& converter,
                SuperLattice<T, DESCRIPTOR>& sLattice,
                std::size_t iT, util::Timer<T> timer)
{
  // Ausgabe alle 10 Zeitschritte
  const std::size_t iTvtk = converter.getLatticeTime(physMaxT/10.);

  SuperVTMwriter2D<T> vtmWriter("sineWave"); // Dateipräfix für die Ausgabe
  SuperLatticePhysVelocity2D velocity(sLattice, converter); // Physikalische Geschwindigkeit
  SuperLatticePhysPressure2D pressure(sLattice, converter); // Physikalische Druck
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  if (iT == 0) {
    vtmWriter.createMasterFile(); // Erstellt die Master-Datei für ParaView
  }

  // Schreibt die VTK-Dateien
  if (iT % iTvtk == 0 && iT >= 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  // Gibt Statistiken aus
  if (iT % (iTvtk * 2) == 0) { // Z.B. alle 20 Zeitschritte
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

// === Hauptprogramm ===
int main(int argc, char* argv[])
{
  // 1. OpenLB initialisieren
  initialize(&argc, &argv);
  OstreamManager clout( std::cout,"main" );

  // Konverter für physikalische und Gittereinheiten
  const UnitConverter<T,DESCRIPTOR> converter (
    physDeltaX, physDeltaT, physLength, physLidVelocity, physViscosity, physDensity
  );
  converter.print();

  // 2. Geometrie vorbereiten
  Vector<T,2> extend{static_cast<T>(N_X * converter.getPhysDeltaX()), static_cast<T>(N_Y * converter.getPhysDeltaX())};
  Vector<T,2> origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extend, origin);
  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());

  // Wichtig: Periodische Randbedingungen für eine "freie Welle"
  cuboidDecomposition.setPeriodicity({true,true});
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer);
  prepareGeometry(converter, superGeometry); // Ruft die oben definierte Funktion auf

  // 3. Lattice vorbereiten
  SuperLattice<T,DESCRIPTOR> sLattice(superGeometry);
  prepareLattice(converter, sLattice, superGeometry); // Ruft die oben definierte Funktion auf

  // 4. Hauptsimulationsschleife
  const std::size_t iTmax = converter.getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    // 5. Initial- und Randbedingungen setzen (nur bei iT=0 relevant)
    setBoundaryValues(converter, sLattice, iT, superGeometry);
    // 6. Kollision und Strömung
    sLattice.collideAndStream();
    // 7. Ergebnisse berechnen und ausgeben
    getResults(converter, sLattice, iT, timer);
  }

  timer.stop();
  timer.printSummary();
  return 0;
}

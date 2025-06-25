
const T BACKGROUND_DENSITY = 1.0;     // Die Dichte des Fluids außerhalb der Welle
const T WAVE_AMPLITUDE = 0.005;       // Amplitude der Dichtewelle (z.B. 0.5% Abweichung von der Hintergrunddichte)
const T WAVE_LENGTH = 20.0;           // Wellenlänge in Gittereinheiten (z.B. 20 Gitterpunkte pro Welle)
const T WAVE_INITIAL_PHASE_X = 0.0;   // Startposition der Welle (Phasenverschiebung in x-Richtung)
const T WAVE_DIRECTION_X = 1.0;       // Ausbreitungsrichtung: 1.0 für rechts, -1.0 für links

// === Hilfsklasse für die sinusförmige Dichteverteilung ===
// Diese Klasse berechnet die Dichte an jedem Gitterpunkt.
struct SineDensityField : public olb::AnalyticalF<D, T, T> {
    T rho0;
    T amplitude;
    T k; // Wellenzahl (2*PI / Wellenlänge)
    T initialPhaseX;
    T directionX;

    SineDensityField(T rho0, T amplitude, T wavelength, T initialPhaseX, T directionX)
        : olb::AnalyticalF<D, T, T>(1), // Gibt einen Skalarwert (Dichte) zurück
          rho0(rho0), amplitude(amplitude), initialPhaseX(initialPhaseX), directionX(directionX) {
        k = 2.0 * M_PI / wavelength; // Wellenzahl berechnen
    }

    bool operator()(T output[], const T input[]) override {
        // input[0] ist die x-Koordinate, input[1] ist die y-Koordinate
        T x_coord = input[0];
        output[0] = rho0 + amplitude * std::sin(k * (x_coord - initialPhaseX) * directionX);
        return true;
    }
};

// === Hilfsklasse für das sinusförmige Geschwindigkeitsfeld ===
// Diese Klasse berechnet den Geschwindigkeitsvektor an jedem Gitterpunkt.
struct SineVelocityField : public olb::AnalyticalF<D, T, olb::Vector<T, D>> {
    T rho0;
    T amplitude;
    T k; // Wellenzahl
    T initialPhaseX;
    T directionX;
    T cs; // Schallgeschwindigkeit

    SineVelocityField(T rho0, T amplitude, T wavelength, T initialPhaseX, T directionX, T cs)
        : olb::AnalyticalF<D, T, olb::Vector<T, D>>(D), // Gibt einen Vektor (Geschwindigkeit) zurück
          rho0(rho0), amplitude(amplitude), initialPhaseX(initialPhaseX), directionX(directionX), cs(cs) {
        k = 2.0 * M_PI / wavelength; // Wellenzahl berechnen
    }

    bool operator()(olb::Vector<T, D> output[], const T input[]) override {
        T x_coord = input[0];
        // Die Geschwindigkeitsamplitude ist proportional zur Dichteamplitude
        // (für kleine Amplituden und viskositätsfreie Wellen)
        T velocity_magnitude = (amplitude / (rho0 * cs)) * std::sin(k * (x_coord - initialPhaseX) * directionX);

        // Die Geschwindigkeit zeigt in die Ausbreitungsrichtung (nur x-Komponente)
        output[0][0] = velocity_magnitude * directionX; // x-Komponente
        output[0][1] = 0.0;                             // y-Komponente ist 0
        return true;
    }
};
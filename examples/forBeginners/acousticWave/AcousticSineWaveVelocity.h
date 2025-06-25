#ifndef ACOUSTIC_SINE_WAVE_VELOCITY_H
#define ACOUSTIC_SINE_WAVE_VELOCITY_H

#include "olb2D.h" // Stellt grundlegende OpenLB-Typen wie AnalyticalF, Vector, etc. bereit
#include <cmath>   // Für std::sin und M_PI (für Pi)

namespace olb {

/**
 * @brief Definiert ein sinusförmiges Geschwindigkeitsfeld, das konsistent mit einer
 * akustischen Dichtewelle ist.
 *
 * Die Geschwindigkeit wird berechnet als: u = (amplitude / (rho0 * cs)) * sin(k * (x - x0) . direction) * direction
 * wobei k die Wellenzahl und cs die Schallgeschwindigkeit ist.
 *
 * @tparam ndim Dimension des Raumes (z.B. 2 für 2D).
 * @tparam T Datentyp (z.B. double).
 */
template <unsigned ndim, typename T>
class AcousticSineWaveVelocity : public AnalyticalF<ndim, T, Vector<T, ndim>> {
protected:
  T               rho0;       ///< Hintergrunddichte
  T               amplitude;  ///< Amplitude der Dichtewelle
  T               wavelength; ///< Wellenlänge der Sinuswelle
  Vector<T, ndim> x0;         ///< Zentrum für die Phasenanpassung der Welle
  Vector<T, ndim> direction;  ///< Normalisierter Vektor der Wellenausbreitungsrichtung
  T               cs;         ///< Schallgeschwindigkeit des Mediums

public:
  /**
   * @brief Konstruktor für AcousticSineWaveVelocity.
   * @param rho0 Die Hintergrunddichte.
   * @param amplitude Die Amplitude der Dichtewelle (dies ist die Dichteamplitude, nicht die Geschwindigkeitsamplitude).
   * @param wavelength Die Wellenlänge der Sinuswelle.
   * @param x0 Das Zentrum, um das die Phase der Welle ausgerichtet wird (in Gittereinheiten).
   * @param direction Der Vektor, der die Ausbreitungsrichtung der Welle angibt.
   * @param cs Die Schallgeschwindigkeit des Mediums.
   */
  AcousticSineWaveVelocity(T rho0, T amplitude, T wavelength, Vector<T, ndim> x0, Vector<T, ndim> direction, T cs)
      : AnalyticalF<ndim, T, Vector<T, ndim>>(ndim) // Diese AnalyticalF gibt ndim Skalarwerte (einen Vektor) zurück
      , rho0(rho0)
      , amplitude(amplitude)
      , wavelength(wavelength)
      , x0(x0)
      , direction(direction.getNormalized()) // Sicherstellen, dass die Richtung normalisiert ist
      , cs(cs)
  {}

  /**
   * @brief Operator zum Berechnen des Geschwindigkeitsvektors an einem gegebenen Punkt.
   * @param output Array, in das der berechnete Geschwindigkeitsvektor geschrieben wird (output[0][d] für Dimension d).
   * @param input Array, das die Koordinaten des Punktes enthält (input[d] für Dimension d).
   * @return true, wenn die Berechnung erfolgreich war.
   */
  bool operator()(Vector<T, ndim> output[], const T input[]) override
  {
    T k = 2.0 * M_PI / wavelength; // Wellenzahl

    // Vektor vom Zentrum x0 zum aktuellen Punkt input
    Vector<T, ndim> x_minus_x0;
    for (unsigned d = 0; d < ndim; d++) {
      x_minus_x0[d] = input[d] - x0[d];
    }

    // Skalarprodukt von (x - x0) und der Ausbreitungsrichtung
    T dot_product = 0.0;
    for (unsigned d = 0; d < ndim; d++) {
      dot_product += x_minus_x0[d] * direction[d];
    }

    // Berechnung der Geschwindigkeitsamplitude für die akustische Welle
    // (Diese Formel ist eine Vereinfachung für kleine Amplituden und viskositätsfreie Wellen)
    T velocity_magnitude = (amplitude / (rho0 * cs)) * std::sin(k * dot_product);

    // Setzen der Komponenten des Geschwindigkeitsvektors
    for (unsigned d = 0; d < ndim; d++) {
      output[0][d] = velocity_magnitude * direction[d];
    }
    return true;
  }
};

} // namespace olb

#endif


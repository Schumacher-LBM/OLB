#ifndef ACOUSTIC_SINE_WAVE_H
#define ACOUSTIC_SINE_WAVE_H

#include "olb2D.h" // Stellt grundlegende OpenLB-Typen wie AnalyticalF, Vector, etc. bereit
#include <cmath>   // Für std::sin und M_PI (für Pi)

namespace olb {

/**
 * @brief Definiert eine sinusförmige Dichteverteilung für die Initialisierung.
 *
 * Die Dichte wird berechnet als: rho = rho0 + amplitude * sin(k * (x - x0) . direction)
 * wobei k der Wellenzahl (2*PI / Wellenlänge) ist.
 *
 * @tparam ndim Dimension des Raumes (z.B. 2 für 2D).
 * @tparam T Datentyp (z.B. double).
 */
template <unsigned ndim, typename T>
class AcousticSineWave : public AnalyticalF<ndim, T, T> {
protected:
  T               rho0;       ///< Hintergrunddichte
  T               amplitude;  ///< Amplitude der Dichtewelle
  T               wavelength; ///< Wellenlänge der Sinuswelle
  Vector<T, ndim> x0;         ///< Zentrum für die Phasenanpassung der Welle
  Vector<T, ndim> direction;  ///< Normalisierter Vektor der Wellenausbreitungsrichtung

public:
  /**
   * @brief Konstruktor für AcousticSineWave.
   * @param rho0 Die Hintergrunddichte.
   * @param amplitude Die Amplitude der Dichtewelle.
   * @param wavelength Die Wellenlänge der Sinuswelle.
   * @param x0 Das Zentrum, um das die Phase der Welle ausgerichtet wird (in Gittereinheiten).
   * @param direction Der Vektor, der die Ausbreitungsrichtung der Welle angibt.
   */
  AcousticSineWave(T rho0, T amplitude, T wavelength, Vector<T, ndim> x0, Vector<T, ndim> direction)
      : AnalyticalF<ndim, T, T>(1) // Diese AnalyticalF gibt 1 Skalarwert (Dichte) zurück
      , rho0(rho0)
      , amplitude(amplitude)
      , wavelength(wavelength)
      , x0(x0)
      , direction(direction.getNormalized()) // Sicherstellen, dass die Richtung normalisiert ist
  {}

  /**
   * @brief Operator zum Berechnen der Dichte an einem gegebenen Punkt.
   * @param output Array, in das die berechnete Dichte geschrieben wird (output[0]).
   * @param input Array, das die Koordinaten des Punktes enthält (input[d] für Dimension d).
   * @return true, wenn die Berechnung erfolgreich war.
   */
  bool operator()(T output[], const T input[]) override
  {
    T k = 2.0 * M_PI / wavelength; // Wellenzahl (k = 2*pi / lambda)

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

    // Berechnung der Dichte nach der Sinusfunktion
    output[0] = rho0 + amplitude * std::sin(k * dot_product);
    return true;
  }
};

} // namespace olb

#endif

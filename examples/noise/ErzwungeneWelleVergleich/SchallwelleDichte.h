/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *Diese Klasse berechn
*/

/* Wellengleichung
 * p(x)=rho0+A*sin(kx-ot+großes Phi)
 * u(x)=1/(rho_0*c_s)*p(x)
*/
// Gerade nur fuer die lokale Welle nutzbar!!!

#ifndef SchallwelleDichte
#define SchallwelleDichte
#include <cmath>


/* Berechnung Dichtefeld*/
namespace olb {

template <unsigned ndim, typename T, typename DESCRIPTOR>
class SchallwelleRho : public AnalyticalF<ndim, T, T> {
private:
    T roh0;
    T amplitude;
    T waveNumber; // k
    T omega;      // ω
    T phase;
    T time;
    int btype;
    UnitConverter<T,DESCRIPTOR> converter;

public:
    SchallwelleRho(T roh0, T amplitude_, T waveNumber_, T omega_, T phase_, T time_,int btype, UnitConverter<T,DESCRIPTOR> converter_)
        : AnalyticalF<ndim, T, T>(1),
          roh0(roh0),amplitude(amplitude_), waveNumber(waveNumber_),
          omega(omega_), phase(phase_), time(time_),btype(btype), converter(converter_) {}

    // bool operator()(T output[], const T input[]) override {
        /* Wellengleichung
        * p(x)=A*sin(kx+großes Phi)
        * u(x)=1/(rho_0*c_s)*p(x)
        */
        // const T x = input[0]; // nur x-Koordinate
        // output[0] = /*T(1)*/ roh0 + amplitude * std::sin(waveNumber * x - omega * time + phase);
        // return true;
    // }

    bool operator()(T output[], const T input[]) override {
        const T phi   = waveNumber*input[0] - omega*time + phase;

        // gleiche Hann-Rampe wie in Gesch-Functor
        const T Tramp = (T(2)*std::numbers::pi_v<T>)/omega;
        const T s     = std::clamp(time/Tramp, T(0), T(1));
        const T env   = T(0.5)*(T(1)-std::cos(std::numbers::pi_v<T>*s));

        // optionale räumliche Glättung (kleine Gaußblase)
        const T sigma = T(2) * converter.getPhysDeltaX();
        T r2 = T(0);
        // wenn du x0 nicht im Rho-Functor gespeichert hast, füge es analog wie im Gesch-Functor hinzu
        for (unsigned d=0; d<ndim; ++d) {
            const T dd = (input[d] /*- x0[d]*/); // ggf. x0 ergänzen
            r2 += dd*dd;
        }
        const T w = std::exp(-r2/(sigma*sigma));

        T p_phys = w * env * amplitude * std::sin(phi);      // [Pa]
        if (btype == 2) { /* local */ }                      // (Logik behältst du)
        // LU-Druck und rho
        const T p_lu  = converter.getLatticePressure(p_phys);
        const T rho_lu = p_lu / (T(1)/T(3));                 // c_s^2=1/3
        output[0] = roh0 + rho_lu;
        return true;
    }
    
};


}
    #endif

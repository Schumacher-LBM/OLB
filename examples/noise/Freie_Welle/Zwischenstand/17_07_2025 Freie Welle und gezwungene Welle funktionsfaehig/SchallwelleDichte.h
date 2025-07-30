/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *Diese Klasse berechn
*/

/* Wellengleichung
 * p(x)=rho0+A*sin(kx-ot+großes Phi)
 * u(x)=1/(rho_0*c_s)*p(x)
*/

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
    UnitConverter<T,DESCRIPTOR> converter;

public:
    SchallwelleRho(T roh0, T amplitude_, T waveNumber_, T omega_, T phase_, T time_, UnitConverter<T,DESCRIPTOR> converter_)
        : AnalyticalF<ndim, T, T>(1),
          roh0(roh0),amplitude(amplitude_), waveNumber(waveNumber_),
          omega(omega_), phase(phase_), time(time_), converter(converter_) {}

    bool operator()(T output[], const T input[]) override {
        /* Wellengleichung
        * p(x)=A*sin(kx+großes Phi)
        * u(x)=1/(rho_0*c_s)*p(x)
        */
        const T x = input[0]; // nur x-Koordinate
        output[0] = /*T(1)*/ roh0 + amplitude * std::sin(waveNumber * x - omega * time + phase);
        return true;
    }
};


}
    #endif

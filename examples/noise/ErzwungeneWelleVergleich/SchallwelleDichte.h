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
        const T x = input[0];                // [m], physisch
        const T t = time;                 // [s], physisch (per converter gegeben!)
        const T phi = waveNumber * x - omega * t + phase;
      
        // z.B. 1 Periodendauer rampen
        const T Tramp = (2*std::numbers::pi_v<T>)/omega;  // Rampdauer in s
        const T envelope = std::sin(std::min(T(1), t/Tramp) * std::numbers::pi_v<T>/T(2));

        T p_phys=0.;
        
        if (btype==2) //local
        {
            p_phys = envelope*amplitude * std::sin(phi);  // [Pa]
        }
        else{ p_phys = amplitude * std::sin(phi);          // [Pa] periodic
            }
        if (p_phys==0.)
        {std::cout << "p_phys ist 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";}
        
        const T p_lu   = converter.getLatticePressure(p_phys);
        // output[0]= p_lu;
        const T rho_lu = p_lu / (T(1)/T(3));                  // c_s^2 = 1/3
    
        output[0] = roh0 + rho_lu;                            // LU-Dichte zurückgeben
        return true;
    }
    
};


}
    #endif

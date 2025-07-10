/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *Diese Klasse berechn
*/

/* Wellengleichung
 * p(x)=A*sin(kx-ot+großes Phi)
 * u(x)=1/(rho_0*c_s)*p(x)
*/

#ifndef SchallwelleDruck
#define SchallwelleDruck
#include <cmath>


/* Berechnung Druckfeld*/
namespace olb {





  // template <unsigned ndim, typename T>
  // class SchallwelleDru : public AnalyticalF<ndim, T, T> {
  // protected:
  //   T amplitude;
  //   T waveNumber;   // k
  //   T frequency;    // ω
  //   T phase;
  //   Vector<T, ndim> center;
    
  
  // public:
  //   SchallwelleDru(T amplitude, T waveNumber, T frequency, T phase, Vector<T, ndim> center)
  //     : AnalyticalF<ndim, T, T>(1), amplitude(amplitude), waveNumber(waveNumber), frequency(frequency),
  //       phase(phase), center(center) {}
  
   
  
  //   bool operator()(T output[], const T input[]) override {
  //     T r2 = 0;
  //     for (unsigned d = 0; d < ndim; ++d) {
  //       T dx = input[d] - center[d];
  //       r2 += dx * dx;
  //     }
  //     T r = std::sqrt(r2) + 1e-8; // vermeidet Division durch 0
  //     output[0] = amplitude * sin(waveNumber * r - frequency + phase) / r;
  //     return true;
  //   }
  // };
  





  // // Meine Loesung. Funktioniert nicht wie gewuenscht  
    //     template <unsigned ndim, typename T>
    //     class SchallwelleDru : public AnalyticalF<ndim, T, T> 
    //     {
    // protected:
   
    //   T               amplitude;
    //   T               Wellenzahl;
    //   T               phase;
    // Vector<T, ndim> x0;
    // public:
    //   SchallwelleDru(T amplitude, T Wellenzahl, T phase, Vector<T, ndim> x0 = Vector<T, ndim>(0.))
    //       : AnalyticalF<ndim, T, T>(1)
    //       , amplitude(amplitude)
    //       , Wellenzahl(Wellenzahl)
    //       , phase(phase)
    //       ,x0(x0) {};
    
    //   bool operator()(T output[], const T input[]) override
    //   {
    //     //T distance = input[0]-x0[0];
    //     T distance = std::sqrt((input[0]-x0[0])* (input[0]-x0[0])+ (input[1]-x0[1])*(input[1]-x0[1])+ (input[2]-x0[2])*(input[2]-x0[2]));
    //     output[0] = amplitude * sin(Wellenzahl*distance+ phase);
    //     return true;
    //   };
    // };

    // /*Geschwindigkeitsfeld*/
    //     template <unsigned ndim, typename T>
    //     class AcousticPulse : public AnalyticalF<ndim, T, T> {
    //     protected:
       
    //       T               amplitude;
    //       T               Wellenzahl;
    //       T               phase;
    //       T               rho0;
    //       T               Geschwindigkeit;
    //       Vector<T, ndim> x0;
    //     //  T               Geschwindigkeit; // Von Descriptor ablesen?
        
        // public:
        //    SchallwelleGeschwindigkeit(T amplitude, T Wellenzahl, T phase,T rho0,T Geschwindigkeit,Vector<T, ndim> x0 = Vector<T,  ndim>(0.))
        //       : AnalyticalF<ndim, T, T>(1)
        //       , amplitude(amplitude)
        //       , Wellenzahl(Wellenzahl)
        //       , phase(phase)
        //       , rho0(rho0)
        //       , Geschwindigkeit(Geschwindigkeit)
        //       , x0(x0) {};
              
        
        //   bool operator()(T output[], const T input[]) override
        //   {
        //     T distance = 0;
        //     for (unsigned d = 0; d < ndim; d++) {
        //       T distance_i = input[d] - x0[d];
        //       distance += distance_i * distance_i; // actually, distance would be square root, but would be squared in exponent
        //     }
        //     output[0] = 1/(rho0*Geschwindigkeit)*(amplitude*sin(Wellennummer*distance+phase));
        //     return true;
        //   };
        // };

       // } // namespace olb


//         template <unsigned ndim, typename T>
//            class SchallwelleDru : public AnalyticalF<ndim, T, T> 
//           {
        
// protected:
//     T amplitude;
//     T Wellenzahl;  // = k
//     T omega;       // Kreisfrequenz ω
//     T phase;
//     T time;
//     Vector<T, ndim> x0;
//     T t = 0.0;     // aktuelle Zeit

// public:
//     SchallwelleDru(T amplitude, T Wellenzahl, T omega, T phase,T time, Vector<T, ndim> x0 = Vector<T, ndim>(0.))
//         : AnalyticalF<ndim, T, T>(1),
//           amplitude(amplitude),
//           Wellenzahl(Wellenzahl),
//           omega(omega),
//           phase(phase),
//           time(time),
//           x0(x0) {}


//     bool operator()(T output[], const T input[]) override {
//         T distance = 0.0;
//         if (time==0)
//         {
//           time=0.1;
//         }
        
//         for (unsigned i = 0; i < ndim; ++i) {
//             distance += (input[i] - x0[i]) * (input[i] - x0[i]);
//         }
//         distance = std::sqrt(distance);
//         output[0] = amplitude * std::sin(Wellenzahl * distance - omega * time + phase);
//         return true;
//     }
// };

template <unsigned ndim, typename T, typename DESCRIPTOR>
class SchallwelleDru : public AnalyticalF<ndim, T, T> {
private:
    T amplitude;
    T waveNumber; // k
    T omega;      // ω
    T phase;
    T time;
    UnitConverter<T,DESCRIPTOR> converter;

public:
    SchallwelleDru(T amplitude_, T waveNumber_, T omega_, T phase_, T time_, UnitConverter<T,DESCRIPTOR> converter_)
        : AnalyticalF<ndim, T, T>(1),
          amplitude(amplitude_), waveNumber(waveNumber_),
          omega(omega_), phase(phase_), time(time_), converter(converter_) {}

    bool operator()(T output[], const T input[]) override {
        /* Wellengleichung
        * p(x)=A*sin(kx+großes Phi)
        * u(x)=1/(rho_0*c_s)*p(x)
        */
        const T x = input[0]; // nur x-Koordinate
        output[0] = converter.getLatticeDensityFromPhysPressure(amplitude * std::sin(waveNumber * x - omega * time + phase));
        return true;
    }
};


}
    #endif

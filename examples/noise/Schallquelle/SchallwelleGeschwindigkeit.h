/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *Diese Klasse berechn
*/

/* Wellengleichung
 * 
 * u(x)=1/(rho_0*c_s)*p(x)
*/

#ifndef SchallwelleGeschwindigkeit
#define SchallwelleGeschwindigkeit

/* Berechnung Druckfeld*/
namespace olb {
  template <unsigned ndim, typename T>
  class SchallwelleGesch : public AnalyticalF<ndim, T, T> 
   {
        protected:
       
          T               amplitude;
          T               Wellenzahl;
          T               phase;
          T               rho0;
          T               Geschwindigkeit;
          Vector<T, ndim> x0;
        //  T               Geschwindigkeit; // Von Descriptor ablesen?
        
        public:
           SchallwelleGesch(T amplitude, T Wellenzahl, T phase,T rho0,T Geschwindigkeit,Vector<T, ndim> x0 = Vector<T,  ndim>(0.))
              : AnalyticalF<ndim, T, T>(1)
              , amplitude(amplitude)
              , Wellenzahl(Wellenzahl)
              , phase(phase)
              , rho0(rho0)
              , Geschwindigkeit(Geschwindigkeit)
              , x0(x0) {};
              
        
          bool operator()(T output[], const T input[]) override
          {
            T distance = input[0]-x0[0];
            output[0] = 1/(rho0*Geschwindigkeit)*(amplitude*sin(Wellenzahl*distance+phase));
            return true;
          };
        };
    



    } // namespace olb
    
    #endif
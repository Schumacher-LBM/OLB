/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *Diese Klasse berechn
*/

/* Wellengleichung
 * p(x)=A*sin(kx+großes Phi)
 * u(x)=1/(rho_0*c_s)*p(x)
*/

#ifndef SchallwelleGeschwindigkeit
#define SchallwelleGeschwindigkeit

/* Berechnung Druckfeld*/
namespace olb {
  template <unsigned ndim, typename T, typename DESCRIPTOR>
  class SchallwelleGesch : public AnalyticalF<ndim, T, T> 
   {
        protected:
       
          T               amplitude;
          T               Wellenzahl;
          T               kreisfrequenz;
          T               phase;
          T               time;
          T               rho0;
          T               cs;
          Vector<T, ndim> x0;
          UnitConverter<T,DESCRIPTOR> converter;
        //  T               Geschwindigkeit; // Von Descriptor ablesen?
        
        public:
           SchallwelleGesch(T amplitude, T Wellenzahl,T kreisfrequenz, T phase,T time, T rho0,T cs, UnitConverter<T,DESCRIPTOR> converter_,
            Vector<T, ndim> x0 = Vector<T,  ndim>(0.))
              : AnalyticalF<ndim, T, T>(1)
              , amplitude(amplitude)
              , Wellenzahl(Wellenzahl)
              , kreisfrequenz (kreisfrequenz)
              , phase(phase)
              , time(time)
              , rho0(rho0)
              , cs(cs)
              , x0(x0)
              , converter(converter_) {};
              
        
          bool operator()(T output[], const T input[]) override
          {
            /* Wellengleichung
             * p(x)=A*sin(kx+großes Phi)
             * u(x)=1/(rho_0*c_s)*p(x)
            */
            T distance = input[0];
            output[0] = 1./(rho0*cs)*util::pressureFromDensity<T,DESCRIPTOR>(rho0 + amplitude*sin(Wellenzahl*distance-kreisfrequenz*time+phase));
            for ( size_t i=1; i<ndim; i++ ) output[i] = 0.;
            return true;
          };
        };
    



    } // namespace olb
    
    #endif
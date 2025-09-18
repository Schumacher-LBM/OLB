/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *Diese Klasse berechn
*/

/* Wellengleichung
 * p(x)=A*sin(kx+großes Phi)
 * u(x)=1/(rho_0*c_s)*p(x)
*/
//MOMENTAN NUR FUER DEN LOKALEN FALL NUTZBAR!!!!!
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
            // T distance = input[0];
            // output[0] = 1./(rho0*cs)*util::pressureFromDensity<T,DESCRIPTOR>(rho0 + amplitude*sin(Wellenzahl*distance-kreisfrequenz*time+phase));
            // for ( size_t i=1; i<ndim; i++ ) output[i] = 0.;
            // return true;
            // Phase
            T phi = T(0);
            {
              T xdotk = T(0);
              // Falls du wirklich nur längs x anregen willst:
              // const T x = input[0];
              // xdotk = Wellenzahl * x;
              // Für allgemeine Richtung (hier nur x), sonst erweitern:
              xdotk = Wellenzahl * input[0];
              phi = xdotk - kreisfrequenz * time + phase;
            }

            // gleiche zeitliche Rampe wie im Dichte-Functor
            const T Tramp = (T(2)*std::numbers::pi_v<T>)/kreisfrequenz;
            const T s     = std::clamp(time/Tramp, T(0), T(1));
            const T env   = T(0.5)*(T(1)-std::cos(std::numbers::pi_v<T>*s)); // Hann

            // optionale räumliche Glättung (kleine Gaußblase)
            // sigma ~ 2-3 Zellen
            const T sigma = T(2) * converter.getPhysDeltaX();
            T r2 = T(0);
            for (unsigned d=0; d<ndim; ++d) {
              const T dd = (input[d] - x0[d]);
              r2 += dd*dd;
            }
            const T w = std::exp(-r2/(sigma*sigma));

            // physikalische Quelle
            const T p_phys = w * env * amplitude * std::sin(phi);  // [Pa]

            // Schallgeschwindigkeit in SI aus LU ableiten (robust):
            const T cs_lat  = std::sqrt(T(1) / descriptors::invCs2<T,DESCRIPTOR>());
            const T cs_phys = cs_lat * converter.getPhysDeltaX()/converter.getPhysDeltaT();

            const T u_phys = p_phys / (rho0 * cs_phys);            // [m/s]
            const T u_lu   = converter.getLatticeVelocity(u_phys); // [LU]

            // nur in x-Richtung anregen (ggf. Achse anpassen)
            output[0] = u_lu;
            for (unsigned d=1; d<ndim; ++d) output[d] = T(0);
            return true;
          }
          };
        };
    



    // namespace olb
    
    #endif
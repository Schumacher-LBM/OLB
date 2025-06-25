/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Adrian Kummerländer,
 *  Mathias J. Krause
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */
// Knudsen pump based on the temperature gradeint along the wall between two vessels
// The thermal creep flow leads to pressure gradient, which forces the fluid to stream backwards
// in the center of the pipe. The total MFR is therefore 0.
// Analytical solution is valid for b = 0 at Kn < 0.5

#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;

using DESCRIPTOR = D3Q19<NORMAL,VELOCITY,TEMPERATURE,TEMPGRADIENT>;
using BulkDynamics = BGKdynamics<T,DESCRIPTOR>::template wrap_collision<collision::SaveVelocity>;

// Parameters for the simulation setup
const int N = 41;
const T Kn = 0.45;
const T pM = 1000.;
const T accomodC = 1.;
const T bCoeff = -0.05;
const T Rspec = 287.;
const T temperature = 300.;
const T kB = 1.38e-23;
const T molarMass = 28.96e-3;
const T mass = molarMass/6.023e23;
const T density = mass*pM/kB/temperature;
const T thermSpeed = util::sqrt(3*8.314*temperature/molarMass);
const T kinD = 364e-12; //kinetic diamter of nitrogen [m]
const T sigma = std::numbers::pi*kinD*kinD;
const T visc = 2./3./util::sqrt(std::numbers::pi)*util::sqrt(molarMass*8.314*temperature)/sigma/6.023e23/density;
const T MFP = visc/thermSpeed;  //DOI: 10.1063/1.2185839
const T height = MFP/Kn;
const T length = 4.*height;
const T slipCoeff = (2. - accomodC)/accomodC/(1. + bCoeff*Kn)*MFP;
const T tempLeft = 299.;
const T tempRight = 301.;
const T creepCoeff = 3./4.*visc*density*Rspec/pM;
const T dP = -6.*visc*visc*density*density*Rspec/ (height/2.)/ (height/2.)/ pM * (tempRight - tempLeft)/(1. + 8.*(T(2) - accomodC)/accomodC*Kn/(1.+bCoeff*Kn));
const T uMax = height*height/T(2)/visc/density*dP/length*(0.25 - 0.5 - (T(2) - accomodC)/accomodC*Kn/(1.+bCoeff*Kn))+creepCoeff*(tempRight-tempLeft)/length;
const T maxPhysT = 2.e-7;
const T physInterval = maxPhysT/50.;
const T residuum = 1e-8;      // residuum for the convergence check

// variables for eoc analysis
T pressureL1RelError = 0;
T pressureL2RelError = 0;
T pressureLinfRelError = 0;

// Anayltical solution for velocity profile
template <typename T, typename S>
class PressureProfile3D : public AnalyticalF3D<T, S>
{
private:
    T radius, length, mu, pM, accomodC, Kn, xStart, xEnd;

public:
    PressureProfile3D(T radius_, T length_, T mu_, T pM_, T accomodC_, T Kn_, T xStart_, T xEnd_) : AnalyticalF3D<T, S>(1),
           radius(radius_), length(length_), mu(mu_), pM(pM_), accomodC(accomodC_), Kn(Kn_), xStart(xStart_), xEnd(xEnd_)
    {
        this->getName() = "PressureProfile3D";
    };

    bool operator()(T output[1], const S x[3])
    {
        T dP = 6.*mu*mu*Rspec/ radius/ radius/ pM * (tempRight - tempLeft)/(1. + 8.*(T(2) - accomodC)/accomodC*Kn/(1.+bCoeff*Kn));
        if(x[0] < xStart) {
          output[0] = pM-0.5*dP;
        }
        else if(x[0] > xEnd) {
          output[0] = pM+0.5*dP;
        }
        else{
          output[0] = pM-0.5*dP + dP/(xEnd - xStart)*(x[0] - xStart);
        }
        return true;
    };
};

// Temperature profile
template <typename T, typename S>
class TemperatureProfile3D : public AnalyticalF3D<T, S>
{
private:
    T tempLeft, tempRight, xStart, xEnd;

public:
    TemperatureProfile3D(T tempLeft_, T tempRight_, T xStart_, T xEnd_) : AnalyticalF3D<T, S>(1),
           tempLeft(tempLeft_), tempRight(tempRight_), xStart(xStart_), xEnd(xEnd_)
    {
        this->getName() = "TemperatureProfile3D";
    };

    bool operator()(T output[1], const S x[3])
    {
        if(x[0] < xStart) {
          output[0] = tempLeft;
        }
        else if(x[0] > xEnd) {
          output[0] = tempRight;
        }
        else{
          output[0] = tempLeft + (tempRight - tempLeft)/(xEnd - xStart)*(x[0] - xStart);
        }

        return true;
    };
};

// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,3>& superGeometry,
                      IndicatorF3D<T>& knudsenPumpLayer,
                      IndicatorF3D<T>& knudsenPump)
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2, knudsenPumpLayer);

  superGeometry.rename(2, 1, knudsenPump);

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
  superGeometry.communicate();

}

// Set up the geometry of the simulation
void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,3>& superGeometry,
                     IndicatorF3D<T>& channel)
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  sLattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator({1}));
  TemperatureProfile3D<T,T> temp(tempLeft, tempRight, 0., length);
  sLattice.defineField<TEMPERATURE>(superGeometry.getMaterialIndicator({0,2}), temp);
  setBouzidiBoundary<T,DESCRIPTOR,KnudsenVelocityPostProcessor<true>>(sLattice, superGeometry.getMaterialIndicator({2}), superGeometry.getMaterialIndicator({1}), channel);
  setBouzidiKnudsenSlipVelocity<T,DESCRIPTOR,true>(sLattice, superGeometry.getMaterialIndicator({2}), superGeometry.getMaterialIndicator({1}), channel);

  // Initial conditions
  T pL0 = converter.getLatticePressure(pM);
  AnalyticalConst3D<T,T> rho(pL0*descriptors::invCs2<T,DESCRIPTOR>()+T(1));
  AnalyticalConst3D<T,T> u0(T(0), T(0), T(0));
  sLattice.defineField<VELOCITY>(superGeometry.getMaterialIndicator({1,2}),u0);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry.getMaterialIndicator({0,1,2}), rho, u0);
  sLattice.iniEquilibrium(superGeometry.getMaterialIndicator({0,1,2}), rho, u0);

  sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  sLattice.template setParameter<KnudsenVelocityPostProcessor<true>::SLIPCOEFF>(slipCoeff/converter.getPhysDeltaX());
  sLattice.template setParameter<KnudsenVelocityPostProcessor<true>::CREEPCOEFF>(creepCoeff/converter.getPhysDeltaX()/converter.getConversionFactorVelocity());

  // Make the lattice ready for simulation
  sLattice.initialize();

  {
    auto& communicator = sLattice.getCommunicator(stage::PostStream());
    communicator.requestField<descriptors::POPULATION>();
    communicator.requestField<descriptors::VELOCITY>();
    communicator.requestOverlap(sLattice.getOverlap());
    communicator.exchangeRequests();
  }
  clout << "Prepare Lattice ... OK" << std::endl;
}

// Compute error norms
void error( SuperGeometry<T,3>& superGeometry,
            SuperLattice<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter )
{
  OstreamManager clout( std::cout,"error" );
  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  int tmp[]= { };
  T result[2]= { };

  // velocity error
  PressureProfile3D<T,T> presSol(height/2., length, density*visc, pM, accomodC, Kn, 0., length);
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pres( sLattice,converter );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  SuperRelativeErrorL1Norm3D<T> relPressureErrorNormL1(pres, presSol, indicatorF);
  relPressureErrorNormL1(result, tmp);
  clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;
  pressureL1RelError = result[0];

  SuperRelativeErrorL2Norm3D<T> relPressureErrorNormL2(pres, presSol, indicatorF);
  relPressureErrorNormL2(result, tmp);
  clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;
  pressureL2RelError = result[0];

  SuperRelativeErrorLinfNorm3D<T> relPressureErrorNormLinf(pres, presSol, indicatorF);
  relPressureErrorNormLinf(result, tmp);
  clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
  pressureLinfRelError = result[0];
}

// Output to console and files
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer, T maxPhysT, bool hasConverged)
{
  OstreamManager clout( std::cout,"getResults" );
  const bool lastTimeStep = ( hasConverged || (iT + 1 == converter.getLatticeTime( maxPhysT )) );
  const int statIter = converter.getLatticeTime( maxPhysT/50. );

  // VTK and image output only if no EOC analysis
  SuperVTMwriter3D<T> vtmWriter( "knudsenPump3d" );
  const int vtmIter  = converter.getLatticeTime( maxPhysT/50. );

  if ( iT==0 ) {
    // Writes the geometry, pipe no. and rank no. as vti file for visualization
    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || lastTimeStep ) {
    sLattice.executePostProcessors(stage::Evaluation{});
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter("knudsenPump3d");
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      SuperLatticePhysPressure3D pressure(sLattice, converter);
      SuperGeometryF<T,3> geometry( superGeometry );
      vtmWriter.addFunctor(velocity);
      vtmWriter.addFunctor(pressure);
      vtmWriter.addFunctor(geometry);
      PressureProfile3D<T,T> presSol(height/2., length, density*visc, pM, accomodC, Kn, 0., length);
      SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalPressureLattice(presSol, sLattice);
      analyticalPressureLattice.getName() = "analytical pressure";
      vtmWriter.addFunctor(analyticalPressureLattice);
      SuperLatticeField3D<T,DESCRIPTOR,descriptors::TEMPERATURE> temperatureF( sLattice );
      temperatureF.getName() = "temperatureF";
      vtmWriter.addFunctor( temperatureF );
      task(vtmWriter, iT);
    });
  }

  // Output on the console
  if ( iT%statIter==0 || lastTimeStep ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Error norms
    if ( lastTimeStep ) {
      error( superGeometry, sLattice, converter);
    }
  }
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");
  CLIreader args(argc, argv);
  const T relTime = args.getValueOrFallback("--relTime", 1.0);

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   relTime,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   height,     // charPhysLength: reference length of simulation geometry
    (T)   uMax,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   visc, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   density,    // physDensity: physical density in __kg / m^3__
    (T)   pM
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("knudsenPump3d");


  // === 2nd Step: Prepare Geometry ===
  Vector<T, 3> center0(0, height/2.-0.5 * converter.getPhysDeltaX(), height/2.-0.5 * converter.getPhysDeltaX());
  Vector<T, 3> center1(length, height/2.-0.5 * converter.getPhysDeltaX(), height/2.-0.5 * converter.getPhysDeltaX());
  std::shared_ptr<IndicatorCylinder3D<T> > pipe( new IndicatorCylinder3D<T>(center0, center1, height/2.));
  Vector<T, 3> centerL(-length/2., height/2.-0.5 * converter.getPhysDeltaX(), height/2.-0.5 * converter.getPhysDeltaX());
  std::shared_ptr<IndicatorCylinder3D<T> > tankL( new IndicatorCylinder3D<T>(centerL, center0, height));
  Vector<T, 3> centerR(1.5*length, height/2.-0.5 * converter.getPhysDeltaX(), height/2.-0.5 * converter.getPhysDeltaX());
  std::shared_ptr<IndicatorCylinder3D<T> > tankR( new IndicatorCylinder3D<T>(center1, centerR, height));
  auto knudsenPump = tankL + pipe + tankR;
  IndicatorLayer3D<T> knudsenPumpLayer( *knudsenPump, 2.*converter.getPhysDeltaX() );

  // Instantiation of a pipeGeometry with weights
  CuboidDecomposition3D<T> knudsenPumpGeometry(knudsenPumpLayer, converter.getPhysDeltaX(), singleton::mpi().getSize() );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( knudsenPumpGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( knudsenPumpGeometry, loadBalancer, 5 );

  prepareGeometry( converter, superGeometry, knudsenPumpLayer, *knudsenPump);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  // Prepare lattice and set boundary conditions
  prepareLattice(converter, sLattice, superGeometry, *knudsenPump);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ),
    superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();

  sLattice.communicate();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      getResults( sLattice, converter, iT, superGeometry, timer, maxPhysT, converge.hasConverged());

      break;
    }

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    sLattice.stripeOffDensityOffset(sLattice.getStatistics().getAverageRho()-(T)1);

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT,superGeometry, timer, maxPhysT, converge.hasConverged());
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), false );
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  timer.stop();
  timer.printSummary();

  singleton::pool().wait();
}

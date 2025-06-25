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

#include "olb3D.h"
#include "olb3D.hh"


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;

using DESCRIPTOR = D3Q19<NORMAL,VELOCITY,TEMPERATURE,TEMPGRADIENT>;
using BulkDynamics = BGKdynamics<T,DESCRIPTOR>::template wrap_collision<collision::SaveVelocity>;

// Parameters for the simulation setup
const T Kn = 0.45;
const T pM = 1000.;
const T dP = 0.5;
const T pOut = pM-dP/2.;
const T accomodC = 1.;
const T bCoeff = 1.;
const T Rspec = 287.;
const T temperature = 300.;
const T pIn = pOut + dP;
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
const T length = 2.*height;
const T slipCoeff = (2. - accomodC)/accomodC/(1. + bCoeff*Kn)*MFP;
const T tempLeft = 298.;
const T tempRight = 302.;
const T creepCoeff = 3./4.*visc*density*Rspec/pM;
const T uMax = height*height/T(2)/visc/density*(-pIn+pOut)/length*(0.25 - 0.5 - (T(2) - accomodC)/accomodC*Kn/(1.+bCoeff*Kn))+creepCoeff*(tempRight-tempLeft)/length;
const T maxPhysT = 2.e-7;       // max. simulation time in s, SI unit
const T physInterval = maxPhysT/50.;  // interval for the convergence check in s
const T residuum = 1e-6;      // residuum for the convergence check

// variables for eoc analysis
T velocityL1AbsError = 0;
T velocityL2AbsError = 0;
T velocityLinfAbsError = 0;
T velocityL1RelError = 0;
T velocityL2RelError = 0;
T velocityLinfRelError = 0;

// Anayltical solution for velocity profile
template <typename T, typename S>
class VelocityKnudsenProfile3D : public AnalyticalF3D<T, S>
{
private:
    T height, length, mu, pIn, pOut, accomodC, Kn;

public:
    VelocityKnudsenProfile3D(T height_, T length_, T mu_, T pIn_, T pOut_, T accomodC_, T Kn_) : AnalyticalF3D<T, S>(3),
           height(height_), length(length_), mu(mu_), pIn(pIn_), pOut(pOut_), accomodC(accomodC_), Kn(Kn_)
    {
        this->getName() = "VelocityKnudsenProfile3D";
    };

    bool operator()(T output[3], const S x[3])
    {
        T radius = util::sqrt((x[1]-height/2.)*(x[1]-height/2.) + (x[2]-height/2.)*(x[2]-height/2.));
        output[0] = height/2.*height/2./T(4)/mu*(pIn-pOut)/length*(1. - radius*radius/height/height*4. + 4.*(T(2) - accomodC)/accomodC*Kn/(T(1) + Kn))+creepCoeff*(tempRight-tempLeft)/length;
        output[1] = T(0);
        output[2] = T(0);

        return true;
    };
};

// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,3>& superGeometry)
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T, 3> center0(-converter.getPhysDeltaX() * 0.2, height/2., height/2.);
  Vector<T, 3> center1(length, height/2., height/2.);
  IndicatorCylinder3D<T> pipe(center0, center1, height/2.);

  superGeometry.rename(0, 2);

  superGeometry.rename(2, 1, pipe);

  superGeometry.clean();
  Vector<T, 3> origin(0, height/2., height/2.);
  Vector<T, 3> extend = origin;

  // Set material number for inflow
  origin[0] = -converter.getPhysDeltaX() * 2;
  extend[0] = converter.getPhysDeltaX() * 2;
  IndicatorCylinder3D<T> inflow(origin, extend, height/2.);
  superGeometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = length - 2 * converter.getPhysDeltaX();
  extend[0] = length + 2 * converter.getPhysDeltaX();
  IndicatorCylinder3D<T> outflow(extend, origin, height/2.);
  superGeometry.rename(2, 4, 1, outflow);

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

  sLattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator({1,3,4}));
  //sLattice.defineDynamics<EquilibriumBoundarySecondOrder>(superGeometry.getMaterialIndicator({3}));

  //boundary::set<boundary::LocalPressure<T,DESCRIPTOR,BulkDynamics>>(sLattice, superGeometry, 3);
  boundary::set<boundary::LocalVelocity<T,DESCRIPTOR,BulkDynamics>>(sLattice, superGeometry, 3);
  boundary::set<boundary::LocalPressure<T,DESCRIPTOR,BulkDynamics>>(sLattice, superGeometry, 4);
  AnalyticalLinear3D<T,T> temp((tempRight-tempLeft)/length, 0, 0, tempLeft);
  sLattice.defineField<TEMPERATURE>(superGeometry.getMaterialIndicator({0,2,5}), temp);
  setBouzidiBoundary<T,DESCRIPTOR,KnudsenVelocityPostProcessor<true>>(sLattice, superGeometry.getMaterialIndicator({2}), superGeometry.getMaterialIndicator({1,3,4}), channel);
  setBouzidiKnudsenSlipVelocity<T,DESCRIPTOR,true>(sLattice, superGeometry.getMaterialIndicator({2}), superGeometry.getMaterialIndicator({1,3,4}), channel);

  // Initial conditions
  T pL0 = converter.getLatticePressure(pOut);
  AnalyticalConst3D<T,T> rho(pL0*descriptors::invCs2<T,DESCRIPTOR>()+T(1));
  AnalyticalConst3D<T,T> u0(T(0), T(0), T(0));
  sLattice.defineField<VELOCITY>(superGeometry.getMaterialIndicator({1,2,3,4}),u0);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry.getMaterialIndicator({0,1,2,3,4}), rho, u0);
  sLattice.iniEquilibrium(superGeometry.getMaterialIndicator({0,1,2,3,4}), rho, u0);

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

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,3>& superGeometry)
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.05 );
  int iTupdate = iTmaxStart/30;

  if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T pLin = converter.getLatticePressure((pIn-pOut)*frac[0] + pOut);
    AnalyticalConst3D<T,T> rho(pLin*descriptors::invCs2<T,DESCRIPTOR>()+T(1));
    VelocityKnudsenProfile3D<T,T> uSol(height, length, density*visc, pIn, pOut, accomodC, Kn);
    AnalyticalConst3D<T,T> conv(frac[0]/converter.getConversionFactorVelocity(), 0., 0.);

    sLattice.defineRho( superGeometry, 3, rho );
    sLattice.defineU( superGeometry, 3, uSol*conv );

    sLattice.setProcessingContext<Array<momenta::FixedDensity::RHO>>(
      ProcessingContext::Simulation);
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
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
  VelocityKnudsenProfile3D<T,T> uSol(height, length, density*visc, pIn, pOut, accomodC, Kn);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( sLattice,converter );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  velocityL1AbsError = result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;
  velocityL1RelError = result[0];

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  velocityL2AbsError = result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;
  velocityL2RelError = result[0];

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  velocityLinfAbsError = result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;
  velocityLinfRelError = result[0];
}

// Output to console and files
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer, bool hasConverged,
                 Gnuplot<T>& gplot)
{
  OstreamManager clout( std::cout,"getResults" );
  const bool lastTimeStep = ( hasConverged || (iT + 1 == converter.getLatticeTime( maxPhysT )) );
  const int statIter = converter.getLatticeTime( maxPhysT/50. );

  // VTK and image output only if no EOC analysis
  SuperVTMwriter3D<T> vtmWriter( "microChannel3d" );
  const int vtmIter  = converter.getLatticeTime( maxPhysT/10. );

  if ( iT==0 ) {
    // Writes the geometry, pipe no. and rank no. as vti file for visualization
    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || lastTimeStep ) {
    sLattice.executePostProcessors(stage::Evaluation{});
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter("microChannel3d");
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      SuperLatticePhysPressure3D pressure(sLattice, converter);
      SuperGeometryF<T,3> geometry( superGeometry );
      vtmWriter.addFunctor(velocity);
      vtmWriter.addFunctor(pressure);
      vtmWriter.addFunctor(geometry);
      VelocityKnudsenProfile3D<T,T> uSol(height, length, density*visc, pIn, pOut, accomodC, Kn);
      SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalVelocityLattice(uSol, sLattice);
      analyticalVelocityLattice.getName() = "analytical solution";
      vtmWriter.addFunctor(analyticalVelocityLattice);
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
      error( superGeometry, sLattice, converter );
    }
  }

  // Gnuplot output
  if ((lastTimeStep)) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    gplot.setData (
      T(converter.getResolution()),
      { velocityL1RelError, velocityL2RelError, velocityLinfRelError},
      { "velocity L1 rel Error","velocity L2 rel Error",
        "velocity Linf rel error"},
      "top right",
      { 'p','p','p' } );
  }
}

void simulateMicroChannel3D( int N, Gnuplot<T>& gplot, T relTime)
{
  OstreamManager clout( std::cout,"simulatePoiseuilleNano" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   relTime,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   height,     // charPhysLength: reference length of simulation geometry
    (T)   uMax,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   visc, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   density,    // physDensity: physical density in __kg / m^3__
    (T)   pOut
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("microChannel3d");


  // === 2nd Step: Prepare Geometry ===
  Vector<T, 3> center0(0, height/2.-0.5 * converter.getPhysDeltaX(), height/2.-0.5 * converter.getPhysDeltaX());
  Vector<T, 3> center1(length + 0.5 * converter.getPhysDeltaX(), height/2.-0.5 * converter.getPhysDeltaX(), height/2.-0.5 * converter.getPhysDeltaX());
  IndicatorCylinder3D<T> pipe(center0, center1, height/2.);
  IndicatorLayer3D<T> pipeLayer( pipe, 3.*converter.getPhysDeltaX() );

  // Instantiation of a pipeGeometry with weights
  CuboidDecomposition3D<T> pipeGeometry(pipeLayer, converter.getPhysDeltaX(), singleton::mpi().getSize() );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( pipeGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( pipeGeometry, loadBalancer, 5 );

  prepareGeometry( converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  // Prepare lattice and set boundary conditions
  Vector<T, 3> center02(-10.*converter.getPhysDeltaX(), height/2., height/2.);
  Vector<T, 3> center12(length + 20 * converter.getPhysDeltaX(), height/2., height/2.);
  IndicatorCylinder3D<T> pipe2(center02, center12, height/2.);
  prepareLattice(converter, sLattice, superGeometry, pipe2);

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
      getResults( sLattice, converter, iT, superGeometry, timer, converge.hasConverged(), gplot);

      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry);

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT,superGeometry, timer, converge.hasConverged(), gplot);
    converge.takeValue( sLattice.getStatistics().getMaxU(), false );
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  timer.stop();
  timer.printSummary();

  singleton::pool().wait();
}

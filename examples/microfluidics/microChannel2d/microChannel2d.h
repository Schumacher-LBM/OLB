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
#include "olb2D.h"
#include "olb2D.hh"


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;

using DESCRIPTOR = D2Q9<NORMAL,VELOCITY,TEMPERATURE,TEMPGRADIENT>;
using BulkDynamics = BGKdynamics<T,DESCRIPTOR>::template wrap_collision<collision::SaveVelocity>;

// Parameters for the simulation setup
const T Kn = 0.45;
const T pM = 1000.;
const T dP = 0.5;
const T pOut = pM-dP/2.;
const T accomodC = 1.;
const T bCoeff = -1.;
const T Rspec = 296.8;
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
const T length = 10.*height;
const T slipCoeff = (2. - accomodC)/accomodC/(1. - bCoeff*Kn)*MFP;
const T tempLeft = 298.;
const T tempRight = 302.;
const T creepCoeff = 3./4.*visc*density*Rspec/pM;
const T uMax = height*height/T(2)/visc/density*(-pIn+pOut)/length*(0.25 - 0.5 - (T(2) - accomodC)/accomodC*Kn/(1.-bCoeff*Kn))+creepCoeff*(tempRight-tempLeft)/length;
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
class VelocityKnudsenProfile2D : public AnalyticalF2D<T, S>
{
private:
    T height, length, mu, pIn, pOut, accomodC, Kn;

public:
    VelocityKnudsenProfile2D(T height_, T length_, T mu_, T pIn_, T pOut_, T accomodC_, T Kn_) : AnalyticalF2D<T, S>(2),
           height(height_), length(length_), mu(mu_), pIn(pIn_), pOut(pOut_), accomodC(accomodC_), Kn(Kn_)
    {
        this->getName() = "VelocityKnudsenProfile2D";
    };

    bool operator()(T output[2], const S x[2])
    {
        output[0] = height*height/T(2)/mu*(-pIn+pOut)/length*(x[1]*x[1]/height/height - x[1]/height - (T(2) - accomodC)/accomodC*Kn/(T(1) - bCoeff*Kn))
                  + creepCoeff*(tempRight-tempLeft)/length;
        output[1] = T(0);

        return true;
    };
};

// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry,
                      IndicatorF2D<T>& cuboidLayer)
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2,cuboidLayer );
  Vector<T,2> extend( length-2.*converter.getPhysDeltaX(), height );
  Vector<T,2> origin(converter.getPhysDeltaX(),0.);
  IndicatorCuboid2D<T> channel( extend, origin );
  superGeometry.rename( 2,1,channel );
  superGeometry.clean();

  Vector<T,2> originBC = superGeometry.getStatistics().getMinPhysR( 2 );
  originBC[1] += converter.getPhysDeltaX()/2.;
  Vector<T,2> extendBC = superGeometry.getStatistics().getMaxPhysR( 2 );
  extendBC[1] = extendBC[1]-originBC[1]-converter.getPhysDeltaX()/2.;

  // Set material number for inflow
  originBC[0] = superGeometry.getStatistics().getMinPhysR( 2 )[0]-converter.getPhysDeltaX();
  extendBC[0] = 2*converter.getPhysDeltaX();
  IndicatorCuboid2D<T> inflow( extendBC,originBC );
  superGeometry.rename( 2,3,1,inflow );

  // Set material number for outflow
  originBC[0] = superGeometry.getStatistics().getMaxPhysR( 2 )[0]-converter.getPhysDeltaX();
  extendBC[0] = 2*converter.getPhysDeltaX();
  IndicatorCuboid2D<T> outflow( extendBC,originBC );
  superGeometry.rename( 2,4,1,outflow );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;

}

// Set up the geometry of the simulation
void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry,
                     IndicatorF2D<T>& channel)
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  sLattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator({1,3,4}));

  boundary::set<boundary::LocalPressure<T,DESCRIPTOR,BulkDynamics>>(sLattice, superGeometry, 3);
  boundary::set<boundary::LocalPressure<T,DESCRIPTOR,BulkDynamics>>(sLattice, superGeometry, 4);
  AnalyticalLinear2D<T,T> temp((tempRight-tempLeft)/length, 0, tempLeft);
  sLattice.defineField<TEMPERATURE>(superGeometry.getMaterialIndicator({0,2}), temp);
  setBouzidiBoundary<T,DESCRIPTOR,KnudsenVelocityPostProcessor<true>>(sLattice, superGeometry.getMaterialIndicator({2}), superGeometry.getMaterialIndicator({1,3,4}), channel);
  setBouzidiKnudsenSlipVelocity<T,DESCRIPTOR,true>(sLattice, superGeometry.getMaterialIndicator({2}), superGeometry.getMaterialIndicator({1,3,4}), channel);

  // Initial conditions
  T pL0 = converter.getLatticePressure(pOut);
  AnalyticalConst2D<T,T> rho(pL0*descriptors::invCs2<T,DESCRIPTOR>()+T(1));
  AnalyticalConst2D<T,T> u0(T(0), T(0));
  sLattice.defineField<VELOCITY>(superGeometry.getMaterialIndicator({0,1,2,3,4}), u0);


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
                        SuperGeometry<T,2>& superGeometry)
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.05 );

  if ( iT<= iTmaxStart ) {
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T pLin = converter.getLatticePressure((pIn-pOut)*frac[0] + pOut);
    AnalyticalConst2D<T,T> rho(pLin*descriptors::invCs2<T,DESCRIPTOR>()+T(1));

    sLattice.defineRho( superGeometry, 3, rho );

    sLattice.setProcessingContext<Array<momenta::FixedDensity::RHO>>(
      ProcessingContext::Simulation);
  }
}

// Compute error norms
void error( SuperGeometry<T,2>& superGeometry,
            SuperLattice<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter )
{
  OstreamManager clout( std::cout,"error" );
  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  int tmp[]= { };
  T result[2]= { };

  // velocity error
  VelocityKnudsenProfile2D<T,T> uSol(height, length, density*visc, pIn, pOut, accomodC, Kn);
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u( sLattice,converter );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm2D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  velocityL1AbsError = result[0];
  SuperRelativeErrorL1Norm2D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;
  velocityL1RelError = result[0];

  SuperAbsoluteErrorL2Norm2D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  velocityL2AbsError = result[0];
  SuperRelativeErrorL2Norm2D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;
  velocityL2RelError = result[0];

  SuperAbsoluteErrorLinfNorm2D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  velocityLinfAbsError = result[0];
  SuperRelativeErrorLinfNorm2D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;
  velocityLinfRelError = result[0];
}

// Output to console and files
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer, bool hasConverged,
                 Gnuplot<T>& gplot)
{
  OstreamManager clout( std::cout,"getResults" );
  const bool lastTimeStep = ( hasConverged || (iT + 1 == converter.getLatticeTime( maxPhysT )) );
  const int statIter = converter.getLatticeTime( maxPhysT/100. );

  // VTK and image output only if no EOC analysis
  SuperVTMwriter2D<T> vtmWriter( "microChannel2d" );
  const int vtmIter  = converter.getLatticeTime( maxPhysT/10. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || lastTimeStep ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperGeometryF<T,2> geom(superGeometry);
    vtmWriter.addFunctor( geom );

    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    vtmWriter.addFunctor( velocity );

    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( pressure );

    VelocityKnudsenProfile2D<T,T> uSol(height, length, density*visc, pIn, pOut, accomodC, Kn);
    SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> analyticalVelocityLattice(uSol, sLattice);
    sLattice.communicate();
    analyticalVelocityLattice.getName() = "analytical solution";
    vtmWriter.addFunctor(analyticalVelocityLattice);

    SuperLatticeField2D<T,DESCRIPTOR,descriptors::TEMPERATURE> temperatureF( sLattice );
    temperatureF.getName() = "temperatureF";
    vtmWriter.addFunctor( temperatureF );
    vtmWriter.write();
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

void simulateMicroChannel2D( int N, Gnuplot<T>& gplot)
{
  OstreamManager clout( std::cout,"simulateMicroChannel2D" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   1.13,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   height,     // charPhysLength: reference length of simulation geometry
    (T)   uMax,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   visc, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   density,    // physDensity: physical density in __kg / m^3__
    (T)   pOut
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("microChannel2d");


  // === 2nd Step: Prepare Geometry ===
  Vector<T,2> extend( length, height );
  Vector<T,2> origin;
  origin[1] -= 0.5*converter.getPhysDeltaX();
  IndicatorCuboid2D<T> cuboid( extend, origin );
  IndicatorLayer2D<T> cuboidLayer( cuboid, converter.getPhysDeltaX() );

  // Instantiation of a cuboidGeometry with weights
  CuboidDecomposition2D<T> cuboidGeometry(cuboidLayer, converter.getPhysDeltaX(), singleton::mpi().getSize() );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry( cuboidGeometry, loadBalancer, 5 );

  prepareGeometry( converter, superGeometry, cuboidLayer );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  // Prepare lattice and set boundary conditions
  origin[0] -= 10.*converter.getPhysDeltaX();
  origin[1] += 0.5*converter.getPhysDeltaX();
  extend[0] += 20.*converter.getPhysDeltaX();
  IndicatorCuboid2D<T> channel( extend, origin );
  prepareLattice(converter, sLattice, superGeometry, channel);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ),
    superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();

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

  timer.stop();
  timer.printSummary();
}

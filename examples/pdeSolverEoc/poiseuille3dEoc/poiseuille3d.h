/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Marc Haußmann, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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


#ifndef POISEUILLE_3D_H
#define POISEUILLE_3D_H

#include <olb.h>


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;

#ifdef ENABLE_MRT
using DESCRIPTOR = D3Q19<tag::MRT,FORCE>;
using BulkDynamics       = MRTdynamics<T,DESCRIPTOR>;
using ForcedBulkDynamics = ForcedMRTdynamics<T,DESCRIPTOR>;
#else
using DESCRIPTOR = D3Q19<FORCE>;
using BulkDynamics       = BGKdynamics<T,DESCRIPTOR>;
using ForcedBulkDynamics = ForcedBGKdynamics<T,DESCRIPTOR>;
#endif

typedef enum {forced, nonForced} FlowType;
typedef enum {bounceBack, local, interpolated, bouzidi, freeSlip, partialSlip} BoundaryType;

// Parameters for the simulation setup
FlowType flowType = nonForced;
BoundaryType boundaryType = interpolated;
bool eoc = false;

const T length  = 2.;         // length of the pipe
const T diameter  = 1.;       // diameter of the pipe
int N = 21;                   // resolution of the model
const T physU = 1.;           // physical velocity
const T Re = 10.;             // Reynolds number
const T physRho = 1.;         // physical density
const T tau = 0.8;            // lattice relaxation time
const T maxPhysT = 20.;       // max. simulation time in s, SI unit
const T residuum = 1e-5;      // residuum for the convergence check
const T tuner = 0.97;         // for partialSlip only: 0->bounceBack, 1->freeSlip

// Scaled Parameters
const T radius  = diameter/2.;            // radius of the pipe
const T physInterval = 0.0125*maxPhysT;   // interval for the convergence check in s

// variables for eoc analysis
T velocityL1AbsError = 0;
T velocityL2AbsError = 0;
T velocityLinfAbsError = 0;
T strainRateL1AbsError = 0;
T strainRateL2AbsError = 0;
T strainRateLinfAbsError = 0;
T pressureL1AbsError = 0;
T pressureL2AbsError = 0;
T pressureLinfAbsError = 0;

// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout(std::cout, "prepareGeometry");

  clout << "Prepare Geometry ..." << std::endl;

  Vector<T, 3> center0(-converter.getPhysDeltaX() * 0.2, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  if (flowType == forced) {
    center0[0] -= 3.*converter.getPhysDeltaX();
    center1[0] += 3.*converter.getPhysDeltaX();
  }
  IndicatorCylinder3D<T> pipe(center0, center1, radius);

  superGeometry.rename(0, 2);

  superGeometry.rename(2, 1, pipe);

  if (flowType == nonForced) {
    // to ensure the in/outlet set correctly
    superGeometry.clean();

    Vector<T, 3> origin(0, radius, radius);
    Vector<T, 3> extend = origin;

    // Set material number for inflow
    origin[0] = -converter.getPhysDeltaX() * 2;
    extend[0] = converter.getPhysDeltaX() * 2;
    IndicatorCylinder3D<T> inflow(origin, extend, radius);
    superGeometry.rename(2, 3, 1, inflow);

    // Set material number for outflow
    origin[0] = length - 2 * converter.getPhysDeltaX();
    extend[0] = length + 2 * converter.getPhysDeltaX();
    IndicatorCylinder3D<T> outflow(extend, origin, radius);
    superGeometry.rename(2, 4, 1, outflow);
  }

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>& sLattice,
                    UnitConverter<T, DESCRIPTOR>const& converter,
                    SuperGeometry<T,3>& superGeometry)
{


  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  if (flowType == nonForced) {
    sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);
  } else {
    sLattice.defineDynamics<ForcedBulkDynamics>(superGeometry, 1);
  }

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);

  std::vector<T> origin = { length, radius, radius};
  std::vector<T> axis = { 1, 0, 0 };

  CirclePoiseuille3D<T> poiseuilleU(origin, axis, converter.getCharLatticeVelocity(), radius);

  if (boundaryType == bounceBack) {
    boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  }
  else if (boundaryType == freeSlip) {
    boundary::set<boundary::FullSlip>(sLattice, superGeometry, 2);
  }
  else if (boundaryType == partialSlip) {
    boundary::set<boundary::PartialSlip>(sLattice, superGeometry, 2);
    sLattice.template setParameter<descriptors::TUNER>(tuner);
  }
  else if (boundaryType == bouzidi) {
    center0[0] -= 0.5*converter.getPhysDeltaX();
    center1[0] += 0.5*converter.getPhysDeltaX();
    if (flowType == forced) {
      center0[0] -= 3.*converter.getPhysDeltaX();
      center1[0] += 3.*converter.getPhysDeltaX();
    }
    IndicatorCylinder3D<T> pipe(center0, center1, radius);
    setBouzidiBoundary<T, DESCRIPTOR, BouzidiPostProcessor>(sLattice, superGeometry, 2, pipe);
  }
  else {
    if (boundaryType == local) {
      boundary::set<boundary::LocalVelocity>(sLattice, superGeometry, 2);
    }
    else {
      boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 2);
    }
  }

  if (flowType == nonForced) {
    if (boundaryType == bouzidi) {
      IndicatorCylinder3D<T> pipe(center0, center1, radius);
      setBouzidiBoundary<T,DESCRIPTOR,BouzidiVelocityPostProcessor>(sLattice, superGeometry, 3, pipe);
    }
    else {
      // Material=3 -->bulk dynamics
      if (boundaryType == local) {
        boundary::set<boundary::LocalVelocity>(sLattice, superGeometry, 3);
      }
      else {
        boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
      }
    }
    // Material=4 -->bulk dynamics
    if (boundaryType == local) {
      boundary::set<boundary::LocalPressure>(sLattice, superGeometry, 4);
    }
    else {
      boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);
    }
  }

  if (flowType == forced) {
    // Initial conditions
    T D = converter.getLatticeLength(diameter);

    std::vector<T> poiseuilleForce(3, T());
    poiseuilleForce[0] = 4. * converter.getLatticeViscosity() * converter.getCharLatticeVelocity() / (D * D / 4. );
    AnalyticalConst3D<T,T> force( poiseuilleForce );

    // Initialize force
    sLattice.defineField<FORCE>(superGeometry, 1, force);
    sLattice.defineField<FORCE>(superGeometry, 2, force );


    AnalyticalConst3D<T, T> rhoF(1);

    sLattice.defineRhoU(superGeometry, 1, rhoF, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 1, rhoF, poiseuilleU);
    sLattice.defineRhoU(superGeometry, 2, rhoF, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 2, rhoF, poiseuilleU);
  }
  else {
    // Initial conditions
    T p0 = 4. * converter.getPhysViscosity() * converter.getCharPhysVelocity() * length / (radius * radius);

    p0 = converter.getLatticePressure(p0);
    AnalyticalLinear3D<T, T> rho(-p0 / length * descriptors::invCs2<T,DESCRIPTOR>(), 0, 0, p0 * descriptors::invCs2<T,DESCRIPTOR>() + 1);

    std::vector<T> velocity(3, T());
    AnalyticalConst3D<T, T> uF(velocity);

    // Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU(superGeometry, 0, rho, uF);
    sLattice.iniEquilibrium(superGeometry, 0, rho, uF);
    sLattice.defineRhoU(superGeometry, 1, rho, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 1, rho, poiseuilleU);
    sLattice.defineRhoU(superGeometry, 2, rho, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 2, rho, poiseuilleU);
    sLattice.defineRhoU(superGeometry, 3, rho, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 3, rho, poiseuilleU);
    if(boundaryType == bouzidi) {
      setBouzidiVelocity(sLattice, superGeometry, 3, poiseuilleU);
    }
    sLattice.defineRhoU(superGeometry, 4, rho, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 4, rho, poiseuilleU);
  }

  // Set relaxation time to omega in all dynamics
  sLattice.setParameter<descriptors::OMEGA>(omega);
  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Compute error norms
void error( SuperGeometry<T,3>& superGeometry,
            SuperLattice<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter)
{
  OstreamManager clout( std::cout,"error" );
  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  int tmp[]= { };
  T result[2]= { };

  // velocity error
  const T maxVelocity = converter.getCharPhysVelocity();
  std::vector<T> axisPoint = {length, radius, radius};
  std::vector<T> axisDirection = { 1, 0, 0 };
  CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( sLattice,converter );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  velocityL1AbsError = result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  velocityL2AbsError = result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  velocityLinfAbsError = result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // strainRate error
  CirclePoiseuilleStrainRate3D<T, DESCRIPTOR> sSol( converter, radius );
  SuperLatticePhysStrainRate3D<T,DESCRIPTOR> s( sLattice,converter );

  SuperAbsoluteErrorL1Norm3D<T> absStrainRateErrorNormL1(s, sSol, indicatorF);
  absStrainRateErrorNormL1(result, tmp);
  clout << "strainRate-L1-error(abs)=" << result[0];
  strainRateL1AbsError = result[0];
  SuperRelativeErrorL1Norm3D<T> relStrainRateErrorNormL1(s, sSol, indicatorF);
  relStrainRateErrorNormL1(result, tmp);
  clout << "; strainRate-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absStrainRateErrorNormL2(s, sSol, indicatorF);
  absStrainRateErrorNormL2(result, tmp);
  clout << "strainRate-L2-error(abs)=" << result[0];
  strainRateL2AbsError = result[0];
  SuperRelativeErrorL2Norm3D<T> relStrainRateErrorNormL2(s, sSol, indicatorF);
  relStrainRateErrorNormL2(result, tmp);
  clout << "; strainRate-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absStrainRateErrorNormLinf(s, sSol, indicatorF);
  absStrainRateErrorNormLinf(result, tmp);
  clout << "strainRate-Linf-error(abs)=" << result[0];
  strainRateLinfAbsError = result[0];
  SuperRelativeErrorLinfNorm3D<T> relStrainRateErrorNormLinf(s, sSol, indicatorF);
  relStrainRateErrorNormLinf(result, tmp);
  clout << "; strainRate-Linf-error(rel)=" << result[0] << std::endl;

  if (flowType == nonForced) {
    // pressure error
    T p0 = 4. * converter.getPhysViscosity() * maxVelocity * length / (radius * radius);
    AnalyticalLinear3D<T, T> pressureSol(-p0 / length, 0, 0, p0);
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);

    SuperAbsoluteErrorL1Norm3D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
    absPressureErrorNormL1(result, tmp);
    clout << "pressure-L1-error(abs)=" << result[0];
    pressureL1AbsError = result[0];
    SuperRelativeErrorL1Norm3D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
    relPressureErrorNormL1(result, tmp);
    clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorL2Norm3D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
    absPressureErrorNormL2(result, tmp);
    clout << "pressure-L2-error(abs)=" << result[0];
    pressureL2AbsError = result[0];
    SuperRelativeErrorL2Norm3D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
    relPressureErrorNormL2(result, tmp);
    clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorLinfNorm3D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    absPressureErrorNormLinf(result, tmp);
    clout << "pressure-Linf-error(abs)=" << result[0];
    pressureLinfAbsError = result[0];
    SuperRelativeErrorLinfNorm3D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    relPressureErrorNormLinf(result, tmp);
    clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
  }
}

// Output to console and files
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer, bool hasConverged,
                 Gnuplot<T>& gplot, bool eoc)
{

  OstreamManager clout( std::cout,"getResults" );
  const bool lastTimeStep = ( hasConverged || (iT + 1 == converter.getLatticeTime( maxPhysT )) );
  const bool noslipBoundary = ((boundaryType != freeSlip) && (boundaryType != partialSlip));
  const int statIter = converter.getLatticeTime( maxPhysT/20. );

  // VTK and image output only if no EOC analysis
  if (! eoc) {

    SuperVTMwriter3D<T> vtmWriter( "poiseuille3d" );
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );

    const T maxVelocity = converter.getCharPhysVelocity();
    std::vector<T> axisPoint = {length, radius, radius};
    std::vector<T> axisDirection = { 1, 0, 0 };
    CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
    SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalVelocityLattice(uSol, sLattice);
    analyticalVelocityLattice.getName() = "analytical solution";
    vtmWriter.addFunctor(analyticalVelocityLattice);

    const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );

    if ( iT==0 ) {
      // Writes the geometry, cuboid no. and rank no. as vti file for visualization
      SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
      SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
      SuperLatticeDiscreteNormal3D<T, DESCRIPTOR> discreteNormal( sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3}) );
      SuperLatticeDiscreteNormalType3D<T, DESCRIPTOR> discreteNormalType( sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3, 4, 5}) );
      vtmWriter.write( cuboid );
      vtmWriter.write( rank );
      vtmWriter.write( discreteNormal );
      vtmWriter.write( discreteNormalType );

      vtmWriter.createMasterFile();
    }

    // Writes the vtm files and profile text file
    if ( iT%vtmIter==0 || lastTimeStep ) {
      sLattice.setProcessingContext(ProcessingContext::Evaluation);

      vtmWriter.write( iT );

      SuperEuklidNorm3D<T> normVel( velocity );
      BlockReduction3D2D<T> planeReduction( normVel, Vector<T,3>({0,0,1}), 600, BlockDataSyncMode::ReduceOnly );
      // write output as JPEG
      heatmap::write(planeReduction, iT);

    }
  }

  // Writes output on the console
  if ( iT%statIter==0 || lastTimeStep ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Error norms
    if (noslipBoundary) {
      if ( (!eoc) || lastTimeStep ) {
        error( superGeometry, sLattice, converter );
      }
    }
  }

  // Gnuplot output
  if ((noslipBoundary) && (lastTimeStep)) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    if ( eoc ) {
      if (flowType == nonForced){
        gplot.setData (
          T(converter.getResolution()),
          { velocityL1AbsError, velocityL2AbsError, velocityLinfAbsError,
            strainRateL1AbsError, strainRateL2AbsError, strainRateLinfAbsError,
            pressureL1AbsError, pressureL2AbsError, pressureLinfAbsError },
          { "velocity L1 abs Error","velocity L2 abs Error",
            "velocity Linf abs error",
            "strain rate L1 abs error", "strain rate L2 abs error",
            "strain rate Linf abs error",
            "pressure L1 abs error", "pressure L2 abs error",
            "pressure Linf abs error" },
          "top right",
          { 'p','p','p','p','p','p','p','p','p' } );
      } else {
        // same as above, but without pressure computation
        gplot.setData (
          T(converter.getResolution()),
          { velocityL1AbsError, velocityL2AbsError, velocityLinfAbsError,
            strainRateL1AbsError, strainRateL2AbsError, strainRateLinfAbsError },
          { "velocity L1 abs Error","velocity L2 abs Error",
            "velocity Linf abs error",
            "strain rate L1 abs error", "strain rate L2 abs error",
            "strain rate Linf abs error",},
          "top right",
          { 'p','p','p','p','p','p' } );
      }
    } else {  // if !eoc
      // plot velocity magnitude over line through the center of the simulation domain
      const T maxVelocity = converter.getCharPhysVelocity();
      T D = converter.getLatticeLength( diameter );
      T dx = 1. / T(converter.getResolution());
      T point[3] { };
      point[0] = length/2.;
      point[2] = ( T )radius;
      std::vector<T> axisPoint {length, radius, radius};
      std::vector<T> axisDirection { 1, 0, 0 };
      CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
      T analytical[3] { };
      SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
      AnalyticalFfromSuperF3D<T> intpolateVelocity( velocity, true, 1 );
      T numerical[3] { };
      for ( int iY=0; iY<=D; ++iY ) {
        point[1] = ( T )converter.getPhysLength(iY);
        uSol( analytical,point );
        intpolateVelocity( numerical,point );
        gplot.setData( iY*dx, {analytical[0],numerical[0]}, {"analytical","numerical"} );
      }
      // Create PNG file
      gplot.writePNG();
    }
  }
}

void simulatePoiseuille(int N, Gnuplot<T>& gplot, bool eoc)
{
  OstreamManager clout( std::cout,"simulatePoiseuille" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},                  // resolution: number of voxels per charPhysL
    (T)   tau,                // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   diameter,           // charPhysLength: reference length of simulation geometry
    (T)   physU,              // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   diameter*physU/Re,  // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physRho             // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("poiseuille3d");


  // === 2nd Step: Prepare Geometry ===

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length + 0.5 * converter.getPhysDeltaX(), radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  IndicatorLayer3D<T> extendedDomain(pipe, converter.getPhysDeltaX());

  // Instantiation of a cuboidDecomposition with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else // ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 1;
#endif // ifdef PARALLEL_MODE_MPI
  CuboidDecomposition3D<T> cuboidDecomposition( extendedDomain, converter.getPhysDeltaX(), noOfCuboids);
  if (flowType == forced) {
    // Periodic boundaries in x-direction
    cuboidDecomposition.setPeriodicity({ true, false, false });
  }

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // Instantiation of a superGeometry
  const int overlap = (flowType == forced) ? 2 : 3;
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, overlap);

  prepareGeometry(converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and setBoundaryConditions
  prepareLattice(sLattice, converter, superGeometry);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      getResults( sLattice, converter, iT, superGeometry, timer,
        converge.hasConverged(), gplot, eoc );

      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer,
      converge.hasConverged(), gplot, eoc  );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), (! eoc) );
  }

  timer.stop();
  timer.printSummary();
}

#endif

/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2021 Stephan Simonis, Davide Dapelo,
 *                2024 Marc Heinzelmann
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

/*  Transient Longitudinal Mixing:
 *  The displacement of an initially homogeneous solute
 *  by the introduction of a solvent.
 *
 *  The numerical setup and analytical solution are taken from
 *  Long Ju, Chunhua Zhang, and Zhaoli Guo. “Local reactive boundary scheme for
 *  irregular geometries in lattice Boltzmann method”. In: International Journal of Heat
 *  and Mass Transfer 150 (2020), p. 119314. doi: 10.1016/j.ijheatmasstransfer.2020.119314.
 *
 *  Numerical solution can be calculated using two approaches.
 *  Error norms are calculated for three norms.
 */


#include <olb.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

// Macro to switch between two approaches
#define BOUNDARY_CONDITION_APPROACH
//#define IN_BULK_APPROACH // IN_BULK_APPROACH only accurate for Peclet<0.1


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using TDESCRIPTOR = D3Q7<VELOCITY,OMEGA>;



using MOMENTA = momenta::AdvectionDiffusionBulkTuple;
using MomentaF = typename momenta::Tuple<
  momenta::SourcedDensity<typename MOMENTA::density>,
  momenta::FixedVelocityMomentum,
  typename MOMENTA::stress,
  typename MOMENTA::definition>;
using ReactionADBulkDynamics = dynamics::Tuple<T, TDESCRIPTOR, MomentaF, equilibria::FirstOrder, collision::BGK>;


const int N0 = 32;                // resolution
const int runs = 4;               // number of simulations with increasing resolution
const T errorTime = 1.;           // time at which errors are calculated in seconds

const T statIter = 0.04;          // time interval between lattice output

const T Pe = 1.;                  // Peclet number
const T physDiffusivity = 0.01;
const T physLength = 1.;          // physical domain length in x direction
const T relaxationTime = 0.53;
const T maxPhysT = 1.01;          // total simulation time in seconds

const T Ceq = 50.;                // equilibrium concentration at left wall

// vectors for error calculation
std::vector<T> L2Errors;
std::vector<T> LinfErrors;
std::vector<T> L1Errors;

/// functor for analytical solution
template<typename T, typename S>
class Solution : public AnalyticalF3D<T,S> {
public:
  T t;
  T u;
  Solution(AdeUnitConverter<T, TDESCRIPTOR> converter, std::size_t iT) : AnalyticalF3D<T,S>(1) {
    t = converter.getPhysTime(iT);
    u = converter.getCharPhysVelocity();
  };
  bool operator()(T output[], const S input[])
  {
    T D = physDiffusivity;
    T x = input[0];
    output[0] = Ceq*( 0.5*erfc((x-u*t)/(2*sqrt(D*t)))
                        + sqrt(u*u*t/(M_PI*D)) * exp(-pow(x-u*t, 2)/(4*D*t))
                        - 0.5*(1+u*x/D+u*u*t/D) * exp(u*x/D) * erfc((x+u*t)/(2*sqrt(D*t))) );
                        return true;
  }
};


void prepareGeometry(SuperGeometry <T,3> &superGeometry,
                     IndicatorF3D <T> &indicator,
                     AdeUnitConverter <T, TDESCRIPTOR> converter) {

  // === MATERIAL NUMBERS === //
  // solid: MN = 0            //
  // bulk: MN = 1             //
  // wall: MN = 2             //
  // inlet: MN = 3            //
  // outlet: MN = 4           //
  // bulk reaction: MN = 5    //
  // ======================== //

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  const T delta = converter.getPhysDeltaX();

  Vector<T,3> extend( 2*delta, physLength+2*delta, physLength+2*delta );
  Vector<T,3> origin(-delta, -delta, -delta);

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,0,0} );


  // Set material number for inflow
  origin[0] = -delta;
  IndicatorCuboid3D<T> inflow( extend, origin );
  superGeometry.rename( 2,3,inflow );

#ifdef IN_BULK_APPROACH
  IndicatorCuboid3D<T> layer( extend, origin );
  superGeometry.rename( 1,5,layer );
#endif

  // Set material number for outflow
  origin[0] = physLength-delta;
  IndicatorCuboid3D<T> outflow( extend, origin );
  superGeometry.rename( 2,4,outflow );

  superGeometry.communicate();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice(SuperLattice <T, TDESCRIPTOR> &ADlattice,
                    SuperGeometry <T,3> &superGeometry,
                    AdeUnitConverter <T, TDESCRIPTOR> converter) {
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

#ifdef BOUNDARY_CONDITION_APPROACH
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4});
  ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(bulkIndicator);
  // reactive Robin boundary at left side MN=3
  boundary::set<boundary::Robin>(ADlattice, superGeometry.getMaterialIndicator({3}));
  T reactionRate = converter.getCharLatticeVelocity(); //determines speed of reaction, here equal to velocity
  AnalyticalConst3D <T,T> coefficients(reactionRate, -converter.getLatticeDiffusivity(), reactionRate * Ceq); //set robin coefficients as found in the paper
  ADlattice.template defineField<descriptors::G>(superGeometry.getMaterialIndicator({3}), (coefficients)); //save coefficient in a field

  // ZeroGradient using RobinBoundary at right side MN=4
  boundary::set<boundary::Robin>(ADlattice, superGeometry.getMaterialIndicator({4}));
  AnalyticalConst3D <T,T> coefficients2(0., 1., 0.);
  ADlattice.template defineField<descriptors::G>(superGeometry.getMaterialIndicator({4}), (coefficients2));
#endif


#ifdef IN_BULK_APPROACH
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,5});
  ADlattice.defineDynamics<SourcedAdvectionDiffusionBGKdynamics>(bulkIndicator);
  // BounceBack for MN=3
  boundary::set<boundary::BounceBack><T, TDESCRIPTOR>(ADlattice, superGeometry.getMaterialIndicator({3}));
  // Zero Gradient at right side MN=4
  setZeroGradientBoundary<T, TDESCRIPTOR>(ADlattice, superGeometry.getMaterialIndicator({4}));
  // constants used in defineField
  AnalyticalConst3D<T,T> mat1(1.);
  AnalyticalConst3D<T,T> mat0(0.);
  // necessary for Zero Gradient Boundary
  ADlattice.template defineField<descriptors::SCALAR>(superGeometry.getMaterialIndicator({1,4}), mat1 );
  ADlattice.template defineField<descriptors::SCALAR>(superGeometry.getMaterialIndicator({0,3,5}), mat0 );
  // use field GAMMA to define where reaction is allowed to happen
  ADlattice.template defineField<descriptors::GAMMA>(superGeometry.getMaterialIndicator({5}), mat1 );
  ADlattice.template defineField<descriptors::GAMMA>(superGeometry.getMaterialIndicator({0,1,3,4}), mat0 );
#endif


  // Initial advection velocity and concentration
  AnalyticalConst3D <T, T> u0(converter.getCharLatticeVelocity(), 0.0, 0.0);
  AnalyticalConst3D<T,T> rho0(0.0);

  ADlattice.defineField<descriptors::VELOCITY>(bulkIndicator, u0);
  ADlattice.defineRho(bulkIndicator, rho0);
  ADlattice.iniEquilibrium(bulkIndicator, rho0, u0);

  ADlattice.setParameter<descriptors::OMEGA>( omega );
  ADlattice.initialize();
  clout << "Prepare Lattice ... OK" << std::endl;
}


// compute absolute and relative errors with different norms
void error(SuperLattice<T, TDESCRIPTOR>& ADlattice,
           SuperGeometry<T,3>& superGeometry,
           std::size_t iT,
           AdeUnitConverter <T, TDESCRIPTOR> converter) {
  OstreamManager clout(std::cout, "error");

  T result[3] = {T(), T(), T()};
  int tmp[] = {int()};

  SuperLatticeDensity3D <T, TDESCRIPTOR> numericalConc(ADlattice);
  Solution<T,T> analyticalConc(converter, iT);

#ifdef IN_BULK_APPROACH
  auto indicatorF = superGeometry.getMaterialIndicator({1,5}); // MN 3 and 4 should not be included for IN_BULK
#endif
#ifdef BOUNDARY_CONDITION_APPROACH
  auto indicatorF = superGeometry.getMaterialIndicator({1,3,4});
#endif

  SuperRelativeErrorL2Norm3D <T> relTemperatureErrorL2Norm(numericalConc, analyticalConc, indicatorF);
  SuperRelativeErrorL1Norm3D<T> relTemperatureErrorL1Norm(numericalConc, analyticalConc, indicatorF);
  SuperRelativeErrorLinfNorm3D<T> relTemperatureErrorLinfNorm(numericalConc, analyticalConc, indicatorF);

  /// Console output of every error norm
  relTemperatureErrorL2Norm(result, tmp);
  L2Errors.push_back(result[0]);
  relTemperatureErrorL1Norm(result, tmp);
  L1Errors.push_back(result[0]);
  relTemperatureErrorLinfNorm(result, tmp);
  LinfErrors.push_back(result[0]);
}



void getResults(SuperLattice<T, TDESCRIPTOR> &ADlattice,
                std::size_t iT,
                AdeUnitConverter <T, TDESCRIPTOR> converter,
                SuperGeometry<T,3> &superGeometry,
                util::Timer <T> &timer) {

  OstreamManager clout(std::cout, "getResults");
  SuperVTMwriter3D <T> vtkWriter("longitudinalMixing3d");

  SuperLatticeDensity3D <T, TDESCRIPTOR> concentration(ADlattice);
  vtkWriter.addFunctor( concentration );
  concentration.getName() = "concentration";

  Solution<T,T> analyticalSolution(converter, iT);
  SuperLatticeFfromAnalyticalF3D<T, TDESCRIPTOR> functor(analyticalSolution, ADlattice);
  functor.getName() = "analyticalSolution";
  vtkWriter.addFunctor( functor );


  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization

    SuperLatticeCuboid3D <T, TDESCRIPTOR> cuboid(ADlattice);
    SuperLatticeRank3D <T, TDESCRIPTOR> rank(ADlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
    vtkWriter.write(iT);

    timer.update(iT);
    timer.printStep();
    ADlattice.getStatistics().print(iT, converter.getPhysTime(iT));
  } else if (iT % converter.getLatticeTime(statIter) == 0) {
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    /// ADLattice statistics console output
    timer.update(iT);
    timer.printStep();
    ADlattice.getStatistics().print(iT, converter.getPhysTime(iT));
    vtkWriter.write(iT);
  }
  if(iT == converter.getLatticeTime(errorTime)){
    error(ADlattice, superGeometry, iT, converter);
  }
}


void simulate(int N) {

  OstreamManager clout(std::cout, "simulate");
  clout << "Executing the simulation with N=" << std::to_string(N) << std::endl;

  AdeUnitConverterFromResolutionAndRelaxationTime<T,TDESCRIPTOR> converter(
      N,                              // Resolution
      relaxationTime,                 // Relaxation Time
      physLength,                     // charPhysLength
      Pe*physDiffusivity/physLength,  // charPhysVelocity from peclet
      physDiffusivity,                // physDiffusivity
      1                               // physDensity
  );
  converter.print();

  // change directory for different runs
  if (runs > 1) {
    singleton::directories().setOutputDir("./tmp/N" + std::to_string(N) + "/");
  }

  /// === 2nd Step: Prepare Geometry ===
  Vector<T,3> extend( physLength, converter.getPhysDeltaX(), converter.getPhysDeltaX() ); //minimal length in y,z direction
  Vector<T,3> origin( 0.0, 0.0, 0.0 );

  IndicatorCuboid3D <T> cuboid(extend, origin);

  /// Instantiation of an empty cuboidDecomposition
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidDecomposition3D <T> cuboidDecomposition(cuboid, converter.getPhysDeltaX(), noOfCuboids);
  cuboidDecomposition.setPeriodicity({false, true, true});

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer <T> loadBalancer(cuboidDecomposition);

  /// Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, 2);
  prepareGeometry(superGeometry, cuboid, converter);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, TDESCRIPTOR> ADlattice(superGeometry);
  prepareLattice(ADlattice, superGeometry, converter);

  // Prepare reaction at left wall for IN_BULK_APPROACH
  SuperLatticeCoupling Reaction(
      LongitudinalMixingReactionCoupling<T>{},
      names::Concentration0{}, ADlattice);
  Reaction.template setParameter<typename LongitudinalMixingReactionCoupling<T>::
    REACTION_CONSTANT>({converter.getCharLatticeVelocity()}); //determines speed of reaction, here equal to velocity
  Reaction.template setParameter<typename LongitudinalMixingReactionCoupling<T>::EQUILIBRIUM>({Ceq}); //equilibrium concentration


  /// === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {

    /// === 6th Step: Collide and Stream Execution ===
#ifdef IN_BULK_APPROACH
    Reaction.execute();
#endif
    ADlattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( ADlattice, iT, converter, superGeometry, timer );
  }
  timer.stop();
  timer.printSummary();
}

int main(int argc, char *argv[]) {
  OstreamManager clout(std::cout, "main");
  initialize(&argc, &argv);

  for (int i = 0; i < runs; ++i) {
    simulate(pow(2, i) * N0);
  }
  for (int i = 0; i < runs; ++i) {
    clout << "Errors with N="<<(pow(2, i) * N0)<<" after "<<errorTime<<"s:   "<<"L1:"<<std::log(L1Errors[i])
          <<" L2:"<<std::log(L2Errors[i])<<" Linf:"<<std::log(LinfErrors[i]) << std::endl;
  }
}

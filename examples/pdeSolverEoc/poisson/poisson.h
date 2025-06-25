/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Fedor Bukreev
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

#include <olb.h>

 /*
  * the expressions used in this file are used for validational part of the paper
  * F. Bukreev, A. Kummerländer, J. Jeßberger, D.Teutscher, S. Ito, S. Simonis, H. Nirschl, M. J. Krause,
  * Hydro-electrochemical Saturation in Nano-scale Porous Media. Part I: Sensitivity Assessing Simulation Approach,
  * submitted to Chemical Engineering Science 2024
  *
  * In current example, Poisson equation with one-dimensional analytical solution is discretized with LBM.
  * This equation is often used for electric potential distribution. The analytical solution for this case is taken from
  * Z. Chai, B. Shi, A novel lattice boltzmann model for the poisson equation, Applied Mathematical Modelling 32 (10) (2008) 2050–2058
  * doi:https://doi.org/10.1016/j.apm.2007.06.033.
  */

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;
using DESCRIPTOR = D3Q7<VELOCITY,SOURCE>;
using BulkDynamics = SourcedAdvectionDiffusionBGKdynamics<T,DESCRIPTOR>;

const T length = 1.;
const T width = 0.25;
int N = 20;
const T maxPhysT = 10.;
const T constK = 27.79;

// the relative error norms
T psiL1RelError = 0;
T psiL2RelError = 0;
T psiLinfRelError = 0;

// Anayltical solution for potential profile
template <typename T, typename S, typename DESCRIPTOR>
    class PotentialProfile1D : public AnalyticalF3D<T, S>
    {
    private:
        T x0, x1, k;
        UnitConverter<T, DESCRIPTOR> const &converter;

    public:
        PotentialProfile1D(T x0_, T x1_, T k_, UnitConverter<T, DESCRIPTOR> const &converter_) : AnalyticalF3D<T, S>(1),
                                                                                                                        x0(x0_), x1(x1_), k(k_), converter(converter_)
        {
            this->getName() = "PotentialProfile1D";
        };

        bool operator()(T output[1], const S x[3])
        {
            T distX = (x[0] - x0) / (x1 - x0);
            T L = converter.getPhysDeltaX();
            output[0] = (util::exp(k) - 1.)/(util::exp(k) - util::exp(-k))*util::exp(-k*distX) + (1. - util::exp(-k))/(util::exp(k) - util::exp(-k))*util::exp(k*distX);
            if(x[0] <= L/2. || x[0] >= x1 - L/2.) output[0] = 1.;

            return true;
        };
    };

// PostProcessor for updating of source term in each time cell in each cell
struct ComputeSourceTerm {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct C_T : public descriptors::FIELD_BASE<1> { };
  struct CONSTK : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<C_T,CONSTK>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    V psi = cell.computeRho();
    V k = parameters.template get<CONSTK>();
    V source = -k*k*psi*parameters.template get<C_T>();
    cell.template setField<descriptors::SOURCE>(source);
  }
};


void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,0,0} );
  superGeometry.clean();

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1,2 -->bulk dynamics
  sLattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator({1, 2}));

  // Material=2 -->Dirichlet boundary
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLattice, superGeometry, 2);

  AnalyticalConst3D<T,T> rhoF( 1. );
  AnalyticalConst3D<T,T> rho0( 0 );
  std::vector<T> velocity( 2,T() );
  AnalyticalConst3D<T,T> uF( velocity );
  sLattice.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({1, 2}),uF);

  sLattice.iniEquilibrium( superGeometry, 2, rhoF, uF );
  sLattice.defineRhoU( superGeometry, 2, rhoF, uF );

  sLattice.iniEquilibrium( superGeometry, 1, rho0, uF );
  sLattice.defineRhoU( superGeometry, 1, rho0, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void error( SuperLattice<T, DESCRIPTOR>& sLattice,
            UnitConverter<T, DESCRIPTOR> const& converter,
            SuperGeometry<T,3>& superGeometry, bool var)
{
  OstreamManager clout( std::cout,"Errors" );
  if( var ){
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    int tmp[]= { };
    T result[2] = { };

    T L = converter.getPhysDeltaX();
    PotentialProfile1D<T,T,DESCRIPTOR> psiSol(0.*L, length - 0.*L, constK, converter);
    SuperLatticeDensity3D<T, DESCRIPTOR> density( sLattice );

    auto material = superGeometry.getMaterialIndicator(1);
    SuperRelativeErrorL1Norm3D<T>   errorPsiL1Norm(density, psiSol, *material);
    SuperRelativeErrorL2Norm3D<T>   errorPsiL2Norm(density, psiSol, *material);
    SuperRelativeErrorLinfNorm3D<T> errorPsiLinfNorm(density, psiSol, *material);

    errorPsiL1Norm(result,tmp);
    clout << "Psi-L1-error: " << result[0] << std::endl;
    psiL1RelError = result[0];

    errorPsiL2Norm(result,tmp);
    clout << "Psi-L2-error: " << result[0] << std::endl;
    psiL2RelError = result[0];

    errorPsiLinfNorm(result,tmp);
    clout << "Psi-Linf-error: " << result[0] << std::endl;
    psiLinfRelError = result[0];
  }
}


void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                 bool hasConverged, Gnuplot<T>& gplot )
{

  OstreamManager clout( std::cout,"getResults" );

  const int vtkIter  = converter.getLatticeTime( .1 );
  const int statIter = converter.getLatticeTime( .1 );

  T L = converter.getPhysDeltaX();
  PotentialProfile1D<T,T,DESCRIPTOR> psiSol(0.*L, length - 0.*L, 27.79, converter);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalSol(psiSol, sLattice);
  analyticalSol.getName() = "analytical solution";

  SuperVTMwriter3D<T> vtmWriter( "poisson1d" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter == 0 ) {

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  // Writes the VTK files
  if ( iT%vtkIter == 0 || hasConverged) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticeDensity3D<T,DESCRIPTOR> potential( sLattice );
    potential.getName() = "potential";

    vtmWriter.addFunctor( potential );
    vtmWriter.addFunctor( analyticalSol );
    vtmWriter.write( iT );
  }

  if(hasConverged){
      timer.update( iT );
      timer.printStep();
      sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  if(hasConverged){
    gplot.setData (
          T(converter.getResolution()),
          { psiL1RelError, psiL2RelError, psiLinfRelError},
          { "psi L1 rel Error","psi L2 rel Error",
            "psi Linf rel error"},
          "top right",
          { 'p','p','p' } );
    clout << "Time Averaged psi-L1-error: " << psiL1RelError << std::endl;
    clout << "Time Averaged psi-L2-error: " << psiL2RelError << std::endl;
    clout << "Time Averaged psi-Linf-error: " << psiLinfRelError << std::endl;
    psiL1RelError = 0;
    psiL2RelError = 0;
    psiLinfRelError = 0;
  }
}



void simulatePoisson3D(int N, Gnuplot<T>& gplot)
{

  // === 1st Step: Initialization ===
  OstreamManager clout( std::cout,"main" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    (int) N,                        // resolution: number of voxels per charPhysL
    (T)   1.,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   width,       // charPhysLength: reference length of simulation geometry
    (T)   1.0,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.0, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                       // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("poisson1d");

  // === 2nd Step: Prepare Geometry ===
  T L = converter.getPhysDeltaX();
  Vector<T,3> extend( length,width, 3*L );
  Vector<T,3> origin( 0,0,0 );
  IndicatorCuboid3D<T> cuboid( extend, origin );

#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize() );
#else
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, converter.getPhysDeltaX(), 1 );
#endif

  cuboidDecomposition.print();

  // This is a one-dimensional case computed on 3D lattice, which is why 2 unneeded directions are set periodic
  cuboidDecomposition.setPeriodicity({false,true,true});

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,3> superGeometry( cuboidDecomposition, loadBalancer );
  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  prepareLattice( converter, sLattice, superGeometry );

  sLattice.addPostProcessor<stage::PreCollide>(meta::id<ComputeSourceTerm>{});
  sLattice.setParameter<ComputeSourceTerm::C_T>(converter.getConversionFactorTime());
  sLattice.setParameter<ComputeSourceTerm::CONSTK>(constK);

  // === 4th Step: Main Loop with Timer ===

  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( 100, 1.e-9 );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {

    sLattice.collideAndStream();

    error( sLattice, converter, superGeometry, converge.hasConverged());
    getResults( sLattice, converter, iT, superGeometry, timer, converge.hasConverged(), gplot );

    if(converge.hasConverged()) { break; }

    converge.takeValue( sLattice.getStatistics().getAverageRho(), false );
  }

  timer.stop();
  timer.printSummary();
}
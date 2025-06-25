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

 /*
  * the expressions used in this file are used for validational part of the paper
  * F. Bukreev, A. Kummerländer, J. Jeßberger, D.Teutscher, S. Ito, S. Simonis, H. Nirschl, M. J. Krause,
  * Hydro-electrochemical Saturation in Nano-scale Porous Media. Part I: Sensitivity Assessing Simulation Approach,
  * submitted to Chemical Engineering Science 2024
  *
  * The current example solves the Poisson-Boltzmann equation describing Gouy-Chapman double electric layer model for one dimension.
  * Example contains 3 lattices:
  * 1 for Poisson equation,
  * 1 for Nernst-Planck equation of cation,
  * 1 for Nernst-Planck equation of anion,
  */

#include <olb.h>


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;
using DESCRIPTOR = D3Q19<VELOCITY,SOURCE>;
using BulkDynamicsPsi = SourcedAdvectionDiffusionBGKdynamics<T,DESCRIPTOR>;
using BulkDynamicsConc = AdvectionDiffusionBGKdynamics<T,DESCRIPTOR>;

const T maxPhysT = 10.;// maximal physical temperature
const T dielectricC = 6.95e-10; //  C/V/m dielectric constant: extent to which a material holds or concentrates electric flux
const T Faraday = 96485.33;  // C/mol  Faraday constant: quotient of the total electric charge (q) by the amount (n) of elementary charge carriers in any given sample of matter
const T valence = 1.;// valence: electrons in the outermost shell of an atom
const T diffusion = 1.e-9;// diffusion: process resulting from random motion of molecules by which there is a net flow of matter from a region of high concentration to a region of low concentration
const T temperature = 298.15;// temperature in Kelvin
const T Boltzmann = 1.38065e-23; // VC/K the Bolzmann constant: the proportionality factor that relates the average relative thermal energy of particles in a gas with the thermodynamic temperature of the gas
const T charge = 1.602177e-19; // C
const T NA = 6.02e23;// Avogadro number: number of constituent particles in a sample
const T psi0 = -0.02; // V electric potential
const T C0 = 0.01; // mol/L  molar concentration
const T Debye = util::sqrt(dielectricC*Boltzmann*temperature/2./charge/charge/valence/valence/C0/NA);// the debye length is a measure of the charge carriers electrostatic effect in a solution and how far the electostatic effect persists
const T width = Debye*13.;// distance between the outer solid bodies and the domain boundaries
const T length = width/4.;

//the relative error norms
T psiL1RelError = 0;
T psiL2RelError = 0;
T psiLinfRelError = 0;

T concL1RelError = 0;
T concL2RelError = 0;
T concLinfRelError = 0;

T conc2L1RelError = 0;
T conc2L2RelError = 0;
T conc2LinfRelError = 0;

// Analytical Profile for electric potential
template <typename T, typename S, typename DESCRIPTOR>
 class PotentialProfile1D : public AnalyticalF3D<T, S>
    {
    private:
        T y0, y1;
        UnitConverter<T, DESCRIPTOR> const &converter;
    public:
        PotentialProfile1D(T y0_, T y1_, UnitConverter<T, DESCRIPTOR> const &converter_) : AnalyticalF3D<T, S>(1),
                                                                          y0(y0_), y1(y1_), converter(converter_)
        {
            this->getName() = "PotentialProfile1D";
        };

        bool operator()(T output[1], const S x[3])
        {
            T distY = x[1] - y0;
            output[0] = psi0 * util::exp( -distY/Debye );
            return true;
        };
    };

// Analytical Profile for ion concentrations
template <typename T, typename S, typename DESCRIPTOR>
    class ConcentrationProfile1D : public AnalyticalF3D<T, S>
    {
    private:
        T valence;
        PotentialProfile1D<T,T,DESCRIPTOR> &potential;

    public:
        ConcentrationProfile1D(T valence_, PotentialProfile1D<T,T,DESCRIPTOR> &potential_) : AnalyticalF3D<T, S>(1),
                                                                          valence(valence_), potential(potential_)
        {
            this->getName() = "ConcentrationProfile1D";
        };

        bool operator()(T output[1], const S x[3])
        {
            T psi[1];
            potential(psi,x);
            output[0] = C0 * util::exp(-charge*valence*psi[0]/Boltzmann/temperature);
            return true;
        };
    };

void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{0,1,0} );
  superGeometry.clean();

  T eps = converter.getPhysDeltaX();
  Vector<T,3> extend(length, 2*eps, 5*eps);
  Vector<T,3> origin;
  IndicatorCuboid3D<T> inlet( extend, origin );
  superGeometry.rename( 2,3,1,inlet );

  origin[1] = width-2*eps;
  IndicatorCuboid3D<T> outlet( extend, origin );
  superGeometry.rename( 2,4,1,outlet );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLatticePoisson( UnitConverter<T,DESCRIPTOR> const& converter,
                            SuperLattice<T, DESCRIPTOR>& sLattice,
                            SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLatticePoisson" );
  //clout << "Prepare Lattice Poisson..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1,3,4 -->bulk dynamics
  sLattice.defineDynamics<BulkDynamicsPsi>(superGeometry.getMaterialIndicator({1, 3, 4}));

  // Material=3,4 -->Dirichlet boundary for electric potential
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLattice, superGeometry, 3);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLattice, superGeometry, 4);

  AnalyticalConst3D<T,T> rhoF( psi0 );
  AnalyticalConst3D<T,T> rho0( 0 );
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> uF( velocity );
  sLattice.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({0, 1, 3, 4}),uF);
  sLattice.defineField<descriptors::SOURCE>(superGeometry.getMaterialIndicator({0, 1, 3, 4}),rho0);

  sLattice.iniEquilibrium( superGeometry, 3, rhoF, uF );
  sLattice.defineRhoU( superGeometry, 3, rhoF, uF );

  auto bulkIndicator = superGeometry.getMaterialIndicator({0, 1, 4});
  sLattice.iniEquilibrium( bulkIndicator, rho0, uF );
  sLattice.defineRhoU( bulkIndicator, rho0, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  sLattice.initialize();
}

void prepareLatticeNernstPlanck( UnitConverter<T,DESCRIPTOR> const& converter,
                                 SuperLattice<T, DESCRIPTOR>& sLattice,
                                 SuperGeometry<T,3>& superGeometry,
                                 T CB, T sign )
{

  OstreamManager clout( std::cout,"prepareLatticeNernstPlanck" );
  clout << "Prepare Lattice Nernst-Planck..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1, 4 -->bulk dynamics
  sLattice.defineDynamics<BulkDynamicsConc>(superGeometry.getMaterialIndicator({1, 4}));

  // Material=3 -->Bounce Back (Zero Gradient)
  sLattice.defineDynamics<BounceBack>(superGeometry, 3);

  // Material=4 -->Dirichlet boundary for concentrations
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLattice, superGeometry, 4);

  AnalyticalConst3D<T,T> rhoF( C0 );
  AnalyticalConst3D<T,T> rho0( 0 );
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> uF( velocity );
  sLattice.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({0, 1, 3, 4}),uF);

  AnalyticalConst3D<T,T> rhoB( CB );
  clout << "Concentration at boundary: " << CB << std::endl;

  sLattice.iniEquilibrium( superGeometry, 4, rhoF, uF );
  sLattice.defineRhoU( superGeometry, 4, rhoF, uF );

  auto bulkIndicator1 = superGeometry.getMaterialIndicator({0,1,3});
  sLattice.iniEquilibrium( bulkIndicator1, rho0, uF );
  sLattice.defineRhoU( bulkIndicator1, rho0, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  sLattice.initialize();

  clout << "Prepare Lattice Nernst-Planck... OK" << std::endl;
}

void error( SuperLattice<T, DESCRIPTOR>& sLatticePoisson,
            SuperLattice<T, DESCRIPTOR>& sLatticeCation,
            SuperLattice<T, DESCRIPTOR>& sLatticeAnion,
            UnitConverter<T, DESCRIPTOR> const& converter,
            SuperGeometry<T,3>& superGeometry, bool var)
{
  OstreamManager clout( std::cout,"Errors" );
  if( var ){
    sLatticePoisson.setProcessingContext(ProcessingContext::Evaluation);
    sLatticeCation.setProcessingContext(ProcessingContext::Evaluation);
    sLatticeAnion.setProcessingContext(ProcessingContext::Evaluation);
    int tmp[]= { };
    T result[2] = { };

    T eps = converter.getPhysDeltaX();
    PotentialProfile1D<T,T,DESCRIPTOR> psiSol(0.*eps, width - 0.*eps, converter);
    ConcentrationProfile1D<T,T,DESCRIPTOR> concSolCation(valence, psiSol);
    ConcentrationProfile1D<T,T,DESCRIPTOR> concSolAnion(-valence, psiSol);
    SuperLatticeDensity3D<T, DESCRIPTOR> psi( sLatticePoisson );
    SuperLatticeDensity3D<T, DESCRIPTOR> cation( sLatticeCation );
    SuperLatticeDensity3D<T, DESCRIPTOR> anion( sLatticeAnion );

    auto material = superGeometry.getMaterialIndicator(1);
    SuperRelativeErrorL1Norm3D<T>   errorPsiL1Norm(psi, psiSol, *material);
    SuperRelativeErrorL2Norm3D<T>   errorPsiL2Norm(psi, psiSol, *material);
    SuperRelativeErrorLinfNorm3D<T> errorPsiLinfNorm(psi, psiSol, *material);

    SuperRelativeErrorL1Norm3D<T>   errorConcL1Norm(cation, concSolCation, *material);
    SuperRelativeErrorL2Norm3D<T>   errorConcL2Norm(cation, concSolCation, *material);
    SuperRelativeErrorLinfNorm3D<T> errorConcLinfNorm(cation, concSolCation, *material);

    SuperRelativeErrorL1Norm3D<T>   errorConc2L1Norm(anion, concSolAnion, *material);
    SuperRelativeErrorL2Norm3D<T>   errorConc2L2Norm(anion, concSolAnion, *material);
    SuperRelativeErrorLinfNorm3D<T> errorConc2LinfNorm(anion, concSolAnion, *material);

    errorPsiL1Norm(result,tmp);
    clout << "Relative Potential-L1-error: " << result[0] << std::endl;
    psiL1RelError = result[0];

    errorPsiL2Norm(result,tmp);
    clout << "Relative Potential-L2-error: " << result[0] << std::endl;
    psiL2RelError = result[0];

    errorPsiLinfNorm(result,tmp);
    clout << "Relative Potential-Linf-error: " << result[0] << std::endl;
    psiLinfRelError = result[0];

    errorConcL1Norm(result,tmp);
    clout << "Relative Cation Conc-L1-error: " << result[0] << std::endl;
    concL1RelError = result[0];

    errorConcL2Norm(result,tmp);
    clout << "Relative Cation Conc-L2-error: " << result[0] << std::endl;
    concL2RelError = result[0];

    errorConcLinfNorm(result,tmp);
    clout << "Relative Cation Conc-Linf-error: " << result[0] << std::endl;
    concLinfRelError = result[0];

    errorConc2L1Norm(result,tmp);
    clout << "Relative Anion Conc-L1-error: " << result[0] << std::endl;
    conc2L1RelError = result[0];

    errorConc2L2Norm(result,tmp);
    clout << "Relative Anion Conc-L2-error: " << result[0] << std::endl;
    conc2L2RelError = result[0];

    errorConc2LinfNorm(result,tmp);
    clout << "Relative Anion Conc-Linf-error: " << result[0] << std::endl;
    conc2LinfRelError = result[0];
  }
}

void getResultsNernstPlanck( SuperLattice<T, DESCRIPTOR>& sLatticeCation, SuperLattice<T, DESCRIPTOR>& sLatticeAnion,
                             SuperLattice<T, DESCRIPTOR>& sLatticePoisson,
                             UnitConverter<T, DESCRIPTOR> const& converter, std::size_t iT,
                             SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                             bool hasConverged, Gnuplot<T>& gplot)
{

  OstreamManager clout( std::cout,"getResultsNernstPlanck" );

  const int vtkIter  = converter.getLatticeTime( .5 );
  const int statIter = converter.getLatticeTime( .5 );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperVTMwriter3D<T> vtmWriterCation( "cation" );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLatticeCation );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLatticeCation );

    vtmWriterCation.write( cuboid );
    vtmWriterCation.write( rank );
    vtmWriterCation.createMasterFile();

    SuperVTMwriter3D<T> vtmWriterPoisson( "poisson" );
    vtmWriterPoisson.createMasterFile();

    SuperVTMwriter3D<T> vtmWriterAnion( "anion" );
    vtmWriterAnion.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter == 0 || hasConverged ) {

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLatticeCation.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  // Writes the VTK files
  if ( iT%vtkIter == 0 || hasConverged) {
    sLatticeCation.setProcessingContext(ProcessingContext::Evaluation);
    sLatticeAnion.setProcessingContext(ProcessingContext::Evaluation);
    sLatticePoisson.setProcessingContext(ProcessingContext::Evaluation);

    sLatticePoisson.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterPoisson( "poisson" );
      SuperLatticeDensity3D psi( sLatticePoisson);
      psi.getName() = "psi";
      T eps = converter.getPhysDeltaX();
      PotentialProfile1D<T,T,DESCRIPTOR> psiSol(0.*eps, width - 0.*eps, converter);
      SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalPsi(psiSol, sLatticePoisson);
      analyticalPsi.getName() = "analytical potential";
      vtmWriterPoisson.addFunctor(psi);
      vtmWriterPoisson.addFunctor(analyticalPsi);
      task(vtmWriterPoisson, iT);
    });

    sLatticeCation.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterCation( "cation" );
      SuperLatticeDensity3D<T,DESCRIPTOR> cation( sLatticeCation );
      cation.getName() = "cation";
      SuperLatticePhysField3D<T,DESCRIPTOR,VELOCITY> velCation( sLatticeCation, converter.getConversionFactorVelocity() );
      velCation.getName() = "velCation";
      T eps = converter.getPhysDeltaX();
      PotentialProfile1D<T,T,DESCRIPTOR> psiSol(0.*eps, width - 0.*eps, converter);
      ConcentrationProfile1D<T,T,DESCRIPTOR> concSol(valence, psiSol);
      SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalConcCation(concSol, sLatticeCation);
      analyticalConcCation.getName() = "analytical concentration cation";
      vtmWriterCation.addFunctor(cation);
      vtmWriterCation.addFunctor(velCation);
      vtmWriterCation.addFunctor(analyticalConcCation);
      task(vtmWriterCation, iT);
    });

    sLatticeAnion.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterAnion( "anion" );
      SuperLatticeDensity3D<T,DESCRIPTOR> anion( sLatticeAnion );
      anion.getName() = "anion";
      SuperLatticePhysField3D<T,DESCRIPTOR,VELOCITY> velAnion( sLatticeAnion, converter.getConversionFactorVelocity() );
      velAnion.getName() = "velAnion";
      T eps = converter.getPhysDeltaX();
      PotentialProfile1D<T,T,DESCRIPTOR> psiSol(0.*eps, width - 0.*eps, converter);
      ConcentrationProfile1D<T,T,DESCRIPTOR> concSol2(-valence, psiSol);
      SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalConcAnion(concSol2, sLatticeAnion);
      analyticalConcAnion.getName() = "analytical concentration anion";
      vtmWriterAnion.addFunctor(anion);
      vtmWriterAnion.addFunctor(velAnion);
      vtmWriterAnion.addFunctor(analyticalConcAnion);
      task(vtmWriterAnion, iT);
    });
   }

  if(hasConverged){
    gplot.setData (
          T(converter.getResolution()),
          { psiL1RelError, psiL2RelError, psiLinfRelError,
            concL1RelError, concL2RelError, concLinfRelError,
            conc2L1RelError, conc2L2RelError, conc2LinfRelError},
          { "potential L1 Rel Error","potential L2 Rel Error",
            "potential Linf Rel error","[cation] L1 Rel Error","[cation] L2 Rel Error",
            "[cation] Linf Rel error","[anion] L1 Rel Error","[anion] L2 Rel Error",
            "[anion] Linf Rel error"},
          "top right",
          { 'p','p','p','p','p','p','p','p','p' } );
    psiL1RelError = 0;
    psiL2RelError = 0;
    psiLinfRelError = 0;
    concL1RelError = 0;
    concL2RelError = 0;
    concLinfRelError = 0;
    conc2L1RelError = 0;
    conc2L2RelError = 0;
    conc2LinfRelError = 0;
  }
}

void simulatePoisson(SuperLattice<T, DESCRIPTOR>& sLattice,
                     UnitConverter<T, DESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry,
                     bool poissonLoop)
{

  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( 10, 1.e-9 );
  timer.start();

  for ( std::size_t iT = 0; iT < (!poissonLoop)*1 + (poissonLoop)*converter.getLatticeTime( maxPhysT ); ++iT ) {

    sLattice.collideAndStream();

    if(converge.hasConverged()) { break; }

    converge.takeValue( sLattice.getStatistics().getAverageRho(), false );
  }

  timer.stop();
}

void simulatePoissonNernstPlanck3D(int N, bool poissonLoop, Gnuplot<T>& gplot)
{

  // === 1st Step: Initialization ===
  OstreamManager clout( std::cout,"main" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converterPoisson(
    (int) N,                        // resolution: number of voxels per charPhysL
    (T)   0.9,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   width,       // charPhysLength: reference length of simulation geometry
    (T)   0.001,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.0, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                       // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converterPoisson.print();
  // Writes the converter log in a file
  converterPoisson.write("poisson");

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converterNernstPlanck(
    (int) N,                        // resolution: number of voxels per charPhysL
    (T)   0.9,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   width,       // charPhysLength: reference length of simulation geometry
    (T)   0.001,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   diffusion, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                       // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converterNernstPlanck.print();
  // Writes the converter log in a file
  converterNernstPlanck.write("nernstPlanck");

  // === 2nd Step: Prepare Geometry ===
  T L = converterNernstPlanck.getPhysDeltaX();
  Vector<T,3> extend( length,width, 5*L );
  Vector<T,3> origin( 0,0,0 );
  IndicatorCuboid3D<T> cuboid( extend, origin );

#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, converterNernstPlanck.getPhysDeltaX(), singleton::mpi().getSize() );
#else
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, converterNernstPlanck.getPhysDeltaX(), 1 );
#endif

  cuboidDecomposition.print();
  cuboidDecomposition.setPeriodicity({true,false,true});

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,3> superGeometry( cuboidDecomposition, loadBalancer );
  prepareGeometry( converterNernstPlanck, superGeometry );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice<T, DESCRIPTOR> sLatticeCation( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLatticePoisson( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLatticeAnion( superGeometry );

  T CationWall = C0 * util::exp(-charge*valence*psi0/Boltzmann/temperature);
  T AnionWall = C0 * util::exp(charge*valence*psi0/Boltzmann/temperature);

  prepareLatticeNernstPlanck( converterNernstPlanck, sLatticeCation, superGeometry, CationWall, 1. );
  prepareLatticeNernstPlanck( converterNernstPlanck, sLatticeAnion, superGeometry, AnionWall, -1. );
  prepareLatticePoisson( converterPoisson, sLatticePoisson, superGeometry );

  T npVelCoeff = charge * valence * diffusion / Boltzmann / temperature / converterNernstPlanck.getConversionFactorVelocity();
  T sourceCoeff = 1./dielectricC * Faraday * valence * converterPoisson.getConversionFactorTime();

  clout << "Debye length [m]: " << Debye << std::endl;

  SuperLatticeCoupling coupling(
    PNPCoupling<T>{},
    names::Concentration0{}, sLatticeCation,
    names::Concentration1{}, sLatticePoisson,
    names::Concentration2{}, sLatticeAnion);
  coupling.setParameter<PNPCoupling<T>::DX>(converterPoisson.getPhysDeltaX());
  coupling.setParameter<PNPCoupling<T>::NPVELCOEFF>(npVelCoeff);
  coupling.setParameter<PNPCoupling<T>::POISSONCOEFF>(sourceCoeff);
  coupling.setParameter<PNPCoupling<T>::OMEGA>(converterPoisson.getLatticeRelaxationFrequency());
  coupling.restrictTo(superGeometry.getMaterialIndicator({1}));

  // === 4th Step: Main Loop with Timer ===

  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converterNernstPlanck.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( 500, 1.e-9 );
  timer.start();

  for ( std::size_t iT = 0; iT < converterNernstPlanck.getLatticeTime( maxPhysT ); ++iT ) {
    // internal Poisson loop. In terms of efficiency choose an appropriate frequency of Poisson equation updates
    int poissonLoopPerTotalSteps = 5;
    if( iT == 0 || iT%((!poissonLoop)*1 + (poissonLoop)*poissonLoopPerTotalSteps) == 0){
      simulatePoisson(sLatticePoisson, converterPoisson, superGeometry, poissonLoop);
    }

    coupling.execute();

    sLatticeCation.collideAndStream();
    sLatticeAnion.collideAndStream();

    error( sLatticePoisson, sLatticeCation, sLatticeAnion, converterNernstPlanck, superGeometry, converge.hasConverged());
    getResultsNernstPlanck( sLatticeCation, sLatticeAnion, sLatticePoisson, converterNernstPlanck, iT, superGeometry, timer, converge.hasConverged(), gplot);

    if(converge.hasConverged()) { break; }
    converge.takeValue( sLatticePoisson.getStatistics().getAverageRho(), false );
  }

  timer.stop();
  timer.printSummary();
  singleton::pool().wait();
}

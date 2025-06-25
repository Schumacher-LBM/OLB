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
 *  GNU General Public License for mor details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this proquotiegram; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

 /*
  * the expressions used in this file are used for validational part of the paper
  * F. Bukreev, A. Kummerländer, J. Jeßberger, D.Teutscher, S. Ito, S. Simonis, H. Nirschl, M. J. Krause,
  * Hydro-electrochemical Saturation in Nano-scale Porous Media. Part I: Sensitivity Assessing Simulation Approach,
  * submitted to Chemical Engineering Science 2024
  *
  * The current example solves the electroosmosis 1D problem
  *[Lattice Boltzmann Simulation of Electroosmotic Flows in Micro- and Nanochannels; Fuzhi Tian1, Baoming Li1,2, and Daniel Y. Kwok1]
  * Example contains 5 lattices:
  * 1 for Poisson equation,
  * 1 for Nernst-Planck equation of cation,
  * 1 for Nernst-Planck equation of anion,
  * 1 for Navier-Stokes equations of the carrier fluid
  */
#ifndef EXAMPLES_ELECTROOSMOSIS_H
#define EXAMPLES_ELECTROOSMOSIS_H

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;
using DESCRIPTOR = D3Q19<VELOCITY,SOURCE>;
using BulkDynamicsPsi = SourcedAdvectionDiffusionBGKdynamics<T,DESCRIPTOR>;
using BulkDynamicsConc = AdvectionDiffusionBGKdynamics<T,DESCRIPTOR>;

using DESCRIPTORNSE = D3Q19<FORCE>;
using BulkDynamicsNSE = ForcedBGKdynamics<T,DESCRIPTORNSE>;

const T maxPhysT = 10.;// maximal physical time
const T dielectricC = 6.95e-10; //  C/V/m dielectric constant: extent to which a material holds or concentrates electric flux
const T Faraday = 9.649e4;  // C/mol Faraday constant: quotient of the total electric charge (q) by the amount (n) of elementary charge carriers in any given sample of matter
const T valence = 1.;// valence: electrons in the outermost shell of an atom
const T viscosity = 1.e-6;// viscosity: a fluids resistance to deformation
const T diffusion = 1.e-8;// diffusion: process resulting from random motion of molecules by which there is a net flow of matter from a region of high concentration to a region of low concentration
const T temperature = 293.15;// temperature in Kelvin
const T Boltzmann = 1.38065e-23; // VC/K the Bolzmann constant: the proportionality factor that relates the average relative thermal energy of particles in a gas with the thermodynamic temperature of the gas
const T charge = 1.602177e-19; // C
const T NA = 6.02e23;// Avogadro number: number of constituent particles in a sample
const T psi0 = -0.02; // V electric potential
const T C0 = 0.01; // m³/mol molar volume
const T Debye = util::sqrt(dielectricC*Boltzmann*temperature/2./charge/charge/valence/valence/C0/NA);// the debye length is a measure of the charge carriers electrostatic effect in a solution and how far the electrostatic effect persists
const T width = Debye*13.;// domain length
const T length = width/4.;//todo why divided by four?
const T eField = 250.;// electric field in V/m
const T density = 1000.;// carrier fluid density

T psiL1RelError = 0;
T psiL2RelError = 0;
T psiLinfRelError = 0;

T concL1RelError = 0;
T concL2RelError = 0;
T concLinfRelError = 0;

T conc2L1RelError = 0;
T conc2L2RelError = 0;
T conc2LinfRelError = 0;

T velL1RelError = 0;
T velL2RelError = 0;
T velLinfRelError = 0;

// Potential profile for 1D
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
    T psi = util::tanh( valence*psi0*charge/4./Boltzmann/temperature) * util::exp( -distY/Debye );
    output[0] = 4.*Boltzmann*temperature/charge/valence * util::atanh(psi);
    return true;
  };
};

// Cation and anion concentration profile for 1D
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

// Velocity of carrier fluid by electroosmosis under electric field
template <typename T, typename S, typename DESCRIPTOR>
class VelocityProfile1D : public AnalyticalF3D<T, S>
{
private:
  T y0, y1;

public:
  VelocityProfile1D(T y0_, T y1_) : AnalyticalF3D<T, S>(3), y0(y0_), y1(y1_)
  {
     this->getName() = "VelocityProfile1D";
  };

bool operator()(T output[3], const S x[3])
  {
     output[0] = -dielectricC*eField*psi0/viscosity/density*(1 - (util::exp((x[1]-y0)/Debye)+util::exp((2*y1-(x[1]-y0))/Debye))/(1+util::exp((2*y1-y0)/Debye)));
     output[1] = 0.;
     output[2] = 0.;
     return true;
   };
};

//preparing the lattices
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

  // Material=3,4 -->Dirichlet BC
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLattice, superGeometry, 3);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLattice, superGeometry, 4);

  AnalyticalConst3D<T,T> rhoF( psi0 );
  AnalyticalConst3D<T,T> rho0( 0 );
  std::vector<T> velocity( 2,T() );
  AnalyticalConst3D<T,T> uF( velocity );
  sLattice.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({1, 3, 4}),uF);
  sLattice.defineField<descriptors::SOURCE>(superGeometry.getMaterialIndicator({0, 3, 4}),rho0);

  auto bulkIndicator = superGeometry.getMaterialIndicator({3});
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );

  auto bulkIndicator1 = superGeometry.getMaterialIndicator({0, 1, 4});
  sLattice.iniEquilibrium( bulkIndicator1, rho0, uF );
  sLattice.defineRhoU( bulkIndicator1, rho0, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  sLattice.initialize();

  //clout << "Prepare Lattice Poisson... OK" << std::endl;
}

void prepareLatticeNernstPlanck( UnitConverter<T,DESCRIPTOR> const& converter,
                                 SuperLattice<T, DESCRIPTOR>& sLattice,
                                 SuperGeometry<T,3>& superGeometry)
{

  OstreamManager clout( std::cout,"prepareLatticeNernstPlanck" );
  clout << "Prepare Lattice Nernst-Planck..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1,4 -->bulk dynamics
  sLattice.defineDynamics<BulkDynamicsConc>(superGeometry.getMaterialIndicator({1, 4}));
  // Material=3 -->Bounce Back (Zero Gradient)
  sLattice.defineDynamics<BounceBack>(superGeometry, 3);

  // Material=4 -->Dirichlet BC
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLattice, superGeometry, 4);

  AnalyticalConst3D<T,T> rhoF( C0 ); //C0
  AnalyticalConst3D<T,T> rho0( 0 );
  std::vector<T> velocity( 2,T() );
  AnalyticalConst3D<T,T> uF( velocity );
  sLattice.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({1, 3, 4}),uF);

  sLattice.iniEquilibrium( superGeometry, 4, rhoF, uF );
  sLattice.defineRhoU( superGeometry, 4, rhoF, uF );

  auto bulkIndicator = superGeometry.getMaterialIndicator({0, 1, 3});
  sLattice.iniEquilibrium( bulkIndicator, rho0, uF );
  sLattice.defineRhoU( bulkIndicator, rho0, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  sLattice.initialize();

  clout << "Prepare Lattice Nernst-Planck... OK" << std::endl;
}

void prepareLatticeNSE( UnitConverter<T,DESCRIPTORNSE> const& converter,
                        SuperLattice<T, DESCRIPTORNSE>& sLattice,
                        SuperGeometry<T,3>& superGeometry)
{

  OstreamManager clout( std::cout,"prepareLatticeNSE" );
  clout << "Prepare Lattice NSE..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1, 4 -->bulk dynamics
  sLattice.defineDynamics<BulkDynamicsNSE>(superGeometry.getMaterialIndicator({1, 4}));
  // Material=3 -->No Slip
  sLattice.defineDynamics<BounceBack>(superGeometry, 3);


  // Material=4 -->Zero Gradient
  setZeroGradientBoundary<T,DESCRIPTORNSE>(sLattice, superGeometry, 4);

  AnalyticalConst3D<T,T> rho1( 1. );
  std::vector<T> velocity( 2,T() );
  AnalyticalConst3D<T,T> uF( velocity );
  sLattice.defineField<descriptors::FORCE>(superGeometry.getMaterialIndicator({0, 1, 3, 4}),uF);

  auto bulkIndicator = superGeometry.getMaterialIndicator({0, 1, 3, 4});
  sLattice.iniEquilibrium( bulkIndicator, rho1, uF );
  sLattice.defineRhoU( bulkIndicator, rho1, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  sLattice.initialize();

  clout << "Prepare Lattice NSE... OK" << std::endl;
}

void error( SuperLattice<T, DESCRIPTOR>& sLatticePoisson,
            SuperLattice<T, DESCRIPTOR>& sLatticeCation,
            SuperLattice<T, DESCRIPTOR>& sLatticeAnion,
            SuperLattice<T, DESCRIPTORNSE>& sLatticeNSE,
            UnitConverter<T, DESCRIPTOR> const& converter,
            UnitConverter<T,DESCRIPTORNSE> const& converterNSE,
            SuperGeometry<T,3>& superGeometry, bool var)
{
  OstreamManager clout( std::cout,"Errors" );
  if( var ){
  sLatticePoisson.setProcessingContext(ProcessingContext::Evaluation);
  sLatticeCation.setProcessingContext(ProcessingContext::Evaluation);
  sLatticeAnion.setProcessingContext(ProcessingContext::Evaluation);
  sLatticeNSE.setProcessingContext(ProcessingContext::Evaluation);
  int tmp[]= { };
  T result[2] = { };

  T eps = converter.getPhysDeltaX();
  PotentialProfile1D<T,T,DESCRIPTOR> psiSol(0.*eps, width - 0.*eps, converter);
  ConcentrationProfile1D<T,T,DESCRIPTOR> cationSol(valence, psiSol);
  ConcentrationProfile1D<T,T,DESCRIPTOR> anionSol(-valence, psiSol);
  VelocityProfile1D<T,T,DESCRIPTORNSE> velSol(0.*eps, width - 0.*eps);
  SuperLatticeDensity3D<T, DESCRIPTOR> psi( sLatticePoisson );
  SuperLatticeDensity3D<T, DESCRIPTOR> cation( sLatticeCation);
  SuperLatticeDensity3D<T, DESCRIPTOR> anion( sLatticeAnion);
  SuperLatticePhysVelocity3D<T, DESCRIPTORNSE> velNSE( sLatticeNSE, converterNSE);

  auto material = superGeometry.getMaterialIndicator(1);
  SuperRelativeErrorL1Norm3D<T>   errorPsiL1Norm(psi, psiSol, *material);
  SuperRelativeErrorL2Norm3D<T>   errorPsiL2Norm(psi, psiSol, *material);
  SuperRelativeErrorLinfNorm3D<T> errorPsiLinfNorm(psi, psiSol, *material);

  SuperRelativeErrorL1Norm3D<T>   errorConcL1Norm(cation, cationSol, *material);
  SuperRelativeErrorL2Norm3D<T>   errorConcL2Norm(cation, cationSol, *material);
  SuperRelativeErrorLinfNorm3D<T> errorConcLinfNorm(cation, cationSol, *material);

  SuperRelativeErrorL1Norm3D<T>   errorConc2L1Norm(anion, anionSol, *material);
  SuperRelativeErrorL2Norm3D<T>   errorConc2L2Norm(anion, anionSol, *material);
  SuperRelativeErrorLinfNorm3D<T> errorConc2LinfNorm(anion, anionSol, *material);

  SuperRelativeErrorL1Norm3D<T>   errorVelL1Norm(velNSE, velSol, *material);
  SuperRelativeErrorL2Norm3D<T>   errorVelL2Norm(velNSE, velSol, *material);
  SuperRelativeErrorLinfNorm3D<T> errorVelLinfNorm(velNSE, velSol, *material);

  errorPsiL1Norm(result,tmp);
  clout << "Relative Psi-L1-error: " << result[0] << std::endl;
  psiL1RelError = result[0];

  errorPsiL2Norm(result,tmp);
  clout << "Relative Psi-L2-error: " << result[0] << std::endl;
  psiL2RelError = result[0];

  errorPsiLinfNorm(result,tmp);
  clout << "Relative Psi-Linf-error: " << result[0] << std::endl;
  psiLinfRelError = result[0];

  errorConcL1Norm(result,tmp);
  clout << "Relative [Cation]-L1-error: " << result[0] << std::endl;
  concL1RelError = result[0];

  errorConcL2Norm(result,tmp);
  clout << "Relative [Cation]-L2-error: " << result[0] << std::endl;
  concL2RelError = result[0];

  errorConcLinfNorm(result,tmp);
  clout << "Relative [Cation]-Linf-error: " << result[0] << std::endl;
  concLinfRelError = result[0];

  errorConc2L1Norm(result,tmp);
  clout << "Relative [Anion]-L1-error: " << result[0] << std::endl;
  conc2L1RelError = result[0];

  errorConc2L2Norm(result,tmp);
  clout << "Relative [Anion]-L2-error: " << result[0] << std::endl;
  conc2L2RelError = result[0];

  errorConc2LinfNorm(result,tmp);
  clout << "Relative [Anion]-Linf-error: " << result[0] << std::endl;
  conc2LinfRelError = result[0];

  errorVelL1Norm(result,tmp);
  clout << "Relative Velocity-L1-error: " << result[0] << std::endl;
  velL1RelError = result[0];

  errorVelL2Norm(result,tmp);
  clout << "Relative Velocity-L2-error: " << result[0] << std::endl;
  velL2RelError = result[0];

  errorVelLinfNorm(result,tmp);
  clout << "Relative Velocity-Linf-error: " << result[0] << std::endl;
  velLinfRelError = result[0];
  }
}

void getResults( SuperLattice<T, DESCRIPTOR>& sLatticeCation, SuperLattice<T, DESCRIPTOR>& sLatticeAnion,
                 SuperLattice<T, DESCRIPTOR>& sLatticePoisson, SuperLattice<T, DESCRIPTORNSE>& sLatticeNSE,
                 UnitConverter<T, DESCRIPTOR> const& converter,
                 UnitConverter<T, DESCRIPTORNSE> const& converterNSE, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                 bool hasConverged, Gnuplot<T>& gplot)
{
  OstreamManager clout( std::cout,"getResults" );

  const int vtkIter  = 100000;
  const int statIter = 100000;

  SuperVTMwriter3D<T> vtmWriter( "electroosmosis" );
  T eps = converter.getPhysDeltaX();
  PotentialProfile1D<T,T,DESCRIPTOR> psiSol(0.*eps, width - 0.*eps, converter);
  VelocityProfile1D<T,T,DESCRIPTOR> velSol(0.*eps, width - 0.*eps);
  ConcentrationProfile1D<T,T,DESCRIPTOR> cationSol(valence, psiSol);
  ConcentrationProfile1D<T,T,DESCRIPTOR> anionSol(-valence, psiSol);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalPotential(psiSol, sLatticePoisson);
  analyticalPotential.getName() = "analytical potential solution";
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalVelocity(velSol, sLatticePoisson);
  analyticalVelocity.getName() = "analytical velocity";
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalCation(cationSol, sLatticeCation);
  analyticalCation.getName() = "analytical cation concentration";
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalAnion(anionSol, sLatticeAnion);
  analyticalAnion.getName() = "analytical anion concentration";

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLatticeCation );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLatticeCation );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter == 0 || hasConverged ) {

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLatticeCation.getStatistics().print( iT,converter.getPhysTime( iT ) );
    sLatticeNSE.getStatistics().print( iT,converterNSE.getPhysTime( iT ) );
  }

  // Writes the VTK files
  if ( iT%vtkIter == 0 || hasConverged) {
    sLatticeCation.setProcessingContext(ProcessingContext::Evaluation);
    sLatticeAnion.setProcessingContext(ProcessingContext::Evaluation);
    sLatticePoisson.setProcessingContext(ProcessingContext::Evaluation);
    sLatticeNSE.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticeDensity3D<T,DESCRIPTOR> cation( sLatticeCation );
    cation.getName() = "cation concentration";
    SuperLatticeDensity3D<T,DESCRIPTOR> anion( sLatticeAnion);
    anion.getName() = "anion concentration";
    SuperLatticePhysField3D<T,DESCRIPTOR,VELOCITY> velCation( sLatticeCation, converter.getConversionFactorVelocity() );
    velCation.getName() = "cation velocity";
    SuperLatticePhysField3D<T,DESCRIPTOR,VELOCITY> velAnion( sLatticeAnion, converter.getConversionFactorVelocity() );
    velAnion.getName() = "anion velocity";
    SuperLatticeDensity3D<T,DESCRIPTOR> psi( sLatticePoisson);
    psi.getName() = "potential";
    SuperLatticePhysVelocity3D<T,DESCRIPTORNSE> physVel( sLatticeNSE, converterNSE );
    physVel.getName() = "physVelNSE";

    vtmWriter.addFunctor( cation );
    vtmWriter.addFunctor( velCation );
    vtmWriter.addFunctor( anion );
    vtmWriter.addFunctor( velAnion );
    vtmWriter.addFunctor( psi );
    vtmWriter.addFunctor( analyticalPotential );
    vtmWriter.addFunctor( analyticalVelocity );
    vtmWriter.addFunctor( analyticalCation );
    vtmWriter.addFunctor( analyticalAnion );
    vtmWriter.addFunctor( physVel );
    vtmWriter.write( iT );
   }

  if(hasConverged){
    sLatticeCation.setProcessingContext(ProcessingContext::Evaluation);
    sLatticeAnion.setProcessingContext(ProcessingContext::Evaluation);
    sLatticePoisson.setProcessingContext(ProcessingContext::Evaluation);
    sLatticeNSE.setProcessingContext(ProcessingContext::Evaluation);
    gplot.setData (
          T(converter.getResolution()),
          { psiL1RelError, psiL2RelError, psiLinfRelError,
            concL1RelError, concL2RelError, concLinfRelError,
            conc2L1RelError, conc2L2RelError, conc2LinfRelError,
            velL1RelError, velL2RelError, velLinfRelError},
          { "psi L1 Rel Error","psi L2 Rel Error",
            "psi Linf Rel error","conc L1 Rel Error","conc L2 Rel Error",
            "conc Linf Rel error","conc2 L1 Rel Error","conc2 L2 Rel Error",
            "conc2 Linf Rel error","vel L1 Rel Error","vel  L2 Rel Error",
            "vel  Linf Rel error"},
          "top right",
          { 'p','p','p','p','p','p','p','p','p','p','p','p' } );
    psiL1RelError = 0;
    psiL2RelError = 0;
    psiLinfRelError = 0;
    concL1RelError = 0;
    concL2RelError = 0;
    concLinfRelError = 0;
    conc2L1RelError = 0;
    conc2L2RelError = 0;
    conc2LinfRelError = 0;
    velL1RelError = 0;
    velL2RelError = 0;
    velLinfRelError = 0;
  }
}

void simulatePoissonNernstPlanck3D(int N, Gnuplot<T>& gplot)
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

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTORNSE> const converterNSE(
    (int) N,                        // resolution: number of voxels per charPhysL
    (T)   0.7,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   width,       // charPhysLength: reference length of simulation geometry
    (T)   0.5,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   viscosity, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   density                      // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converterNSE.print();
  // Writes the converter log in a file
  converterNSE.write("nse");

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converterNernstPlanck(
    (int) N,                        // resolution: number of voxels per charPhysL
    (T)   0.7, //1./converterNSE.getLatticeRelaxationFrequencyFromDiffusivity<DESCRIPTOR>(diffusion),                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   width,       // charPhysLength: reference length of simulation geometry
    (T)   0.5,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   diffusion, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.                       // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converterNernstPlanck.print();
  // Writes the converter log in a file
  converterNernstPlanck.write("nernstPlanck");

  // === 2nd Step: Prepare Geometry ===
  T L = converterNernstPlanck.getPhysDeltaX();
  Vector<T,3> extend( length,width, 3*L );
  Vector<T,3> origin( 0,0,0 );
  IndicatorCuboid3D<T> cuboid( extend, origin );

#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, converterNernstPlanck.getPhysDeltaX(), singleton::mpi().getSize() );
#else
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, converterNernstPlanck.getPhysDeltaX(), 2 );
#endif

  cuboidDecomposition.print();
  cuboidDecomposition.setPeriodicity({true,false,true});

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,3> superGeometry( cuboidDecomposition, loadBalancer );
  prepareGeometry( converterNernstPlanck, superGeometry );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice<T, DESCRIPTORNSE> sLatticeNSE( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLatticeCation( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLatticePoisson( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLatticeAnion( superGeometry );

  prepareLatticeNernstPlanck( converterNernstPlanck, sLatticeCation, superGeometry);
  prepareLatticeNernstPlanck( converterNernstPlanck, sLatticeAnion, superGeometry);
  prepareLatticePoisson( converterPoisson, sLatticePoisson, superGeometry );
  prepareLatticeNSE( converterNSE, sLatticeNSE, superGeometry );

  T npVelCoeff = charge * valence * diffusion / Boltzmann / temperature / converterNernstPlanck.getConversionFactorVelocity();
  T sourceCoeff = 1./dielectricC * Faraday * valence * converterPoisson.getConversionFactorTime();
  T forceCoeff = eField * Faraday * valence/density * converterNSE.getConversionFactorMass() / converterNSE.getConversionFactorForce();
  clout << "Debye length: " << Debye << std::endl;

  SuperLatticeCoupling coupling(
    NSPNPCoupling<T>{},
    names::Concentration0{}, sLatticeCation,
    names::Temperature{}, sLatticePoisson,
    names::Concentration1{}, sLatticeAnion,
    names::NavierStokes{}, sLatticeNSE);
  coupling.setParameter<NSPNPCoupling<T>::DX>(converterPoisson.getPhysDeltaX());
  coupling.setParameter<NSPNPCoupling<T>::NPVELCOEFF>(npVelCoeff);
  coupling.setParameter<NSPNPCoupling<T>::POISSONCOEFF>(sourceCoeff);
  coupling.setParameter<NSPNPCoupling<T>::FORCECOEFF>(forceCoeff);
  coupling.setParameter<NSPNPCoupling<T>::DTADE>(converterNernstPlanck.getConversionFactorTime());
  coupling.setParameter<NSPNPCoupling<T>::DTNSE>(converterNSE.getConversionFactorTime());
  coupling.setParameter<NSPNPCoupling<T>::OMEGA>(converterPoisson.getLatticeRelaxationFrequency());
  coupling.restrictTo(superGeometry.getMaterialIndicator({1,3}));

  // === 4th Step: Main Loop with Timer ===

  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converterNernstPlanck.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( 500, 1.e-9 );
  timer.start();

  for ( std::size_t iT = 0; iT < converterNernstPlanck.getLatticeTime( maxPhysT ); ++iT ) {

    sLatticePoisson.collideAndStream();

    coupling.execute();

    // === 6th Step: Collide and Stream Execution ===
    sLatticeNSE.collideAndStream();
    sLatticeCation.collideAndStream();
    sLatticeAnion.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    error( sLatticePoisson, sLatticeCation, sLatticeAnion, sLatticeNSE, converterNernstPlanck, converterNSE, superGeometry, converge.hasConverged());
    getResults( sLatticeCation, sLatticeAnion, sLatticePoisson, sLatticeNSE, converterNernstPlanck, converterNSE, iT, superGeometry, timer, converge.hasConverged(), gplot);

    if(converge.hasConverged()) { break; }
    converge.takeValue( sLatticePoisson.getStatistics().getAverageRho(), false );
  }

  timer.stop();
  timer.printSummary();
}

#endif
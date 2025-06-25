/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Shota Ito
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

/* 3D force-driven Poiseuille flow with different non-Newtonian viscosity models
 * Implemented options are Newtonian, PowerLaw, Casson, and Carreau-Yasuda models.
 * The first three models are validated for the same pressure drop and dynamic
 * viscosity used for the unit-conversion. The pressure drop is computed according
 * the analytical solution for the Newtonian case.
 * Carreau-Yasuda model is provided as a popular alternative model but still
 * remains untested yet.
*/

#ifndef NON_NEWTONIAN_POISEUILLE_3D_H
#define NON_NEWTONIAN_POISEUILLE_3D_H

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using T = double;
using DESCRIPTOR = D3Q19<OMEGA,FORCE>;

int N = 31;                                        // resolution of the model
const T length  = 0.003;                           // length of the pie
const T diameter  = 0.001;                         // diameter of the pipe
const T tau = 0.64265;                             // relaxation time
const T physRho = 1060.;                           // physical density
const T maxPhysT = 0.3;                            // max. simulation time in s, SI unit
const T radius  = diameter/2.;                     // radius of the pipe
const Vector<T, 3> center0(0, radius, radius);
const Vector<T, 3> center1(length, radius, radius);
const std::vector<T> origin = { length, radius, radius };
const std::vector<T> axis = { 1, 0, 0 };
const T Re = 100.;
const T dynamicViscosity = 0.0040755;              // dynamic viscosity in Pa*s, SI unit
const T viscosity = dynamicViscosity / physRho;    // kinematic viscosity in m^2/s, SI unit
const T pressureDrop = (16.*util::pow(dynamicViscosity,2.)*Re)/
  (physRho*util::pow(diameter,3.));                // pressure drop in Pa / m, SI unit

// PowerLaw model parameter
const T n = 0.65;
const T consistencyIndex = util::pow(dynamicViscosity,n)*
  util::pow((radius*pressureDrop/8.)*
  util::pow(n/(n+1),n),1-n);                       // consistency index in Pa s^n, SI unit

// Casson model parameter
const T k0 = 0.07;                                 // k0 constant (Pa)^0.5, k0^2 corresponds to yield stress
const T k1 = util::sqrt(dynamicViscosity);         // k1 constant (Pa*s)^0.5, k1^2 corresponds to Casson viscosity

// CarreauYasuda model parameter
const T n_CY = 0.9;                                // CY fluid index
const T a = 1.5;                                   // model constant
const T lambda = 3.313;                            // characteristic time constant (s)
const T mu_zero = dynamicViscosity*1.45;           // zero viscosity (Pa*s)
const T mu_inf = dynamicViscosity/10.;             // infinity viscosity (Pa*s)
const T n_PL = 0.708;                              // PowerLaw index used for comparison

struct Newtonian
{
  using dynamics = ForcedBGKdynamics<T,DESCRIPTOR>;

  static T charVelocity() {
    return  Re * viscosity / diameter;
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePoiseuille3D<T>>(origin, axis, converter.getCharLatticeVelocity(), radius);
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePoiseuille3D<T>>(origin, axis, converter.getCharPhysVelocity(), radius);
  }

  static void setModelParameter(UnitConverter<T,DESCRIPTOR> const& converter, SuperLattice<T,DESCRIPTOR>& sLattice){};
};

struct PowerLaw
{
  using dynamics = PowerLawForcedBGKdynamics<T,DESCRIPTOR>;

  static T charVelocity() {
    return ((n)/(n+1))*(util::pow(pressureDrop/(2.*consistencyIndex),1./n))*util::pow(radius,(n+1)/n);
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePowerLaw3D<T>>(origin, axis, converter.getCharLatticeVelocity(), radius, n);
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePowerLaw3D<T>>(origin, axis, converter.getCharPhysVelocity(), radius, n);
  }

  static void setModelParameter(UnitConverter<T,DESCRIPTOR> const& converter, SuperLattice<T,DESCRIPTOR>& sLattice){
    T conversionConsistency = util::pow(converter.getPhysDeltaT(),n-2)*util::pow(converter.getPhysDeltaX(),2.);
    sLattice.setParameter<powerlaw::M>(consistencyIndex / (conversionConsistency * physRho));
    sLattice.setParameter<powerlaw::N>(n);
    const T nuMin = 2.9686e-3;
    const T nuMax = 3.1667;
    sLattice.setParameter<powerlaw::OMEGA_MIN>(1./(nuMax*descriptors::invCs2<T,DESCRIPTOR>() + 0.5));
    sLattice.setParameter<powerlaw::OMEGA_MAX>(1./(nuMin*descriptors::invCs2<T,DESCRIPTOR>() + 0.5));
  };
};

struct Casson
{
  using dynamics = CassonForcedBGKdynamics<T,DESCRIPTOR>;

  static T charVelocity() {
    return  Re * viscosity / diameter;
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CircleCasson3D<T>>(origin, axis, radius, k1*k1, pressureDrop, k0*k0, 1.0/converter.getConversionFactorVelocity());
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CircleCasson3D<T>>(origin, axis, radius, k1*k1, pressureDrop, k0*k0, 1.0);
  }

  static void setModelParameter(UnitConverter<T,DESCRIPTOR> const& converter, SuperLattice<T,DESCRIPTOR>& sLattice) {
    const T conversionK1 = util::sqrt( converter.getConversionFactorViscosity());
    const T conversionK0 = conversionK1 / util::sqrt(converter.getConversionFactorTime());
    sLattice.setParameter<visco::K_ZERO>(k0 / (util::sqrt(physRho) * conversionK0));
    sLattice.setParameter<visco::K_ONE>(k1 / (util::sqrt(physRho) * conversionK1));
  };
};

struct CarreauYasuda
{
  using dynamics = CarreauYasudaForcedBGKdynamics<T,DESCRIPTOR>;

  static T charVelocity() {
    return  Re * mu_zero / (diameter * physRho);
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePowerLaw3D<T>>(origin, axis, converter.getCharLatticeVelocity(), radius, n_PL);
  }

  static std::shared_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(UnitConverter<T,DESCRIPTOR> const& converter) {
    return std::make_shared<CirclePowerLaw3D<T>>(origin, axis, converter.getCharPhysVelocity(), radius, n_PL);
  }

  static void setModelParameter(UnitConverter<T,DESCRIPTOR> const& converter, SuperLattice<T,DESCRIPTOR>& sLattice) {
    sLattice.setParameter<visco::N>(n_CY);
    sLattice.setParameter<visco::A>(a);
    sLattice.setParameter<visco::LAMBDA>(lambda / converter.getConversionFactorTime());
    const T conversionMU = converter.getConversionFactorViscosity();
    sLattice.setParameter<visco::MU_ZERO>(mu_zero / (physRho * conversionMU));
    sLattice.setParameter<visco::MU_INF>(mu_inf / (physRho * conversionMU));
  }
};

template<typename MODEL>
using DYNAMICS = MODEL::dynamics;

// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T, 3> center0(-converter.getPhysDeltaX() * 0.2, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  center0[0] -= 3.*converter.getPhysDeltaX();
  center1[0] += 3.*converter.getPhysDeltaX();
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, pipe);
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
template<typename MODEL>
void prepareLattice(SuperLattice<T, DESCRIPTOR>& sLattice,
                    UnitConverter<T, DESCRIPTOR>const& converter,
                    SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();
  sLattice.defineDynamics<DYNAMICS<MODEL>>(superGeometry.getMaterialIndicator({1}));

  std::shared_ptr<AnalyticalF3D<T,T>> profileU = MODEL::getLatticeVelocityProfile(converter);

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  center0[0] -= 0.5*converter.getPhysDeltaX();
  center1[0] += 0.5*converter.getPhysDeltaX();
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  setBouzidiBoundary<T, DESCRIPTOR, BouzidiPostProcessor>(sLattice, superGeometry, 2, pipe);

  Vector<T,3> forceTerm = {0,0,0};
  T conversionF_force = converter.getConversionFactorLength() / util::pow(converter.getConversionFactorTime(), 2.);
  forceTerm[0] = pressureDrop / physRho / conversionF_force;

  AnalyticalConst3D<T,T> force (forceTerm);
  sLattice.defineField<FORCE>(superGeometry, 1, force);
  sLattice.defineField<FORCE>(superGeometry, 2, force);

  AnalyticalConst3D<T,T> rho(1.);

  sLattice.defineRhoU(superGeometry, 1, rho, *profileU);
  sLattice.iniEquilibrium(superGeometry, 1, rho, *profileU);
  sLattice.defineRhoU(superGeometry, 2, rho, *profileU);
  sLattice.iniEquilibrium(superGeometry, 2, rho,*profileU);

  sLattice.setParameter<descriptors::OMEGA>(omega);
  MODEL::setModelParameter(converter, sLattice);

  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Compute error norms
template<typename MODEL>
std::vector<T> error( SuperGeometry<T,3>& superGeometry,
            SuperLattice<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter )
{
  OstreamManager clout( std::cout,"error" );
  std::vector<T> errors;
  int tmp[]= { };
  T result[2]= { };

  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  std::shared_ptr<AnalyticalF3D<T,T>> uSol = MODEL::getPhysVelocityProfile(converter);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( sLattice,converter );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, *uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, *uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;
  errors.push_back(result[0]);

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, *uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, *uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;
  errors.push_back(result[0]);

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, *uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, *uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;
  errors.push_back(result[0]);

  return errors;
}

// Output to console and files
template<typename MODEL>
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer )
{
  OstreamManager clout( std::cout,"getResults" );
  const bool lastTimeStep = (iT + 1 == converter.getLatticeTime( maxPhysT ));
  const int statIter = converter.getLatticeTime( maxPhysT / 10. );

  SuperVTMwriter3D<T> vtmWriter( "nonNewtonianPoiseuille3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  std::shared_ptr<AnalyticalF3D<T,T>> uSol = MODEL::getPhysVelocityProfile(converter);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalVelocityLattice(*uSol, sLattice);
  analyticalVelocityLattice.getName() = "analytical solution";
  vtmWriter.addFunctor(analyticalVelocityLattice);

  SuperLatticeField3D<T,DESCRIPTOR,OMEGA> omega(sLattice);
  omega.getName() = "omega";
  vtmWriter.addFunctor(omega);

  const int vtmIter  = converter.getLatticeTime( maxPhysT/100. );

  if ( iT==0 ) {
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );

    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || lastTimeStep ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write( iT );
  }

  // Writes output on the console
  if ( iT%statIter==0 || lastTimeStep ) {
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    error<MODEL>( superGeometry, sLattice, converter );
  }
}

template<typename MODEL>
std::vector<T> simulatePoiseuilleWith()
{
  OstreamManager clout( std::cout,"simulatePoiseuille" );
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    (T)   N,
    (T)   tau,
    (T)   diameter,
    (T)   MODEL::charVelocity(),
    (T)   viscosity,
    (T)   physRho
  );
  converter.print();
  converter.write("nonNewtonianPoiseuille3d");

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

  cuboidDecomposition.setPeriodicity({true, false, false});

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // Instantiation of a superGeometry
  const int overlap = 3;
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, overlap);

  prepareGeometry(converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and setBoundaryConditions
  prepareLattice<MODEL>(sLattice, converter, superGeometry);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  std::vector<T> errors;
  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    sLattice.collideAndStream();
    getResults<MODEL>( sLattice, converter, iT, superGeometry, timer );
  }
  timer.stop();
  timer.printSummary();

  return error<MODEL>(superGeometry, sLattice, converter);
}

template<typename MODEL>
void eocPoiseuilleWith() {
  OstreamManager clout( std::cout,"eocPoiseuille" );
  std::vector<int> res {41, 61, 81, 101};
  std::vector<std::vector<T>> errors;
  for (int iter = 0; iter < (int) res.size(); ++iter) {
    N = res[iter];
    errors.push_back(simulatePoiseuilleWith<MODEL>());
  }

  std::vector<T> eoc_L1, eoc_L2, eoc_Linf;
  for (int i = 0; i < (int) res.size() - 1; ++i) {
    eoc_L1.push_back( (util::log(errors.at(i).at(0)) - util::log(errors.at(i+1).at(0)))/(util::log(res[i]) - util::log(res[i+1])) );
    eoc_L2.push_back( (util::log(errors.at(i).at(1)) - util::log(errors.at(i+1).at(1)))/(util::log(res[i]) - util::log(res[i+1])) );
    eoc_Linf.push_back( (util::log(errors.at(i).at(2)) - util::log(errors.at(i+1).at(2)))/(util::log(res[i]) - util::log(res[i+1])) );
  }

  for (int i = 0; i < (int) res.size() -1; ++i) {
    clout << "Error L1: " << errors.at(i).at(0) << std::endl;
    clout << "Error L2: " << errors.at(i).at(1) << std::endl;
    clout << "Error Linf: " << errors.at(i).at(2) << std::endl;
  }
  clout << "EOC with N: " << res << std::endl;
  clout << "EOC L1: " << eoc_L1 << std::endl;
  clout << "EOC L2: " << eoc_L2 << std::endl;
  clout << "EOC Linf: " << eoc_Linf << std::endl;
}

int main(int argc, char* argv[])
{
  initialize( &argc, &argv );
  CLIreader args(argc, argv);
  N = args.getValueOrFallback("--resolution", 31);
  const bool eoc = args.contains("--eoc");
  const std::string viscosityModel = args.getValueOrFallback<std::string>("--model", "PowerLaw");

  if (!eoc) {
    if (viscosityModel == "Newtonian") {
      simulatePoiseuilleWith<Newtonian>();
    } else if (viscosityModel == "PowerLaw") {
      simulatePoiseuilleWith<PowerLaw>();
    } else if (viscosityModel == "Casson") {
      simulatePoiseuilleWith<Casson>();
    } else if (viscosityModel == "CarreauYasuda") {
      simulatePoiseuilleWith<CarreauYasuda>();
    } else {
      throw std::runtime_error(viscosityModel + " is not a valid viscosity model");
    }
  } else {
    if (viscosityModel == "Newtonian") {
      eocPoiseuilleWith<Newtonian>();
    } else if (viscosityModel == "PowerLaw") {
      eocPoiseuilleWith<PowerLaw>();
    } else if (viscosityModel == "Casson") {
      eocPoiseuilleWith<Casson>();
    } else if (viscosityModel == "CarreauYasuda") {
      eocPoiseuilleWith<CarreauYasuda>();
    } else {
      throw std::runtime_error(viscosityModel + " is not a valid viscosity model");
    }
  }
}

#endif

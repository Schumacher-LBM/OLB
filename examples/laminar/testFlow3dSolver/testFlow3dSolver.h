/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2021 Mathias J. Krause, Fabian Klemens,
 *  Julius Jessberger
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

#ifndef TEST_FLOW_SOLVER_H
#define TEST_FLOW_SOLVER_H

/** \file A fluid flow simulation test case with analytical solution - implementation
 * cf. https://publikationen.bibliothek.kit.edu/1000019768.
 */

#include <olb.h>

#include "analyticalSolutionTestFlow3D.h"

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

using Lattices = meta::map<
  NavierStokes, descriptors::D3Q19<descriptors::FORCE>
>;

enum BoundaryType {LOCAL, INTERPOLATED, BOUZIDI};


// =============== Base class for TestFlow3d implementation ===================
/// define struct that keeps the simulation parameters
template<typename T, typename LATTICES>
struct TfSimulationParams : public parameters::XmlSimulation<T,LATTICES>
{
  // the simulation domain can either be the unit sphere of a unit cube
  bool                   sphericDomain;
  // different LBM boundary conditions can be applied for the (fixed velocity) BC
  BoundaryType           bcType;
};

/// read simulation parameters from xml file
template<typename T, typename LATTICES, typename TAG>
struct parameters::Reader<TfSimulationParams<T,LATTICES>, TAG>
 : public parameters::ReaderBase<TfSimulationParams<T,LATTICES>>
{
  using parameters::ReaderBase<TfSimulationParams<T,LATTICES>>::ReaderBase;

  void read(XMLreader const& xml) {
    parameters::Reader<parameters::XmlSimulation<T,LATTICES>, TAG>(this->params).read(xml);

    std::string helper("cube");
    xml.readOrWarn<std::string>("Application", "Domain", "",
      helper, true, true, true);
    if (helper=="cube") {
      this->params->sphericDomain = false;
    } else if (helper=="sphere") {
      this->params->sphericDomain = true;
    } else {
      throw std::invalid_argument("invalid argument in xml [Application][Domain]");
    }

    std::string helper2("interpolated");
    xml.readOrWarn<std::string>("Application", "BoundaryCondition", "",
      helper2, true, false, true);
    if (helper2=="interpolated") {
      this->params->bcType = INTERPOLATED;
    } else if (helper2=="local") {
      this->params->bcType = LOCAL;
    } else if (helper2=="bouzidi") {
      this->params->bcType = BOUZIDI;
    } else {
      throw std::invalid_argument("invalid argument in xml [Application][BoundaryCondition]");
    }
  }
};

/** struct where the simulation results are written
 * (error norms w.r.t. analytical solution)
 */
template<typename T>
struct SimulationErrors : public parameters::ResultsBase
{
  T velocityAbsL1Error   {0};
  T velocityAbsL2Error   {0};
  T velocityAbsLinfError {0};

  T pressureAbsL1Error   {0};
  T pressureAbsL2Error   {0};
  T pressureAbsLinfError {0};

  T strainRateAbsL1Error   {0};
  T strainRateAbsL2Error   {0};
  T strainRateAbsLinfError {0};

  T dissipationAbsL1Error   {0};
  T dissipationAbsL2Error   {0};
  T dissipationAbsLinfError {0};
};

/** solver implementation: base class
 * Three variants will be implemented via inheritance:
 * - direct simulation
 * - parameter optimization: identify the prefactors of the force field
 *   component-wise (3 parameters) or cuboid-wise for arbitrary number of cuboids
 * - parameter optimization: identify the force field point-wise via adjoint
 *   optimization
 * In any case, we inherit from AdjointLbSolver which makes no difference
 * per default but allows adjoint optimization as well.
 */
template<
  typename T,
  typename PARAMETERS,
  typename LATTICES>
class TestFlowBase : virtual public LbSolver<T,PARAMETERS,LATTICES>
{
private:
  mutable OstreamManager            clout {std::cout, "TestFlowSolver"};

public:
  TestFlowBase(utilities::TypeIndexedSharedPtrTuple<PARAMETERS> params)
   : TestFlowBase::LbSolver(params)
  { }

  using descriptor = typename LATTICES::values_t::template get<0>;

protected:
  virtual void prepareGeometry() override
  {
    // construct sphere (only used for spherical domain)
    std::vector<T> origin(3,T());
    std::vector<T> extend(3, this->converter().getCharPhysLength());
    std::vector<T> center(origin);
    center[0]+=extend[0]/2.;
    center[1]+=extend[1]/2.;
    center[2]+=extend[2]/2.;
    IndicatorSphere3D<T> sphere(center, this->converter().getCharPhysLength()/2.);
    IndicatorLayer3D<T> sphereLayer(sphere, this->converter().getPhysLength(1));

    if (this->parameters(Simulation()).sphericDomain) {
      this->_cGeometry = std::make_shared<CuboidDecomposition3D<T>>(
        sphereLayer,
        this->converter().getPhysLength(1),
        this->parameters(Simulation()).noC * singleton::mpi().getSize() );
    }
    else {
      this->_cGeometry = std::make_shared<CuboidDecomposition3D<T>> (
        origin,
        this->converter().getPhysLength(1),
        std::vector<int>(3, 1/this->converter().getPhysLength(1)+1),
        this->parameters(Simulation()).noC * singleton::mpi().getSize() );
    }

    this->_loadBalancer = std::make_shared<HeuristicLoadBalancer<T>> (
      *this->_cGeometry);

    this->_sGeometry = std::make_shared<SuperGeometry<T,3>> (
      *this->_cGeometry,
      *this->_loadBalancer,
      this->parameters(Simulation()).overlap);

    if (this->parameters(Simulation()).sphericDomain) {
      this->geometry().rename(0,2,sphereLayer);
      this->geometry().rename(2,1,sphere);
    }
    else {
      this->geometry().rename(0,2);
      this->geometry().rename(2,1,{1,1,1});
    }
  }

  virtual void prepareLattices() override
  {
    setDynamics();
    setForceField();
  }

  void setDynamics()
  {
    auto& lattice = this->lattice();
    auto bulk = this->geometry().getMaterialIndicator({1});
    auto boundary = this->geometry().getMaterialIndicator({2});

    // material=0 -->do nothing

    // material=1 -->bulk dynamics
    using BulkDynamics = ForcedBGKdynamics<T,descriptor>;
    lattice.template defineDynamics<BulkDynamics>(bulk);

    // material=2 -->Dirichlet non-zero
    switch (this->parameters(Simulation()).bcType) {
      case LOCAL:
        boundary::set<boundary::LocalVelocity>(
          lattice, boundary);
        break;
      case INTERPOLATED:
        boundary::set<boundary::InterpolatedVelocity>(lattice, boundary);
        break;
      case BOUZIDI:
        std::vector<T> center(3,0.5);
        IndicatorSphere3D<T> sphere(center,
          this->converter().getCharPhysVelocity()/2.);
        setBouzidiBoundary<T,descriptor,BouzidiVelocityPostProcessor>(
          lattice, this->geometry(), 2, sphere);
        break;
    }

    const T omega = this->converter().getLatticeRelaxationFrequency();
    lattice.template setParameter<descriptors::OMEGA>(omega);
  }

  void setForceField()
  {
    // set force field (on whole domain)
    ForceTestFlow3D<T,T,descriptor> forceSol(this->converter());
    const T latticeScaling(this->converter().getConversionFactorMass()
      / this->converter().getConversionFactorForce());
    AnalyticalScaled3D<T,T> force(forceSol, latticeScaling);  // conversion to lattice units
    this->lattice().template defineField<descriptors::FORCE>(
      this->geometry().getMaterialIndicator({1}),
      force);
  }

  virtual void setInitialValues() override
  {
    // initialize lattice population
    AnalyticalConst3D<T,T> rhoF(1);
    std::vector<T> velocity(3,T());
    AnalyticalConst3D<T,T> uF(velocity);
    this->lattice().defineRhoU(
      this->geometry().getMaterialIndicator({1,2}),rhoF, uF);
    this->lattice().iniEquilibrium(
      this->geometry().getMaterialIndicator({1,2}), rhoF, uF);
  }

  virtual void setBoundaryValues(std::size_t iT) override
  {
    // smoothly scale inflow velocity (in order to reduce pressure waves)
    const std::size_t itMaxStart = this->converter().getLatticeTime(
      this->parameters(Simulation()).startUpTime);
    if (iT <= itMaxStart) {
      // compute scaling factor
      PolynomialStartScale<T,T> startScale(itMaxStart, T(1));
      T iTvec[1] = {T(iT)};
      T frac[1] = {};
      startScale(frac, iTvec);

      // analytical solution
      VelocityTestFlow3D<T,T,descriptor> velocity(this->converter());
      AnalyticalScaled3D<T,T> uBoundaryStart(velocity,
        frac[0] / this->converter().getConversionFactorVelocity());
      switch (this->parameters(Simulation()).bcType) {
        case BOUZIDI:
          setBouzidiVelocity(this->lattice(), this->geometry(), 2, uBoundaryStart);
          break;
        default:
          this->lattice().defineU(this->geometry(), 2, uBoundaryStart);
      }
    }
  }

  virtual void writeVTK(std::size_t iT) const override
  {
    SuperVTMwriter3D<T> writer(this->parameters(VisualizationVTK()).filename);

    writeFunctorsToVTK(writer, iT);
  }

  void writeFunctorsToVTK(SuperVTMwriter3D<T>& writer, int index) const
  {
    SuperLatticePhysVelocity3D<T,descriptor> velocity(this->lattice(), this->converter());
    SuperLatticePhysPressure3D<T,descriptor> pressure(this->lattice(), this->converter());
    SuperLatticeField3D<T,descriptor,descriptors::FORCE> force(
      this->lattice());

    std::shared_ptr<AnalyticalF3D<T,T>> solution
     = std::make_shared<VelocityTestFlow3D<T,T,descriptor>>(this->converter());
    SuperLatticeFfromAnalyticalF3D<T,descriptor> solution_lattice(solution, this->lattice());
    solution_lattice.getName() = "analytical_solution";
    writer.addFunctor(solution_lattice);

    SuperLatticeDensity3D<T,descriptor> density(this->lattice());
    writer.addFunctor(density);

    DensityTestFlow3D<T,T,descriptor> analytical_density(this->converter());
    SuperLatticeFfromAnalyticalF3D<T,descriptor> solution_density_lattice(analytical_density, this->lattice());
    solution_density_lattice.getName() = "analytical_density";
    writer.addFunctor(solution_density_lattice);

    writer.addFunctor(force);
    writer.addFunctor(velocity);
    writer.addFunctor(pressure);
    writer.write(index);
  }

  virtual void computeResults() override
  {
    // error computation
    if constexpr (PARAMETERS::keys_t::template contains<Errors>())
    {
      auto& converter = this->converter();
      auto& lattice = this->lattice();
      auto& superGeometry = this->geometry();
      auto& results = this->parameters(Errors());

      OstreamManager clout(std::cout,"error");

      int tmp[1] = {0}; T result[1]; T result1[1];
      VelocityTestFlow3D<T,T,descriptor> uSol(converter);
      SuperLatticePhysVelocity3D<T,descriptor> u(lattice,converter);
      SuperLatticeFfromAnalyticalF3D<T,descriptor> uSolLattice(uSol,lattice);
      PressureTestFlow3D<T,T,descriptor> pSol(converter);
      SuperLatticePhysPressure3D<T,descriptor> p(lattice,converter);
      SuperLatticeFfromAnalyticalF3D<T,descriptor> pSolLattice(pSol,lattice);
      StrainRateTestFlow3D<T,T,descriptor> sSol(converter);
      SuperLatticePhysStrainRate3D<T,descriptor> s(lattice,converter);
      SuperLatticeFfromAnalyticalF3D<T,descriptor> sSolLattice(sSol,lattice);
      DissipationTestFlow3D<T,T,descriptor> dSol(converter);
      SuperLatticePhysDissipation3D<T,descriptor> d(lattice,converter);
      SuperLatticeFfromAnalyticalF3D<T,descriptor> dSolLattice(dSol,lattice);

      // velocity error
      SuperL1Norm3D<T> uL1Norm(uSolLattice-u,superGeometry,1);
      SuperL1Norm3D<T> uSolL1Norm(uSolLattice,superGeometry,1);
      uL1Norm(result,tmp); uSolL1Norm(result1,tmp);
      results.velocityAbsL1Error = result[0];
      clout << "velocity-L1-error(abs)=" << result[0] << "; velocity-L1-error(rel)=" << result[0]/result1[0] << std::endl;

      SuperL2Norm3D<T> uL2Norm(uSolLattice-u,superGeometry,1);
      SuperL2Norm3D<T> uSolL2Norm(uSolLattice,superGeometry,1);
      uL2Norm(result,tmp); uSolL2Norm(result1,tmp);
      results.velocityAbsL2Error = result[0];
      clout << "velocity-L2-error(abs)=" << result[0] << "; velocity-L2-error(rel)=" << result[0]/result1[0] << std::endl;

      SuperLinfNorm3D<T> uLinfNorm(uSolLattice-u,superGeometry,1);
      SuperLinfNorm3D<T> uSolLinfNorm(uSolLattice,superGeometry,1);
      uLinfNorm(result,tmp); uSolLinfNorm(result1,tmp);
      results.velocityAbsLinfError = result[0];
      clout << "velocity-Linf-error(abs)=" << result[0] << "; velocity-Linf-error(rel)=" << result[0]/result1[0] << std::endl;

      // pressure error
      SuperL1Norm3D<T> pL1Norm(pSolLattice-p,superGeometry,1);
      SuperL1Norm3D<T> pSolL1Norm(pSolLattice,superGeometry,1);
      pL1Norm(result,tmp); pSolL1Norm(result1,tmp);
      results.pressureAbsL1Error = result[0];
      clout << "pressure-L1-error(abs)=" << result[0] << "; pressure-L1-error(rel)=" << result[0]/result1[0] << std::endl;

      SuperL2Norm3D<T> pL2Norm(pSolLattice-p,superGeometry,1);
      SuperL2Norm3D<T> pSolL2Norm(pSolLattice,superGeometry,1);
      pL2Norm(result,tmp); pSolL2Norm(result1,tmp);
      results.pressureAbsL2Error = result[0];
      clout << "pressure-L2-error(abs)=" << result[0] << "; pressure-L2-error(rel)=" << result[0]/result1[0] << std::endl;

      SuperLinfNorm3D<T> pLinfNorm(pSolLattice-p,superGeometry,1);
      SuperLinfNorm3D<T> pSolLinfNorm(pSolLattice,superGeometry,1);
      pLinfNorm(result,tmp); pSolLinfNorm(result1,tmp);
      results.pressureAbsLinfError = result[0];
      clout << "pressure-Linf-error(abs)=" << result[0] << "; pressure-Linf-error(rel)=" << result[0]/result1[0] << std::endl;

      // strain rate error
      SuperL1Norm3D<T> sL1Norm(sSolLattice-s,superGeometry,1);
      SuperL1Norm3D<T> sSolL1Norm(sSolLattice,superGeometry,1);
      sL1Norm(result,tmp); sSolL1Norm(result1,tmp);
      results.strainRateAbsL1Error = result[0];
      clout << "strainRate-L1-error(abs)=" << result[0] << "; strainRate-L1-error(rel)=" << result[0]/result1[0] << std::endl;

      SuperL2Norm3D<T> sL2Norm(sSolLattice-s,superGeometry,1);
      SuperL2Norm3D<T> sSolL2Norm(sSolLattice,superGeometry,1);
      sL2Norm(result,tmp); sSolL2Norm(result1,tmp);
      results.strainRateAbsL2Error = result[0];
      clout << "strainRate-L2-error(abs)=" << result[0] << "; strainRate-L2-error(rel)=" << result[0]/result1[0] << std::endl;

      SuperLinfNorm3D<T> sLinfNorm(sSolLattice-s,superGeometry,1);
      sLinfNorm(result,tmp);
      results.strainRateAbsLinfError = result[0];
      clout << "strainRate-Linf-error(abs)=" << result[0] << std::endl;

      // dissipation error
      SuperL1Norm3D<T> dL1Norm(dSolLattice-d,superGeometry,1);
      SuperL1Norm3D<T> dSolL1Norm(dSolLattice,superGeometry,1);
      dL1Norm(result,tmp); dSolL1Norm(result1,tmp);
      results.dissipationAbsL1Error = result[0];
      clout << "dissipation-L1-error(abs)=" << result[0] << "; dissipation-L1-error(rel)=" << result[0]/result1[0] << std::endl;

      SuperL2Norm3D<T> dL2Norm(dSolLattice-d,superGeometry,1);
      SuperL2Norm3D<T> dSolL2Norm(dSolLattice,superGeometry,1);
      dL2Norm(result,tmp); dSolL2Norm(result1,tmp);
      results.dissipationAbsL2Error = result[0];
      clout << "dissipation-L2-error(abs)=" << result[0] << "; dissipation-L2-error(rel)=" << result[0]/result1[0] << std::endl;

      SuperLinfNorm3D<T> dLinfNorm(dSolLattice-d,superGeometry,1);
      dLinfNorm(result,tmp);
      results.dissipationAbsLinfError = result[0];
      clout << "dissipation-Linf-error(abs)=" << result[0] << std::endl;
    }
  }
};


// =================== Variant A: standard simulation =========================

template<typename T>
using Params_TfBasic = meta::map<
  Simulation,       TfSimulationParams<T,Lattices>,
  Output,           parameters::OutputGeneral<T>,
  VisualizationVTK, parameters::OutputPlot<T>,
  Errors,           SimulationErrors<T>
>;

template<typename T>
using TestFlowSolver = TestFlowBase<T,Params_TfBasic<T>,Lattices>;

// Further variants suitable for optimization can be found in
// examples/optimization/testFlowOpti3d.h

#endif

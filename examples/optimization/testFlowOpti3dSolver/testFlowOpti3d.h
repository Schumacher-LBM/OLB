/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2024 Mathias J. Krause, Fabian Klemens,
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

#ifndef TEST_FLOW_OPTI_H
#define TEST_FLOW_OPTI_H


/** \file A fluid flow simulation test case with analytical solution - implementation
 * of parameter identification using an optimization approach,
 * cf. https://publikationen.bibliothek.kit.edu/1000019768.
 */

// There, the standard setup of this simulation is implemented
#include "../../laminar/testFlow3dSolver/testFlow3dSolver.h"

using DualLattices = meta::map<
  NavierStokes, descriptors::DualForcedD3Q19Descriptor
>;

enum TestFlowOptiMode {VELOCITY, DISSIPATION};
enum OptiReferenceSolution {ANALYTICAL, DISCRETE};



// =================== Variant B: optimizable solver ==========================

/// Parameters for optimizable TestFlow: has to inherit interface from OptiSimulation
/** cuboidWise=true: a cuboid decomposition of the domain will be created and the
 * control values are applied as force factors for each cuboid.
 * cuboidWise=false: 3 control values are expected and applied globally as
 * force factors (component-wise).
 */
template<typename T>
struct TfDirectOpti : public parameters::DirectOptiSimulation<T> {

  // cf. above class description
  bool                      cuboidWise        {false};
  /// scales the force field (comes from control variable)
  std::vector<T>            forceFactor       {};
  /// material number, where the objective is computed
  int                       objectiveMaterial {1};
  /// define whether VELOCITY or DISSIPATION error is used for objective computation
  TestFlowOptiMode          optiMode          {VELOCITY};
  /// define whether DISCRETE or ANALYTICAL solution is used to compute the objective
  OptiReferenceSolution     optiReferenceMode {DISCRETE};
  bool                      computeObjective  {true};

  /// here, we can store a discrete reference solution for objective computation
  using BT = BaseType<T>;
  std::shared_ptr<SuperF3D<BT>> referenceState;

  /// define how the control variable enters the simulation
  void applyControl(const std::vector<T>& control) override {
    forceFactor = control;
  }
};

/// XML interface for TestFlowOptiSimulation parameters
template<typename T, typename TAG>
struct parameters::Reader<TfDirectOpti<T>, TAG>
 : public parameters::ReaderBase<TfDirectOpti<T>>
{
  using parameters::ReaderBase<TfDirectOpti<T>>::ReaderBase;

  void read(XMLreader const& xml) {
    xml.readOrWarn<bool>("Optimization", "CuboidWiseControl", "",
      this->params->cuboidWise, true, false, true);
    xml.readOrWarn<int>("Optimization", "ObjectiveMaterial", "",
      this->params->objectiveMaterial, true, false, false);

    std::string helper("velocity");
    xml.readOrWarn<std::string>("Optimization", "TestFlowOptiMode", "",
      helper, true, false, true);
    if (helper=="velocity") {
      this->params->optiMode = VELOCITY;
    } else if (helper=="dissipation") {
      this->params->optiMode = DISSIPATION;
    } else {
      throw std::invalid_argument("invalid argument in xml [Optimization][TestFlowOptiMode]");
    }

    helper = "analytical";
    xml.readOrWarn<std::string>("Optimization", "OptiReferenceMode", "",
      helper, true, false, true);
    if (helper=="discrete") {
      this->params->optiReferenceMode = DISCRETE;
    } else if (helper=="analytical") {
      this->params->optiReferenceMode = ANALYTICAL;
    } else {
      throw std::invalid_argument("invalid argument in xml [Optimization][OptiReferenceMode]");
    }
  }
};

/** Optimization-related results of test-flow simulation
 * Additionally save velocity or dissipation functor to allow usage of
 * discrete objective
 */
template<typename T>
struct TfDirectOptiResults : public parameters::DirectOptiResults<T> {

  std::shared_ptr<SuperF3D<T>> solution;
};


template<typename T>
using Params_TfDirectOpti = meta::map<
  Opti,             TfDirectOpti<T>,
  Output,           parameters::OutputGeneral<T>,
  OutputOpti,       parameters::OptiOutput<T>,
  Results,          TfDirectOptiResults<T>,
  Simulation,       TfSimulationParams<T,Lattices>,
  VisualizationVTK, parameters::OutputPlot<T>
>;


template<typename T>
class TestFlowSolverDirectOpti : public TestFlowBase<T,Params_TfDirectOpti<T>,Lattices> {
private:
  mutable OstreamManager            clout {std::cout, "TestFlowSolverDirectOpti"};
  using descriptor = typename Lattices::values_t::template get<0>;

public:
  TestFlowSolverDirectOpti(
    utilities::TypeIndexedSharedPtrTuple<Params_TfDirectOpti<T>> params)
   : TestFlowSolverDirectOpti::LbSolver(params),
    TestFlowSolverDirectOpti::TestFlowBase(params)
  { }

protected:
  void prepareLattices() override
  {
    // set dynamics as in standard simulation
    this->setDynamics();

    // scale force field according to control
    std::shared_ptr<AnalyticalF<3,T,T>> controlF;
    if (this->parameters(Opti()).cuboidWise) {
      controlF = std::make_shared<AnalyticalCuboidwiseConst<3,T,T>>(
        this->geometry(), this->parameters(Opti()).forceFactor);
    } else {
      OLB_ASSERT((this->parameters(Opti()).forceFactor.size()==3),
        "A control vector of size 3 is required");
      controlF = std::make_shared<AnalyticalConst3D<T,T>>(
        this->parameters(Opti()).forceFactor);
    }
    OLB_ASSERT((controlF->getTargetDim()==3),
      "Dimension of force scaling functor must be 3");

    ForceTestFlow3D<T,T,descriptor> forceSol(this->converter());
    const T latticeScaling(this->converter().getConversionFactorMass()
      / this->converter().getConversionFactorForce());
    std::shared_ptr<AnalyticalF<3,T,T>> force
     = std::make_shared<AnalyticalScaled3D<T,T>>(forceSol, latticeScaling);
    force = force * controlF;

    this->lattice().template defineField<descriptors::FORCE>
    (this->geometry().getMaterialIndicator({1, 2}), *force);
  }

  void computeResults() override
  {
    // objective computation
    if (this->parameters(Opti()).computeObjective) {
      this->parameters(Results()).objective = computeObjective();
    }

    // save simulated solution as a functor
    switch (this->parameters(Opti()).optiMode)
    {
      case DISSIPATION:
      {
        this->parameters(Results()).solution
        = std::make_shared<SuperLatticePhysDissipation3D<T,descriptor>>(
          this->lattice(), this->converter());
      } break;

      case VELOCITY:
      {
        this->parameters(Results()).solution
        = std::make_shared<SuperLatticePhysVelocity3D<T,descriptor>>(
          this->lattice(), this->converter());
      } break;
    }
  }

  /** compute velocity or dissipation L2-error
   * objective = 0.5 * error^2
   */
  T computeObjective() const
  {
    int tmp[1];
    T output[1];

    using BT = BaseType<T>;
    std::shared_ptr<AnalyticalF<3,BT,BT>> help;
    std::shared_ptr<AnalyticalF<3,T,T>> analytical;
    std::shared_ptr<SuperF3D<T,T>>      simulated;

    switch (this->parameters(Opti()).optiMode)
    {
      case DISSIPATION:
      {
        // assign analytical and simulated dissipation for error computation
        switch (this->parameters(Opti()).optiReferenceMode)
        {
          case ANALYTICAL:
          {
            analytical = std::make_shared<DissipationTestFlow3D<T,T,descriptor>>(
              this->converter());
          } break;

          case DISCRETE:
          {
            // user has to guarantee that reference solution is computed
            // cf. node Optimization -> ReferenceSolution in xml file
            help = std::make_shared<AnalyticalFfromSuperF3D<BT>>(
              *(this->parameters(Opti()).referenceState));
            // cast reference solution to T-valued functor
            // for opti with difference quotients, this does nothing
            // for opti with automatic differentiation, we cast double to ADf
            analytical = std::make_shared<AnalyticalTypecast<3,BT,T,BT,T>>(*help);
          } break;
        }

        simulated = std::make_shared<SuperLatticePhysDissipation3D<T,descriptor>>(
          this->lattice(), this->converter());
      } break;

      case VELOCITY:
      {
        // assign analytical and simulated velocity field for error computation
        switch (this->parameters(Opti()).optiReferenceMode)
        {
          case ANALYTICAL:
          {
            analytical = std::make_shared<VelocityTestFlow3D<T,T,descriptor>>(
              this->converter());
          } break;

          case DISCRETE:
          {
            // user has to guarantee that reference solution is computed
            // cf. node Optimization -> ReferenceSolution in xml file
            help = std::make_shared<AnalyticalFfromSuperF3D<BT>>(
              *(this->parameters(Opti()).referenceState));
            analytical = std::make_shared<AnalyticalTypecast<3,BT,T,BT,T>>(*help);
          } break;
        }

        simulated = std::make_shared<SuperLatticePhysVelocity3D<T,descriptor>>(
          this->lattice(), this->converter());
      } break;
    }

    RelativeDifferenceObjective3D<T,descriptor>(
      this->lattice(),
      simulated,
      analytical,
      this->geometry().getMaterialIndicator(this->parameters(Opti()).objectiveMaterial)
    ).operator()(output, tmp);

    return output[0];
  }
};


// ================ Variant C: Adjoint optimizable solver =====================

/// Parameters for optimizable TestFlow: has to inherit interface from OptiSimulation
template<typename T, opti::SolverMode MODE>
struct TestFlowAdjointOpti
 : public parameters::DistributedOptiSimulation<T,DualLattices,MODE> {

  int                       objectiveMaterial {1};
  TestFlowOptiMode          optiMode          {VELOCITY};
  OptiReferenceSolution     optiReferenceMode {DISCRETE};
};

/// XML interface for TestFlowOptiSimulation parameters
template<typename T, opti::SolverMode MODE, typename TAG>
struct parameters::Reader<TestFlowAdjointOpti<T,MODE>, TAG>
 : public parameters::ReaderBase<TestFlowAdjointOpti<T,MODE>>
{
  using parameters::ReaderBase<TestFlowAdjointOpti<T,MODE>>::ReaderBase;

  void read(XMLreader const& xml) {
    Reader<DistributedOptiSimulation<T,DualLattices,MODE>, TAG>(this->params).read(xml);

    xml.readOrWarn<int>("Optimization", "ObjectiveMaterial", "",
      this->params->objectiveMaterial, true, false, false);

    std::string helper("velocity");
    xml.readOrWarn<std::string>("Optimization", "TestFlowOptiMode", "",
      helper, true, false, true);
    if (helper=="velocity") {
      this->params->optiMode = VELOCITY;
    } else if (helper=="dissipation") {
      this->params->optiMode = DISSIPATION;
    } else {
      throw std::invalid_argument("invalid argument in xml [Optimization][TestFlowOptiMode]");
    }

    helper = "analytical";
    xml.readOrWarn<std::string>("Optimization", "OptiReferenceMode", "",
      helper, true, false, true);
    if (helper=="discrete") {
      this->params->optiReferenceMode = DISCRETE;
    } else if (helper=="analytical") {
      this->params->optiReferenceMode = ANALYTICAL;
    } else {
      throw std::invalid_argument("invalid argument in xml [Optimization][OptiReferenceMode]");
    }
  }
};


template<typename T, SolverMode MODE>
using Params_TfOptiAdjoint = meta::map<
  Opti,             TestFlowAdjointOpti<T,MODE>,
  OutputOpti,       parameters::OptiOutput<T>,
  Output,           parameters::OutputGeneral<T>,
  Results,          parameters::DistributedOptiSimulationResults<T,DualLattices>,
  Simulation,       TfSimulationParams<T,DualLattices>,
  VisualizationVTK, parameters::OutputPlot<T>
>;

template<typename T, SolverMode MODE>
class TestFlowSolverOptiAdjoint
 : virtual public TestFlowBase<T,Params_TfOptiAdjoint<T,MODE>,DualLattices>,
   virtual public AdjointLbSolver<T,Params_TfOptiAdjoint<T,MODE>,DualLattices,MODE>
{
private:
  mutable OstreamManager            clout {std::cout, "TestFlowSolverOptiAdjoint"};

public:
  using descriptor = typename TestFlowSolverOptiAdjoint::AdjointLbSolver::DESCRIPTOR;

  TestFlowSolverOptiAdjoint(
    utilities::TypeIndexedSharedPtrTuple<Params_TfOptiAdjoint<T,MODE>> params)
   : TestFlowSolverOptiAdjoint::LbSolver(params),
     TestFlowSolverOptiAdjoint::TestFlowBase(params),
     TestFlowSolverOptiAdjoint::AdjointLbSolver(params)
  {
    const std::string name = (MODE == SolverMode::Reference) ? "Reference_Solution"
      : ((MODE == SolverMode::Primal) ? "Flow_Simulations"
      : "Dual_Simulations");
    this->parameters(VisualizationVTK()).filename = name;
  }

  void prepareGeometry() override
  {
    TestFlowSolverOptiAdjoint::TestFlowBase::prepareGeometry();

    // store some data for optimization
    // this needs to be done in every adjoint optimization
    this->parameters(Opti()).bulkIndicator = std::move(this->geometry().getMaterialIndicator(1));
    this->parameters(Opti()).designDomain = this->parameters(Opti()).bulkIndicator;
    this->parameters(Opti()).objectiveDomain = this->parameters(Opti()).bulkIndicator;
    this->parameters(Results()).geometry = this->_sGeometry;
  }

  void prepareLattices() override
  {
    auto& lattice = this->lattice();

    if constexpr (MODE != SolverMode::Dual) {
      // dynamics from standard simulation
      this->setDynamics();
    }
    else {
      const T omega = this->converter().getLatticeRelaxationFrequency();

      // material=1 -->bulk dynamics
      auto bulk = this->geometry().getMaterialIndicator({1});
      lattice.template defineDynamics<DualForcedBGKDynamics<T,descriptor>>(bulk);

      // material=2 -->Dirichlet zero (bounce back)
      auto boundary = this->geometry().getMaterialIndicator({2});
      boundary::set<boundary::BounceBack>(lattice, boundary);

      lattice.template setParameter<descriptors::OMEGA>(omega);
    }

    if constexpr (MODE == SolverMode::Reference) {
      // standard execution: use analytical force function
      this->setForceField();
    }
    else {
      // optimization: get force from control variables
      lattice.template defineField<descriptors::FORCE>(
        this->geometry().getMaterialIndicator({1}),
        *(this->parameters(Opti()).controlledField));
    }
  }

  void setInitialValues() override
  {
    TestFlowSolverOptiAdjoint::TestFlowBase::setInitialValues();

    // load primal data
    // this needs to be done in every adjoint optimization
    if constexpr ( MODE == SolverMode::Dual ) {
      this->loadPrimalPopulations();
    }
  }

  void setBoundaryValues(std::size_t iT) override
  {
    if constexpr (MODE != SolverMode::Dual) {
      TestFlowSolverOptiAdjoint::TestFlowBase::setBoundaryValues(iT);
    }
    // dual simulation has homogeneous Dirichlet BC, so nothing to do in that case
  }

  void prepareVTK() const override {
    const bool prepare
     = ((MODE == SolverMode::Reference) || (this->parameters(OutputOpti()).counterOptiStep == 0));
    if (prepare) {
      TestFlowSolverOptiAdjoint::LbSolver::prepareVTK();
    }
  }

  void writeVTK(std::size_t iT) const override
  {
    // reference simulation: write solution (unsteady)
    // optimization: write stationary flow field
    SuperVTMwriter3D<T> writer(this->parameters(VisualizationVTK()).filename);
    const int index = (MODE == SolverMode::Reference) ?
      iT : this->parameters(OutputOpti()).counterOptiStep;

    this->writeFunctorsToVTK(writer, index);
  }

  void computeResults() override
  {
    if constexpr (MODE == SolverMode::Reference) {
      // error computation
      TestFlowSolverOptiAdjoint::TestFlowBase::computeResults();
    }
    // store some data for optimization
    TestFlowSolverOptiAdjoint::AdjointLbSolver::computeResults();
  }
};


/// @brief Objective and derivative computation for adjoint optimization
/// Objective = 0.5 * norm(u-uReference)^2
/// @tparam T floating point type
template <typename T>
class TestFlowObjective : public DistributedObjective<T,TestFlowSolverOptiAdjoint>
{
public:
  using descriptor = typename TestFlowSolverOptiAdjoint<T,SolverMode::Reference>::AdjointLbSolver::DESCRIPTOR;
  std::shared_ptr<TestFlowSolverOptiAdjoint<T,SolverMode::Reference>>  _referenceSolver;
  std::shared_ptr<SuperF3D<T,T>>                    _referenceState;

  TestFlowObjective(XMLreader const& xml) {
    // define reference solution
    _referenceSolver = createLbSolver <TestFlowSolverOptiAdjoint<T,SolverMode::Reference>> (xml);

    switch (_referenceSolver->parameters(Opti()).optiMode)
    {
      case DISSIPATION:
      {
        switch (_referenceSolver->parameters(Opti()).optiReferenceMode)
        {
          case ANALYTICAL:
          {
            _referenceSolver->buildAndReturn();  // need to construct the SuperLattice
            std::shared_ptr<AnalyticalF3D<T,T>> analytical
             = std::make_shared<DissipationTestFlow3D<T,T,descriptor>>(
              _referenceSolver->converter());
            _referenceState = std::make_shared<SuperLatticeFfromAnalyticalF3D<T,descriptor>>(
              analytical,
              _referenceSolver->lattice());
          } break;

          case DISCRETE:
          {
            _referenceSolver->solve();
            _referenceState
             = std::make_shared<SuperLatticePhysDissipation3D<T,descriptor>>(
              _referenceSolver->lattice(), _referenceSolver->converter());
          } break;
        }
      } break;

      case VELOCITY:
      {
        switch (_referenceSolver->parameters(Opti()).optiReferenceMode)
        {
          case ANALYTICAL:
          {
            _referenceSolver->buildAndReturn();
            std::shared_ptr<AnalyticalF3D<T,T>> analytical
             = std::make_shared<VelocityTestFlow3D<T,T,descriptor>>(
              _referenceSolver->converter());
            _referenceState = std::make_shared<SuperLatticeFfromAnalyticalF3D<T,descriptor>>(
              analytical,
              _referenceSolver->lattice());
          } break;

          case DISCRETE:
          {
            _referenceSolver->solve();
            _referenceState
             = std::make_shared<SuperLatticePhysVelocity3D<T,descriptor>>(
              _referenceSolver->lattice(), _referenceSolver->converter());
          } break;
        }
      } break;
    }

    // define objective domain
    this->_objectiveDomain = _referenceSolver->parameters(Opti()).designDomain;
    this->_designDomain = _referenceSolver->parameters(Opti()).designDomain;
  }

  T evaluate() override {
    // define functor for simulated solution
    std::shared_ptr<SuperF3D<T,T>>      simulatedState;
    switch (_referenceSolver->parameters(Opti()).optiMode)
    {
      case DISSIPATION:
      {
        simulatedState = std::make_shared<SuperLatticePhysDissipation3D<T,descriptor>>(
          this->_primalSolver->lattice(), _referenceSolver->converter());
      } break;

      case VELOCITY:
      {
        simulatedState = std::make_shared<SuperLatticePhysVelocity3D<T,descriptor>>(
          this->_primalSolver->lattice(), _referenceSolver->converter());
      } break;
    }

    // compare simulated with reference solution
    T output[1]; int tmp[4];
    RelativeDifferenceObjective3D<T,descriptor>(
      this->_primalSolver->lattice(),
      simulatedState,
      _referenceState,
      this->_objectiveDomain
    ).operator()(output, tmp);

    return output[0];
  }

  std::shared_ptr<SuperF3D<T,T>> derivativeByPopulations() override {
    // define functor for simulated state and its derivative by populations
    std::shared_ptr<SuperF3D<T,T>>      simulatedState;
    std::shared_ptr<SuperF3D<T,T>>      simulatedDerivative;

    switch (_referenceSolver->parameters(Opti()).optiMode)
    {
      case DISSIPATION:
      {
        simulatedState = std::make_shared<SuperLatticePhysDissipation3D<T,descriptor>>(
          this->_primalSolver->lattice(),
          this->_primalSolver->converter());
        simulatedDerivative = std::make_shared<SuperLatticeDphysDissipationDf<T,descriptor>> (
          this->_primalSolver->lattice(),
          this->_primalSolver->converter());
      } break;

      case VELOCITY:
      {
        simulatedState = std::make_shared<SuperLatticePhysVelocity3D<T,descriptor>>(
          this->_primalSolver->lattice(),
          this->_primalSolver->converter());
        simulatedDerivative = std::make_shared<SuperLatticeDphysVelocityDf3D<T,descriptor>> (
          this->_primalSolver->lattice(),
          this->_primalSolver->converter());
      } break;
    }

    // define derivative of objective by populations
    return std::make_shared<DrelativeDifferenceObjectiveDf3D<T,descriptor>>(
      this->_primalSolver->lattice(),
      simulatedState,
      simulatedDerivative,
      _referenceState,
      this->_objectiveDomain);
  }

  std::shared_ptr<SuperF3D<T,T>> derivativeByControl() override {
    // objective does not depend (directly) on control
    const Vector<T,3> v(0);
    return std::make_shared<SuperConst3D<T,T>>(this->_primalSolver->lattice(), v);
  }
};

template <typename T>
class RelativeDifferenceVelocityObjectiveGeneric
{
public:
  using Descriptor = descriptors::DualForcedD3Q19Descriptor;
  std::shared_ptr<TestFlowSolverOptiAdjoint<T,SolverMode::Reference>> _referenceSolver;
  std::shared_ptr<SuperF3D<T,T>>                           _referenceState;
protected:
  T _uRefNorm;
  T _conversionVelocity;
  T _regAlpha {0};

public:
  // these indicators are used by the GenericObjective class
  std::shared_ptr<SuperIndicatorF<T,3>>           _objectiveDomain;
  std::shared_ptr<SuperIndicatorF<T,3>>           _designDomain;

  RelativeDifferenceVelocityObjectiveGeneric(XMLreader const& xml) {
    // create reference solver, compute reference solution
    _referenceSolver = createLbSolver <TestFlowSolverOptiAdjoint<T,SolverMode::Reference>> (xml);

    switch (_referenceSolver->parameters(Opti()).optiMode)
    {
      case DISSIPATION:
      {
        switch (_referenceSolver->parameters(Opti()).optiReferenceMode)
        {
          case ANALYTICAL:
          {
            _referenceSolver->buildAndReturn();  // need to construct the SuperLattice
            std::shared_ptr<AnalyticalF3D<T,T>> analytical
             = std::make_shared<DissipationTestFlow3D<T,T,Descriptor>>(
              _referenceSolver->converter());
            _referenceState = std::make_shared<SuperLatticeFfromAnalyticalF3D<T,Descriptor>>(
              analytical,
              _referenceSolver->lattice());
          } break;

          case DISCRETE:
          {
            _referenceSolver->solve();
            _referenceState
             = std::make_shared<SuperLatticePhysDissipation3D<T,Descriptor>>(
              _referenceSolver->lattice(), _referenceSolver->converter());
          } break;
        }
      } break;

      case VELOCITY:
      {
        switch (_referenceSolver->parameters(Opti()).optiReferenceMode)
        {
          case ANALYTICAL:
          {
            _referenceSolver->buildAndReturn();
            std::shared_ptr<AnalyticalF3D<T,T>> analytical
             = std::make_shared<VelocityTestFlow3D<T,T,Descriptor>>(
              _referenceSolver->converter());
            _referenceState = std::make_shared<SuperLatticeFfromAnalyticalF3D<T,Descriptor>>(
              analytical,
              _referenceSolver->lattice());
          } break;

          case DISCRETE:
          {
            _referenceSolver->solve();
            _referenceState
             = std::make_shared<SuperLatticePhysVelocity3D<T,Descriptor>>(
              _referenceSolver->lattice(), _referenceSolver->converter());
          } break;
        }
      } break;
    }

    // define objective domain
    _objectiveDomain = _referenceSolver->parameters(Opti()).designDomain;
    _designDomain = _referenceSolver->parameters(Opti()).designDomain;

    // define constants
    const int dummy[1]{0};
    SuperL2Norm3D<T>(_referenceState, _objectiveDomain)(&_uRefNorm, dummy);
    _conversionVelocity = _referenceSolver->converter().getConversionFactorVelocity();
    // todo dissipation
    xml.readOrWarn<T>("Optimization", "RegAlpha", "", _regAlpha);
  }

  /// Implements the population-dependent term of objective functional
  // This interface is expected by the GenericObjective class. It uses different
  // cell/data types CELL/ V. Hence, beware of typecasts from V to T!
  template <typename CELL, typename V>
  void j(V* output, CELL& cell, int iC, const int* coords) {
    // application-specific implementation
    V diff(0);
    const int iCglob = _referenceSolver->lattice().getLoadBalancer().glob(iC);
    switch (_referenceSolver->parameters(Opti()).optiMode)
    {
      case DISSIPATION:
      {
        const V diss = dissipation(cell);
        T dissRef;
        (*_referenceState)(&dissRef, iCglob, coords[0], coords[1], coords[2]);
        diff = (diss-dissRef) * (diss-dissRef);
      } break;

      case VELOCITY:
      {
        Vector<V,3> u;
        cell.computeU(u.data());
        T uRef[3];
        (*_referenceState)(uRef, iCglob, coords[0], coords[1], coords[2]);
        diff = norm_squared(_conversionVelocity * u - Vector<V,3>(uRef));
      } break;
    }
    // scale and return
    output[0] = diff / (T(2)*_uRefNorm*_uRefNorm);
  }

  /// Implements the control-dependent term of objective functional
  // This interface is expected by GenericObjective, with different data types V.
  // Hence, do not typecasts from V to T!
  template <typename V>
  void r(V* output, int iC, const int* coords, const V* control) {
    // application-specific implementation: compute regularization
    output[0] = V(0);  // T(0.5) * _regAlpha * util::norm2(control, 3);
  }

private:
  template <typename CELL, typename V = typename CELL::value_t>
  V dissipation(CELL& cell) {
    V rho, uTemp[Descriptor::d], pi[util::TensorVal<Descriptor >::n];
    cell.computeAllMomenta(rho, uTemp, pi);

    V PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
    if (util::TensorVal<Descriptor >::n == 6) {
      PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                      + pi[5] * pi[5];
    }

    const auto& converter = _referenceSolver->converter();
    const V dt = converter.getConversionFactorTime();
    const V omega = 1. / converter.getLatticeRelaxationTime();
    const V nuLattice = converter.getLatticeViscosity();

    return PiNeqNormSqr * nuLattice
            * util::pow(omega * descriptors::invCs2<T,Descriptor>() / rho, 2) / 2.
            * converter.getPhysViscosity() / nuLattice / dt / dt;
  }
};

#endif

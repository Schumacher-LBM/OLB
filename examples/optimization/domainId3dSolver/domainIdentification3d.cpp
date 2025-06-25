/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020-24 Fabian Klemens, Mathias J. Krause, Benjamin
 *  Förster, Julius Jeßberger
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

/** \file
 * In this example, a domain identification problem is solved with the help of
 * adjoint LBM. Cf. https://doi.org/10.1016/j.camwa.2018.07.010 for details.
 */

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;
using namespace olb::opti;

using Descriptor = D3Q19<POROSITY,opti::F,opti::DJDF>;
using Lattices = meta::map<
                 NavierStokes, Descriptor
                 >;

template <typename T>
class RelativeDifferenceVelocityObjective;

/// Collect parameters which are relevant for optimization.
// Inherit/ extend DistributedOptiSimulation which has the adjoint-optimization-
// specific data.
template <typename T, SolverMode MODE>
struct TestDomainOptiParameters
 : public parameters::DistributedOptiSimulation<T,Lattices,MODE>
{
  std::shared_ptr<IndicatorF3D<T>>      designDomainAnalytical;
  std::shared_ptr<IndicatorF3D<T>>      objectiveDomainAnalytical;
};

/// Enable reading the parameters in DistributedOptiSimulation from xml
template<typename T, SolverMode MODE, typename TAG>
struct parameters::Reader<TestDomainOptiParameters<T,MODE>, TAG>
 : public parameters::ReaderBase<TestDomainOptiParameters<T,MODE>>
{
  using ReaderBase<TestDomainOptiParameters<T,MODE>>::ReaderBase;

  void read(XMLreader const& xml)
  {
    using namespace parameters;
    Reader<DistributedOptiSimulation<T,Lattices,MODE>, TAG>(this->params).read(xml);
    this->params->objectiveDomainAnalytical = createIndicatorF3D<T> (xml["ObjectiveDomain"], false);
    this->params->designDomainAnalytical = createIndicatorF3D<T> (xml["DesignDomain"]);
  }
};

/// Parameters which are relevant for simulation
template <typename T>
struct TestDomainSimulationParameters
 : public parameters::XmlSimulation<T,Lattices>
{
  std::shared_ptr<IndicatorF3D<T>> simulationObject;
};

/// XML interface for TestDomainSimulationParameters
template<typename T, typename TAG>
struct parameters::Reader<TestDomainSimulationParameters<T>, TAG>
 : public parameters::ReaderBase<TestDomainSimulationParameters<T>>
{
  using ReaderBase<TestDomainSimulationParameters<T>>::ReaderBase;

  void read(XMLreader const& xml) {
    using namespace parameters;
    Reader<XmlSimulation<T,Lattices>, TAG>(this->params).read(xml);
    this->params->simulationObject = createIndicatorF3D<T> (xml["SimulationObject"]);
  }
};



/// Pack all parameter structs together in a map
template<typename T, opti::SolverMode MODE>
using Parameters = meta::map<
  Opti,             TestDomainOptiParameters<T,MODE>,
  Simulation,       TestDomainSimulationParameters<T>,
  Stationarity,     parameters::Stationarity<T>,
  OutputOpti,       parameters::OptiOutput<T>,
  Results,          parameters::DistributedOptiSimulationResults<T,Lattices>,
  Output,           parameters::OutputGeneral<T>,
  VisualizationVTK, parameters::OutputPlot<T>
>;

/// Solver class implements the simulation
// Inherit from AdjointLbSolver since most of the adjoint-specific details
// are implemented there
template<typename T, SolverMode MODE>
class DomainIdSolver : public AdjointLbSolver<T,Parameters<T,MODE>,Lattices,MODE>
{
private:
  mutable OstreamManager                clout {std::cout, "DomainIdSolver"};

public:
  using LATTICES = Lattices;
  using descriptor = typename Lattices::values_t::template get<0>;

  DomainIdSolver(utilities::TypeIndexedSharedPtrTuple<Parameters<T,MODE>> params)
   : DomainIdSolver::LbSolver(params),  // has to be called "by hand" because AdjointLbSolver inherits virtually
    DomainIdSolver::AdjointLbSolver(params)
  {
    const std::string name = (MODE == SolverMode::Reference) ? "Reference_Solution"
      : ((MODE == SolverMode::Primal) ? "Flow_Simulations"
      : "Dual_Simulations");
    this->parameters(VisualizationVTK()).filename = name;
  }

protected:
  /// Set up a cubic geometry and mark the design domain inside
  void prepareGeometry() override
  {
    auto& xml = *(this->parameters(Simulation()).xml);
    // Read mesh geometry from XML file
    std::shared_ptr<IndicatorF3D<T>> cube
                                  = createIndicatorF3D<T>(xml["Application"]["Mesh"], false);
    // Create CuboidDecomposition3D from mesh indicator
    this->_cGeometry = std::make_shared<CuboidDecomposition3D<T>>(
                         *cube,
                         this->converter().getPhysDeltaX(),
                         singleton::mpi().getSize() * this->parameters(Simulation()).noC);
    // Create HeuristicLoadBalancer from CuboidDecomposition3D
    this->_loadBalancer = std::make_shared<HeuristicLoadBalancer<T>>(*this->_cGeometry);
    // Create SuperGeometry
    this->_sGeometry = std::make_shared<SuperGeometry<T,3>> (
                         *this->_cGeometry,
                         *this->_loadBalancer,
                         this->parameters(Simulation()).overlap);

    // Material 2 everywhere
    this->geometry().rename(0,2);
    // Material 1 inside of border
    this->geometry().rename(2,1, {1,1,1});

    // Outer boundaries
    Vector<T,3> origin = {0.,0.,0.};
    origin[0] += this->converter().getPhysDeltaX()/2.;
    origin[1] -= this->converter().getPhysDeltaX()/2.;
    origin[2] += this->converter().getPhysDeltaX()/2.;

    Vector<T,3> extend = {1.,0.,1.};
    extend[0] -= 2.*this->converter().getPhysDeltaX()/2.;
    extend[1] +=    this->converter().getPhysDeltaX()/2.;
    extend[2] -= 2.*this->converter().getPhysDeltaX()/2.;

    IndicatorCuboid3D<T> inflow( extend,origin );
    this->geometry().rename( 2,3,inflow );

    origin = {0.,1.,0.};
    origin[0] += this->converter().getPhysDeltaX()/2.;
    origin[1] -= this->converter().getPhysDeltaX()/2.;
    origin[2] += this->converter().getPhysDeltaX()/2.;

    extend = {1.,1.,1.};
    extend[0] -= 2.*this->converter().getPhysDeltaX()/2.;
    extend[1] +=    this->converter().getPhysDeltaX()/2.;
    extend[2] -= 2.*this->converter().getPhysDeltaX()/2.;

    IndicatorCuboid3D<T> outflow(extend, origin);
    this->geometry().rename( 2,4,outflow );

    // Indicators for material 6 (designDomain)
    this->geometry().rename(1, 6, this->parameters(Opti()).designDomainAnalytical);

    // store some data for optimization
    // this needs to be done in every adjoint optimization
    this->parameters(Opti()).designDomain
     = std::make_shared<SuperIndicatorFfromIndicatorF3D<T>>(
      this->parameters(Opti()).designDomainAnalytical, this->geometry());
    this->parameters(Opti()).bulkIndicator = std::move(this->geometry().getMaterialIndicator({1,6}));
    this->parameters(Opti()).objectiveDomain
     = std::make_shared<SuperIndicatorFfromIndicatorF3D<T>>(
      this->parameters(Opti()).objectiveDomainAnalytical, this->geometry());
    this->parameters(Results()).geometry = this->_sGeometry;

    if constexpr (MODE == SolverMode::Reference) {
      const unsigned nSimulation = this->geometry().getStatistics().getNvoxel();
      clout << "Number of fluid voxels in simulation domain: " << nSimulation << std::endl;
      const unsigned nDesign = this->geometry().getStatistics().getNvoxel(6);
      clout << "Number of voxels in design domain: " << nDesign << std::endl;
    }
  }

  void prepareLattices() override
  {
    auto& lattice = this->lattice(NavierStokes());

    // PRIMAL
    if constexpr ( (MODE == SolverMode::Reference) || (MODE == SolverMode::Primal) ) {
      auto bulkIndicator = this->geometry().getMaterialIndicator({1, 3, 4, 6});
      lattice.template defineDynamics<PorousBGKdynamics>(bulkIndicator);

      boundary::set<boundary::BounceBack>(lattice, this->geometry(), 2);

      boundary::set<boundary::LocalVelocity>(
        lattice,
        this->geometry().getMaterialIndicator({3}));
      boundary::set<boundary::LocalPressure>(
        lattice,
        this->geometry().getMaterialIndicator({4}));
    }
    else {   // DUAL
      auto bulkIndicator = this->geometry().getMaterialIndicator({1, 6});
      lattice.template defineDynamics<DualPorousBGKDynamics>(bulkIndicator);

      auto bounceBackIndicator = this->geometry().getMaterialIndicator({2, 3, 4});
      boundary::set<boundary::BounceBack>(lattice, bounceBackIndicator);
    }
    const T omega = this->converter().getLatticeRelaxationFrequency();
    lattice.template setParameter<descriptors::OMEGA>(omega);

    AnalyticalConst3D<T,T> omegaF(omega);
    lattice.template defineField<OMEGA>(
      this->geometry().getMaterialIndicator({1,2,3,4,6}),
      omegaF);
  }

  void setInitialValues() override
  {
    std::unique_ptr<SuperIndicatorF3D<T>> bulkIndicator;
    if constexpr ( (MODE == SolverMode::Reference) || (MODE == SolverMode::Primal) ) {
      bulkIndicator = this->geometry().getMaterialIndicator({1, 3, 4, 6});
    }
    else {
      bulkIndicator = this->geometry().getMaterialIndicator({1, 2, 3, 4, 6});
    }

    AnalyticalConst3D<T,T> rhoF(1);
    std::vector<T> velocity(3,T());
    AnalyticalConst3D<T,T> uF(velocity);

    this->lattice(NavierStokes()).defineRhoU( bulkIndicator, rhoF, uF );
    this->lattice(NavierStokes()).iniEquilibrium( bulkIndicator, rhoF, uF);


    AnalyticalConst3D<T,T> one(1.);
    this->lattice(NavierStokes()).template defineField<POROSITY>(bulkIndicator, one);

    if constexpr ( MODE == SolverMode::Reference ) {
      // solid cube is the middle
      AnalyticalConst3D<T,T> zero(0.);
      this->lattice(NavierStokes()).template defineField<POROSITY>(
        this->geometry(),
        *this->parameters(Simulation()).simulationObject,
        zero);
    } else {
      this->lattice(NavierStokes()).template defineField<POROSITY>(
        this->parameters(Opti()).designDomain,
        *(this->parameters(Opti()).controlledField));
    }

    // load primal data
    // this needs to be done in every adjoint optimization
    if constexpr ( MODE == SolverMode::Dual ) {
      this->loadPrimalPopulations();
    }
  }

  void setBoundaryValues(std::size_t iT) override
  {
    const std::size_t itStart = this->converter().getLatticeTime(this->parameters(Simulation()).startUpTime);

    if (iT <= itStart) {
      PolynomialStartScale<T,std::size_t> StartScale( itStart, T( 1 ) );
      std::size_t iTvec[1] = {iT};
      T frac[1] = {};
      StartScale( frac,iTvec );

      if constexpr ( (MODE == SolverMode::Reference) || (MODE == SolverMode::Primal) ) {
        AnalyticalConst3D<T,T> uF(0., frac[0] * this->converter().getCharLatticeVelocity(), 0.);
        this->lattice(NavierStokes()).defineU(this->geometry(), 3, uF);
      }
    }
  }

  // prepare unsteady vtk output for reference simulation
  void prepareVTK() const override
  {
    if constexpr (MODE == SolverMode::Reference) {
      DomainIdSolver::LbSolver::prepareVTK();
    }
  }

  // (unsteady) vtk output for reference simulation
  void writeVTK(std::size_t iT) const override
  {
    if constexpr (MODE == SolverMode::Reference) {
      // ARTIFICIAL DATA: write solution (unsteady, until convergence in time is obtained)
      SuperVTMwriter3D<T> writer(this->parameters(VisualizationVTK()).filename);

      SuperLatticePhysVelocity3D<T,descriptor> velocity(this->lattice(NavierStokes()), this->converter());
      SuperLatticePhysPressure3D<T,descriptor> pressure(this->lattice(NavierStokes()), this->converter());
      SuperLatticePorosity3D<T,descriptor> porosity(this->lattice(NavierStokes()));
      writer.addFunctor(porosity);
      writer.addFunctor(velocity);
      writer.addFunctor(pressure);
      writer.write(iT);
    }
  }

public:
  /// Optimization: write stationary flow field once after each simulation
  // This method is called by the optimizer once at each successful optimization step
  // (only for the primal simulation)
  void postProcess() override
  {
    if constexpr (MODE == SolverMode::Primal) {
      if (this->parameters(names::VisualizationVTK()).output != "off") {
        if (this->parameters(OutputOpti()).counterOptiStep == 0) {
          DomainIdSolver::LbSolver::prepareVTK();
        }

        SuperVTMwriter3D<T> writer(this->parameters(VisualizationVTK()).filename);
        SuperLatticePhysVelocity3D<T,descriptor> velocity(this->lattice(NavierStokes()), this->converter());
        SuperLatticePhysPressure3D<T,descriptor> pressure(this->lattice(NavierStokes()), this->converter());
        SuperLatticePorosity3D<T,descriptor> porosity(this->lattice(NavierStokes()));
        writer.addFunctor(porosity);
        writer.addFunctor(velocity);
        writer.addFunctor(pressure);

        writer.write(this->parameters(OutputOpti()).counterOptiStep);
      }
    }
  }
};


/// @brief Objective and derivative computation
/// Objective = 0.5 * norm(u-uReference)^2 / norm(u)^2
/// @tparam T floating point type
template <typename T>
class RelativeDifferenceVelocityObjective : public DistributedObjective<T,DomainIdSolver>
{
protected:
  std::shared_ptr<DomainIdSolver<T,SolverMode::Reference>> _referenceSolver;
  std::shared_ptr<SuperF3D<T,T>>                           _referenceState;
  T _uRefNorm;
  T _cellVolume;
  T _regAlpha {0};

public:
  RelativeDifferenceVelocityObjective(XMLreader const& xml) {
    // create reference solver, compute reference solution
    _referenceSolver = createLbSolver <DomainIdSolver<T,SolverMode::Reference>> (xml);
    _referenceSolver->solve();
    _referenceState = std::make_shared<SuperLatticePhysVelocity3D<T,Descriptor>>(
      _referenceSolver->lattice(),
      _referenceSolver->converter());

    // define indicators
    this->_objectiveDomain = _referenceSolver->parameters(Opti()).objectiveDomain;
    this->_designDomain = _referenceSolver->parameters(Opti()).designDomain;

    // define constants
    const int dummy[1]{0};
    SuperL2Norm3D<T>(_referenceState, this->_objectiveDomain)(&_uRefNorm, dummy);
    _cellVolume = util::pow(
      this->_objectiveDomain->getSuperGeometry().getCuboidDecomposition().getDeltaR(), 3);
    xml.readOrWarn<T>("Optimization", "RegAlpha", "", _regAlpha);
  }

  T evaluate() override {
    // compare primal velocity with reference solution
    SuperLatticePhysVelocity3D<T,Descriptor> primalVelocity(
      this->_primalSolver->lattice(),
      this->_primalSolver->converter());
    RelativeDifferenceObjective3D<T,Descriptor> J(
      this->_primalSolver->lattice(),
      &primalVelocity,
      _referenceState,
      this->_objectiveDomain);
    int dummy[4]{0}; T res{0};
    J(&res, dummy);

    // regularization
    SuperLatticeSerialDataF<T,Descriptor> r(
      this->_primalSolver->lattice(),
      *(this->_controller),
      1,
      this->_serializer,
      [](T x) { return x*x; });
    SuperIntegral3D<T> R(&r, this->_designDomain);  // contains weighting by cell volume
    T reg{0};
    R(&reg, dummy);

    return res + 0.5*_regAlpha*reg;
  }

  std::shared_ptr<SuperF3D<T,T>> derivativeByPopulations() override {
    std::shared_ptr<SuperF3D<T,T>> primalVelocity
     = std::make_shared<SuperLatticePhysVelocity3D<T,Descriptor>>(
      this->_primalSolver->lattice(),
      this->_primalSolver->converter());
    std::shared_ptr<SuperF3D<T,T>> dPhysVelocityDf
     = std::make_shared<SuperLatticeDphysVelocityDf3D<T,Descriptor>> (
      this->_primalSolver->lattice(),
      this->_primalSolver->converter());
    std::shared_ptr<SuperF3D<T,T>> dObjectiveDf
     = std::make_shared<DrelativeDifferenceObjectiveDf3D<T,Descriptor>>(
      this->_primalSolver->lattice(),
      primalVelocity,
      dPhysVelocityDf,
      _referenceState,
      this->_objectiveDomain);
    return dObjectiveDf;
  }

  std::shared_ptr<SuperF3D<T,T>> derivativeByControl() override {
    return std::make_shared<SuperLatticeSerialDataF<T,Descriptor>>(
      this->_primalSolver->lattice(),
      *(this->_controller),
      1,
      this->_serializer,
      [this](T x) { return _regAlpha*_cellVolume*x; });
  }
};

/// @brief Objective and derivative computation
/// Objective = integral over 0.5 * norm(u-uReference)^2 / norm(u)^2 + 0.5*alpha*ctrl^2
///           = integral_{objectiveDomain} j(f,x) dx + integral_{designDomain} r(ctrl,x) dx
/// @tparam T floating point type
// This represents an alternative approach to the above class
// RelativeDifferenceVelocityObjective. Here, it is sufficient to define the
// objective via (pointwise) functions j and r. The spatial integration and the
// partial derivatives are then computed by the class GenericObjective
template <typename T>
class RelativeDifferenceVelocityObjectiveGeneric
{
protected:
  std::shared_ptr<DomainIdSolver<T,SolverMode::Reference>> _referenceSolver;
  std::shared_ptr<SuperF3D<T,T>>                           _referenceState;
  T _uRefNorm;
  T _conversion;
  T _regAlpha {0};

public:
  // these indicators are used by the GenericObjective class
  std::shared_ptr<SuperIndicatorF<T,3>>           _objectiveDomain;
  std::shared_ptr<SuperIndicatorF<T,3>>           _designDomain;

  RelativeDifferenceVelocityObjectiveGeneric(XMLreader const& xml) {
    // create reference solver, compute reference solution
    _referenceSolver = createLbSolver <DomainIdSolver<T,SolverMode::Reference>> (xml);
    _referenceSolver->solve();
    _referenceState = std::make_shared<SuperLatticePhysVelocity3D<T,Descriptor>>(
      _referenceSolver->lattice(),
      _referenceSolver->converter());

    // define indicators
    _objectiveDomain = _referenceSolver->parameters(Opti()).objectiveDomain;
    _designDomain = _referenceSolver->parameters(Opti()).designDomain;

    // define constants
    const int dummy[1]{0};
    SuperL2Norm3D<T>(_referenceState, _objectiveDomain)(&_uRefNorm, dummy);
    _conversion = _referenceSolver->converter().getConversionFactorVelocity();
    xml.readOrWarn<T>("Optimization", "RegAlpha", "", _regAlpha);
  }

  /// Implements the population-dependent term of objective functional
  // This interface is expected by the GenericObjective class. It uses different
  // cell/data types CELL/ V. Hence, beware of typecasts from V to T!
  template <typename CELL, typename V>
  void j(V* output, CELL& cell, int iC, const int* coords) {
    // application-specific implementation
    // compute primal velocity
    Vector<V,3> u;
    cell.computeU(u.data());
    // compute reference velocity
    T uRef[3];
    auto refCell = _referenceSolver->lattice().getBlock(iC).get(coords);  // expects local iC
    refCell.computeU(uRef);
    // compare
    const V diff = norm_squared(u - Vector<V,3>(uRef));
    // scale and return
    output[0] = diff * _conversion * _conversion / (T(2)*_uRefNorm*_uRefNorm);
  }

  /// Implements the control-dependent term of objective functional
  // This interface is expected by GenericObjective, with different data types V.
  // Hence, do not typecasts from V to T!
  template <typename V>
  void r(V* output, int iC, const int* coords, const V* control) {
    // application-specific implementation: compute regularization
    output[0] = T(0.5) * _regAlpha * control[0] * control[0];
  }
};



int main(int argc, char **argv)
{
  initialize(&argc, &argv);
  XMLreader config("parameter.xml");
  OstreamManager clout("main");

  using T = FLOATING_POINT_TYPE;

  OptiCaseDual<T,DomainIdSolver,descriptors::POROSITY,PorousBGKdynamics> optiCase(config);

  // classical objective from functors
  //auto objective = std::make_shared<RelativeDifferenceVelocityObjective<T>>(config);
  // alternative: generic objective computation
  auto objectiveHelp = std::make_shared<RelativeDifferenceVelocityObjectiveGeneric<T>>(config);
  auto objective = std::make_shared<GenericObjective<T,DomainIdSolver,
     RelativeDifferenceVelocityObjectiveGeneric<T>,
     PorousBGKdynamics, 1>>(objectiveHelp);
  optiCase.setObjective(objective);

  auto optimizer = createOptimizerLBFGS<T>(config, optiCase._dimCtrl);

  T startValue;
  config.readOrWarn<T>("Optimization", "StartValue", "", startValue);
  startValue = projection::getInitialControl(startValue, optiCase);
  optimizer->setStartValue(startValue);

  optimizer->optimize(optiCase);
}

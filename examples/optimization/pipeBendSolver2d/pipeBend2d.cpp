/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Shota Ito, Felix Schuhmann, Julius Jeßberger
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
 * In this example, the pipebend optimization benchmark case is implemented,
 * cf. e.g. https://doi.org/10.1016/j.buildenv.2024.112508.
 * For given in- and outlets at the left and at the bottom side of a cavity,
 * the (curved) shape of a pipe has to be found, s.t. either the overall
 * pressure drop or the overall dissipation is minimized.
 *
 * We use adjoint optimization for this. Hence, we implement
 * - a solver class, which implements a (primal or dual) simulation
 * - parameter structures, which keep simulation- and otimization-related parameters
 * - two objective classes, which implement the objective computation (one for
 *   dissipation minimization, one for pressure drop minimization)
 *
 * Since, for low Reynolds numbers, leaving the complete design domain open leads to
 * a low dissipation, we penalize the volume of the pipe in order to prevent this.
 * So, we can establish a geometrically meaningful geometry.
 *
 * The values of simulation- and optimization-related parameters are set in the
 * corresponding file parameters.xml. Feel free to change the Reynolds number and
 * see how the optimal pipe shape changes. This may of course necessitate adjustments
 * of e.g. the resolution.
 */

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

using descriptor = descriptors::DualPorousD2Q9Descriptor;
using Lattices = meta::map<
  NavierStokes, descriptor
>;

// select optimization mode: if true, then dissipation is minimized, else the
// pressure drop is minimized
constexpr bool dissipationObjective = false;

// keeps all simulation-related parameters
template <typename T, typename LATTICES>
struct PipeBendSimulationParameters
 : public parameters::XmlSimulation<T,LATTICES>
{
  T L;
  T lengthX,lengthY;
  T inflowY,inflowRadius,inflowLength;
  T outflowX,outflowRadius,outflowLength;
  T startValue;
};

// read the above PipeBendSimulationParameters from xml file
template<typename T, typename LATTICES, typename TAG>
struct parameters::Reader<PipeBendSimulationParameters<T,LATTICES>, TAG>
 : public parameters::ReaderBase<PipeBendSimulationParameters<T,LATTICES>>
{
  using ReaderBase<PipeBendSimulationParameters<T,LATTICES>>::ReaderBase;

  void read(XMLreader const& xml) {
    using namespace parameters;
    Reader<XmlSimulation<T,LATTICES>, TAG>(this->params).read(xml);
    auto params = this->params;
    auto converter = params->converter;

    // Fixed geometry parameters
    params->L = converter->getPhysDeltaX();
    params->lengthX = converter->getCharPhysLength() + 2*params->L;
    params->lengthY = params->lengthX;
    params->inflowY = 0.8*converter->getCharPhysLength() + params->L;
    params->inflowRadius = 0.1*converter->getCharPhysLength();
    params->inflowLength = 0.5*params->lengthX;
    params->outflowX = params->inflowY;
    params->outflowRadius = params->inflowRadius;
    params->outflowLength = params->inflowLength;
    params->startValue = xml["StartValue"].get<T>();
  }
};

/// bundle all parameter structs in a map
template<typename T, SolverMode MODE>
using Parameters = meta::map<
  Simulation,       PipeBendSimulationParameters<T,Lattices>,
  Output,           parameters::OutputGeneral<T>,
  VisualizationVTK, parameters::OutputPlot<T>,
  Opti,             parameters::DistributedOptiSimulation<T,Lattices,MODE>,
  Results,          parameters::DistributedOptiSimulationResults<T,Lattices>,
  OutputOpti,       parameters::OptiOutput<T>
>;

/// @brief Implements the simulation
/// @tparam T floating point type
/// @tparam MODE either Primal of Dual (adjoint)
template<typename T,SolverMode MODE>
class PipeBendOptiSolver : public AdjointLbSolver<T,Parameters<T,MODE>,Lattices,MODE>
{
private:
  mutable OstreamManager            clout {std::cout, "PipeBendOptiSolver"};

public:
  PipeBendOptiSolver(utilities::TypeIndexedSharedPtrTuple<Parameters<T,MODE>> params)
   : PipeBendOptiSolver::LbSolver(params),
    PipeBendOptiSolver::AdjointLbSolver(params)
  {
    const std::string name = (MODE == SolverMode::Primal) ? "Flow_Simulations" : "Dual_Simulations";
    this->parameters(VisualizationVTK()).filename = name;

    if constexpr (MODE == SolverMode::Dual) {
      // since Dirichlet BC are set for the dual problem, we need to re-calibrate the pressure there
      this->parameters(Simulation()).pressureFilter = true;
    }
  }

protected:
  void prepareGeometry() override
  {
    // Create necessary instances for geometry
    auto& params = this->parameters(Simulation());

    Vector<T,2> origin(-params.inflowLength, -params.outflowLength);
    Vector<T,2> extend(params.lengthX + params.inflowLength, params.lengthY + params.outflowLength);
    IndicatorCuboid2D<T> extendedPipeBendDomain(extend, origin);
    this->_cGeometry = std::make_shared<CuboidDecomposition2D<T>>(extendedPipeBendDomain, params.L, params.noC);
    this->_loadBalancer = std::make_shared<HeuristicLoadBalancer<T>>(*this->_cGeometry);
    this->_sGeometry = std::make_shared<SuperGeometry<T,2>>(*this->_cGeometry,
                                                            *this->_loadBalancer,
                                                            params.overlap);
    auto& superGeometry = this->geometry();
    origin = Vector<T,2>{0., 0.};
    extend = Vector<T,2>{ params.lengthX, params.lengthY};
    auto square = std::make_shared<IndicatorCuboid2D<T>> ( extend, origin );
    // cut away two corners in order to reduce the number of control variables
    // the size of the corners should be adjusted to the Reynolds number
    auto t1 = std::make_shared<IndicatorTriangle2D<T>> (Vector<T,2>( 0, 0 ),
      Vector<T,2>( 0.5 * params.lengthX, 0 ),
      Vector<T,2>( 0, 0.5 * params.lengthY ));
    auto t2 = std::make_shared<IndicatorTriangle2D<T>> (Vector<T,2>( params.lengthX, params.lengthY ),
      Vector<T,2>( 0.5 * params.lengthX, params.lengthY ),
      Vector<T,2>( params.lengthX, 0.5 * params.lengthY ));
    std::shared_ptr<IndicatorF2D<T>> pipeBendDomain = square - t1 - t2;

    origin = Vector<T,2>( -params.inflowLength -.5*params.L, params.inflowY - params.inflowRadius );
    extend = Vector<T,2>( params.L, 2. * params.inflowRadius);
    IndicatorCuboid2D<T> inflow( extend, origin );
    origin[0] = -.5*params.L;
    std::shared_ptr<IndicatorF2D<T>> objectiveDomainInflow
     = std::make_shared<IndicatorCuboid2D<T>>( extend, origin );
    origin = Vector<T,2>( -params.inflowLength -.5*params.L, params.inflowY - params.inflowRadius - params.L );
    extend = Vector<T,2>( params.inflowLength + params.L, 2.*params.inflowRadius + 2*params.L );
    IndicatorCuboid2D<T> extendedInflow( extend, origin );

    origin = Vector<T,2>( params.outflowX - params.outflowRadius, -params.outflowLength -.5*params.L );
    extend = Vector<T,2>( 2. * params.outflowRadius, params.L );
    IndicatorCuboid2D<T> outflow( extend, origin );
    origin[1] = -.5*params.L;
    std::shared_ptr<IndicatorF2D<T>> objectiveDomainOutflow
     = std::make_shared<IndicatorCuboid2D<T>>( extend, origin );
    origin = Vector<T,2>( params.outflowX - params.outflowRadius - params.L, -params.outflowLength -.5*params.L );
    extend = Vector<T,2>( 2.*params.outflowRadius + 2*params.L, params.outflowLength + params.L );
    IndicatorCuboid2D<T> extendedOutflow( extend, origin );

    // Material number for outer channel wall = 2
    superGeometry.rename( 0, 2, pipeBendDomain );
    superGeometry.rename( 0, 2, extendedInflow );
    superGeometry.rename( 0, 2, extendedOutflow );

    // Material number for flow domain = 1
    superGeometry.rename( 2, 1, {1,1} );

    // Material number for inflow = 3, outflow = 4
    superGeometry.rename( 2,3, inflow );
    superGeometry.rename( 2,4, outflow );

    // Material number for design domain (where the porosity is variable) = 6
    superGeometry.rename( 1,6, pipeBendDomain );
    superGeometry.rename( 6,7, objectiveDomainInflow );  // for pressure drop objective
    superGeometry.rename( 6,8, objectiveDomainOutflow );

    // store some geometric data for optimization
    this->parameters(Results()).geometry = this->_sGeometry;
    this->parameters(Opti()).bulkIndicator = std::move(this->geometry().getMaterialIndicator({1,6,7,8}));
    this->parameters(Opti()).designDomain = std::move(this->geometry().getMaterialIndicator({6}));
    if constexpr (dissipationObjective) {
      this->parameters(Opti()).objectiveDomain = std::move(this->geometry().getMaterialIndicator({6}));
    } else {
      this->parameters(Opti()).objectiveDomain = std::move(this->geometry().getMaterialIndicator({7,8}));
    }
  }

  void prepareLattices() override
  {
    auto& superGeometry = this->geometry();
    auto& sLattice = this->lattice();

    if constexpr (MODE == SolverMode::Primal) {
      auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,6,7,8});
      sLattice.template defineDynamics<PorousBGKdynamics<T,descriptor>>(bulkIndicator);
      boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
      boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
      boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);
    }
    else if (MODE == SolverMode::Dual){
      auto bulkIndicator = superGeometry.getMaterialIndicator({1,6,7,8});
      sLattice.template defineDynamics<DualPorousBGKDynamics<T,descriptor>>(bulkIndicator);

      auto wallIndicator = superGeometry.getMaterialIndicator({2,3,4});
      boundary::set<boundary::BounceBack>(sLattice, wallIndicator);
    }
    const T omega = this->converter().getLatticeRelaxationFrequency();
    sLattice.template setParameter<descriptors::OMEGA>(omega);
  }

  void setInitialValues() override
  {
    auto& superGeometry = this->geometry();
    auto& sLattice = this->lattice();
    auto& optiParams = this->parameters(Opti());

    // Initial conditions
    AnalyticalConst2D<T,T> rhoF(1);
    Vector<T,2> velocity;
    AnalyticalConst2D<T,T> uF(velocity);

    // Initialize all values of distribution functions to their local equilibrium
    auto all = superGeometry.getMaterialIndicator({0,1,2,3,4,6,7,8});
    sLattice.defineRhoU(all, rhoF, uF);
    sLattice.iniEquilibrium(all, rhoF, uF);

    // Set porosity field
    auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,6,7,8});
    auto solidIndicator = superGeometry.getMaterialIndicator({0,2});
    AnalyticalConst2D<T,T> zero(0.);
    AnalyticalConst2D<T,T> one(1.);

    sLattice.template defineField<descriptors::POROSITY>(bulkIndicator, one);
    sLattice.template defineField<descriptors::POROSITY>(solidIndicator, zero);
    sLattice.template defineField<descriptors::POROSITY>(optiParams.designDomain, *(optiParams.controlledField));

    if constexpr (MODE == SolverMode::Dual) {
      this->loadPrimalPopulations();  // coupling from primal to dual simulation
    }
  }

  void setBoundaryValues(std::size_t iT) override
  {
    auto& converter = this->converter();
    auto& params = this->parameters(Simulation());

    // No of time steps for smooth start-up
    std::size_t iTmaxStart = converter.getLatticeTime(params.startUpTime);
    std::size_t iTupdate = converter.getLatticeTime(params.physBoundaryValueUpdateTime);
    if constexpr(MODE == SolverMode::Primal) {
      if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
        // Smooth start curve, polynomial
        PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );

        // Creates and sets the Poiseuille inflow profile using functors
        T iTvec[1] = {T( iT )};
        T frac[1] = {};
        StartScale( frac,iTvec );
        T maxVelocity = converter.getCharLatticeVelocity()*3./2.*frac[0];
        T distance2Wall = params.L/2.;
        Poiseuille2D<T> poiseuilleU(this->geometry(), 3, maxVelocity, distance2Wall);

        this->lattice().defineU(this->geometry(), 3, poiseuilleU);
        this->lattice().template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
          ProcessingContext::Simulation);
      }
    }
  }

  void prepareVTK() const override
  {
    if (this->parameters(OutputOpti()).counterOptiStep == 0) {
      PipeBendOptiSolver::LbSolver::prepareVTK();
    }
  }

  void writeVTK(std::size_t iT) const override
  {
    SuperVTMwriter2D<T> writer(this->parameters(VisualizationVTK()).filename);
    const int index = this->parameters(OutputOpti()).counterOptiStep;

    auto& converter = this->converter();
    auto& superGeometry = this->geometry();
    auto& sLattice = this->lattice(NavierStokes());

    SuperLatticePhysVelocity2D<T,descriptor> velocity(sLattice, converter);
    SuperLatticePhysPressure2D<T,descriptor> pressure(sLattice, converter);
    SuperLatticePorosity2D<T,descriptor> porosity(sLattice, superGeometry, 1, converter);
    SuperGeometryF2D<T> geometry(superGeometry);

    SuperLatticeFfromCallableF<T,descriptor> viscousDissipation(sLattice, [&](T* output, auto cell){
      T rho, uTemp[descriptor::d], pi[util::TensorVal<descriptor>::n];
      cell.computeAllMomenta(rho, uTemp, pi);
      const T porosity = cell.template getField<descriptors::POROSITY>();

      T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
      if (util::TensorVal<descriptor>::n == 6) {
        PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                        + pi[5] * pi[5];
      }
      const T omega = T(1) / converter.getLatticeRelaxationTime();
      const T dt = converter.getPhysDeltaT();
      if (porosity > 0) {
        output[0] = PiNeqNormSqr * util::pow(omega * descriptors::invCs2<T,descriptor>() / rho / dt, 2) / 2.
              * converter.getPhysViscosity();
      } else {
        output[0] = 0;
      }
    });
    viscousDissipation.getName() = "viscousDissipation";
    SuperLatticeFfromCallableF<T,descriptor> porousDissipation(sLattice, [&](T* output, auto cell){
      T uTemp[descriptor::d];
      cell.computeU(uTemp);
      const T porosity = cell.template getField<descriptors::POROSITY>();

      const T invPermeability = projection::porosityToInvPermeability(porosity, converter);
      const T uNormSq = util::euklidN2(uTemp, descriptor::d);
      output[0] = converter.getPhysViscosity() * invPermeability * uNormSq;
    });
    porousDissipation.getName() = "porousDissipation";

    writer.addFunctor(velocity);
    writer.addFunctor(pressure);
    writer.addFunctor(porosity);
    writer.addFunctor(geometry);
    writer.addFunctor(viscousDissipation);
    writer.addFunctor(porousDissipation);
    writer.write(index);
  }

  void computeResults() override
  {
    // some additional evaluation
    if constexpr (MODE == SolverMode::Primal) {
      if constexpr (! dissipationObjective) {
        computeDissipation();  // print dissipation

        // print fraction of cavity which is occupied by pipe
        T help {0}; int tmp[1] = {0};
        SuperLatticePorosity2D<T,descriptor> porosity(
          this->lattice(), this->geometry(), 6, this->converter());
        SuperL1Norm2D<T> totalVolume(porosity, this->geometry(), 6);
        totalVolume(&help, tmp);
        clout << "volume-fraction=" << help / 0.25 << std::endl;
      }
      computePressureDrop();  // print pressure drop
    }

    // store super lattice for optimization processing
    PipeBendOptiSolver::AdjointLbSolver::computeResults();
  }

protected:
  void computePressureDrop()
  {
    // pressure drop
    SuperLatticePhysPressure2D<T,descriptor> pressure( this->lattice(), this->converter() );
    SuperAverage2D<T,T> objectiveInflow( pressure, this->geometry(), 7 );
    SuperAverage2D<T,T> objectiveOutflow( pressure, this->geometry(), 8 );

    int input[2] = {}; T pressure1[1], pressure2[1];
    objectiveInflow( pressure1, input );
    objectiveOutflow( pressure2, input );
    const T pressureDrop = pressure1[0] - pressure2[0];

    clout << "pressure-in=" << pressure1[0];
    clout << "; pressure-out=" << pressure2[0];
    clout << "; pressure-drop=" << pressureDrop << std::endl;
  }

  void computeDissipation()
  {
    // viscous dissipation
    int tmp[1] = {0}; T result[1];
    SuperLatticePhysDissipation2D<T,descriptor> vDiss(this->lattice(), this->converter());

    SuperL1Norm2D<T> vdL1Norm(vDiss, this->geometry(), 6);
    vdL1Norm(result, tmp);
    clout << "viscous-dissipation=" << result[0] << std::endl;

    // porous dissipation
    const T convVelocity = this->converter().getConversionFactorVelocity();
    SuperLatticeFfromCallableF<T,descriptor> pDiss(
      this->lattice(),
      [&](T* output, auto cell){
        T uTemp[descriptor::d];
        cell.computeU(uTemp);
        const T porosity = cell.template getField<descriptors::POROSITY>();
        const T invPermeability = projection::porosityToInvPermeability(porosity, this->converter());
        const T uNormSq = util::euklidN2(uTemp, descriptor::d) * convVelocity * convVelocity;
        output[0] = this->converter().getPhysViscosity() * invPermeability * uNormSq;
      });
    SuperL1Norm2D<T> pdL1Norm(pDiss, this->geometry(), 6);
    pdL1Norm(result, tmp);
    clout << "porous-dissipation=" << result[0] << std::endl;
  }
};



/// @brief Objective = pressure drop + weight * pipeVolume
/// @tparam T floating point type
template<typename T>
class PressureObjectiveGeneric
{
private:
  // plain instance of a solver, it helps to instance material numbers etc.
  std::shared_ptr<PipeBendOptiSolver<T,SolverMode::Primal>> _refSolver;

  std::shared_ptr<SuperIndicatorF<T,descriptor::d>> _objectiveDomainIn;
  std::shared_ptr<SuperIndicatorF<T,descriptor::d>> _objectiveDomainOut;

  T _convPressure;  // conversion into physical units
  T _scaling;       // transforms integral of pressure into average pressure at cross-section (approximately)
  T _regAlpha {0};  // regularization factor

public:
  std::shared_ptr<SuperIndicatorF<T,descriptor::d>> _objectiveDomain;
  std::shared_ptr<SuperIndicatorF<T,descriptor::d>> _designDomain;

  PressureObjectiveGeneric(XMLreader const& xml) {
    // create reference solver
    _refSolver = createLbSolver<PipeBendOptiSolver<T,SolverMode::Primal>>(xml);
    _refSolver->initialize();

    // set indicators
    _objectiveDomainIn = std::make_shared<SuperIndicatorMaterial2D<T>>(
      _refSolver->geometry(), std::vector<int>{7});
    _objectiveDomainOut = std::make_shared<SuperIndicatorMaterial2D<T>>(
      _refSolver->geometry(), std::vector<int>{8});
    _objectiveDomain = std::make_shared<SuperIndicatorMaterial2D<T>>(
      _refSolver->geometry(), std::vector<int>{7, 8});
    _designDomain = _refSolver->parameters(Opti()).designDomain;

    // set constants
    _convPressure = _refSolver->converter().getConversionFactorPressure();
    _scaling = T(0.5) / (_refSolver->parameters(Simulation()).inflowRadius
      * util::pow(_refSolver->converter().getPhysDeltaX(), descriptor::d-1));
    xml.readOrWarn<T>("Optimization", "RegAlpha", "", _regAlpha);
  }

  template<typename CELL, typename V>
  void j(V* output, CELL& cell, int iC, const int* coords) {
    const int coordsHelp[descriptor::d+1] = {_refSolver->geometry().getLoadBalancer().glob(iC), coords[0], coords[1]};
    if ((*_objectiveDomainIn)(coordsHelp)) {  // at inlet
      output[0] = pressure(cell);
    } else if ((*_objectiveDomainOut)(coordsHelp)) {  // at outlet
      output[0] = - pressure(cell);
    } else {
      output[0] = V(0);
    }
  }

  template<typename V>
  void r(V* output, int iC, const int* coords, const V* control) {
    // cell-wise contribution to the volume of the pipe
    // weighted by volume of the cavity, which is appr. 0.25
    // porosity = Sigmoid(control)
    output[0] = V(_regAlpha / 0.25) * projection::Sigmoid<V>().project(control[0]);
  }

private:
  template <typename CELL, typename V = typename CELL::value_t>
  V pressure(CELL& cell) {
    const V latticePressure = ( cell.computeRho() - V(1) ) / descriptors::invCs2<V,descriptor>();
    return latticePressure * V(_convPressure * _scaling);
  }
};


/// Objective = total dissipation + weight * pipeVolume
/// @tparam T floating point type
// Remark: using the GenericObjective helper class is difficult here, since we
// have terms which involve both state and control. So we have to implement the
// partial derivatives by hand.
template <typename T>
class DissipationObjective : public DistributedObjective<T,PipeBendOptiSolver>
{
private:
  std::shared_ptr<PipeBendOptiSolver<T,SolverMode::Primal>> _refSolver;
  std::shared_ptr<UnitConverter<T,descriptor>> _converter;
  T _regAlpha {0};
  T _viscosity;     // phys. viscosity (constant)
  T _convVelocity;  // conversion phys. velocity (constant)

public:
  DissipationObjective(XMLreader const& xml) {
    // create reference solver
    _refSolver = createLbSolver <PipeBendOptiSolver<T,SolverMode::Primal>> (xml);
    _refSolver->initialize();

    // define indicators
    this->_objectiveDomain = _refSolver->parameters(Opti()).objectiveDomain;
    this->_designDomain = _refSolver->parameters(Opti()).designDomain;

    // define constants
    _converter = std::shared_ptr<UnitConverter<T,descriptor>>(createUnitConverter<T,descriptor>(xml));
    xml.readOrWarn<T>("Optimization", "RegAlpha", "", _regAlpha);
    _viscosity = _converter->getPhysViscosity();
    _convVelocity = _converter->getConversionFactorVelocity();
  }

  T evaluate() override {
    OstreamManager clout("Objective");
    T result {0}; T help {0};
    int tmp[1] = {0};
    SuperLatticeFfromCallableF<T,descriptor> dissipation(
      this->_primalSolver->lattice(),
      [&](T* output, auto cell){
      T rho, uTemp[descriptor::d], pi[util::TensorVal<descriptor>::n];
      cell.computeAllMomenta(rho, uTemp, pi);
      const T porosity = cell.template getField<descriptors::POROSITY>();

      T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
      if (util::TensorVal<descriptor>::n == 6) {
        PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                        + pi[5] * pi[5];
      }
      const T omega = T(1) / _converter->getLatticeRelaxationTime();
      const T dt = _converter->getPhysDeltaT();
      const T viscousDiss = (porosity > 0) ?
        PiNeqNormSqr * util::pow(omega * descriptors::invCs2<T,descriptor>() / rho / dt, 2) / T(2)
              * _converter->getPhysViscosity()
        : T(0);

      const T invPermeability = projection::porosityToInvPermeability(porosity, *_converter);
      const T uNormSq = util::euklidN2(uTemp, descriptor::d) * _convVelocity * _convVelocity;
      const T porousDiss = _viscosity * invPermeability * uNormSq;
      output[0] = viscousDiss + porousDiss;
    });
    SuperL1Norm2D<T> dL1Norm(dissipation,this->_primalSolver->geometry(),6);
    dL1Norm(&result,tmp);
    clout << "total-dissipation=" << result << std::endl;

    // regularization: penalize volume of pipe
    SuperLatticePorosity2D<T,descriptor> porosity(
      this->_primalSolver->lattice(),this->_primalSolver->geometry(),6,*_converter);
    SuperL1Norm2D<T> totalVolume(porosity,this->_primalSolver->geometry(),6);
    totalVolume(&help,tmp);
    clout << "volume-fraction=" << help / 0.25 << std::endl;
    result += _regAlpha * help / 0.25;
    clout << "objective=" << result << std::endl;
    return result;
  }

  std::shared_ptr<SuperF2D<T,T>> derivativeByPopulations() override {
    std::shared_ptr<SuperF2D<T,T>> dPhysDissipationDf
     = std::make_shared<SuperLatticeDphysDissipationDf<T,descriptor>>(
      this->_primalSolver->lattice(),
      *_converter);
    std::shared_ptr<SuperF2D<T,T>> dPorousDissipationDf
     = std::make_shared<SuperLatticeFfromCallableF<T,descriptor>>(
      this->_primalSolver->lattice(),
      [&](T* output, auto cell, int iCloc, const int* coords){
        // get permeability factor
        const T porosity = cell.template getField<descriptors::POROSITY>();
        const T invPermeability = projection::porosityToInvPermeability(porosity, *_converter);

        // get u-dependent term
        Vector<T,descriptor::d> u;
        cell.computeU(u.data());
        Vector<T,descriptor::d> dudf;
        for (unsigned jPop=0; jPop<descriptor::q; ++jPop) {
          opti::dualLbMomentaHelpers<descriptor>::dUDf(cell, dudf, jPop);
          output[jPop] = _viscosity * invPermeability * T(2) * (u * dudf) * _convVelocity * _convVelocity;
        }
      });
    return dPhysDissipationDf + dPorousDissipationDf;
  }

  std::shared_ptr<SuperF2D<T,T>> derivativeByControl() override {
    return std::make_shared<SuperLatticeFfromCallableF<T,descriptor>>(
      this->_primalSolver->lattice(),
      [&](T* output, auto cell, int iCloc, const int* coords){
        // get control
        const int iCglob = _refSolver->geometry().getLoadBalancer().glob(iCloc);
        const int globalLatticeR[3] = {iCglob, coords[0], coords[1]};
        const auto index = this->_serializer->getSerializedCellIndex(globalLatticeR);
        const T control = this->_controller->getControl(index);

        // get other factors
        T uTemp[descriptor::d];
        cell.computeU(uTemp);
        const T uNormSq = util::euklidN2(uTemp, descriptor::d) * _convVelocity * _convVelocity;
        const T dInvKdPorosity = T(-1) / projection::gridTerm<T,descriptor>(*_converter);
        const T dPorosityDalpha = projection::Sigmoid<T>().derivative(control);

        // multiply (chain rule)
        output[0] = _viscosity * uNormSq * dInvKdPorosity * dPorosityDalpha
          + _regAlpha / T(0.25) * dPorosityDalpha;  // derivative of pipe volume term
      });
  }
};


int main(int argc, char **argv)
{
  initialize(&argc, &argv);
  XMLreader config("parameter.xml");
  using S = double;

  // create instances for optimization
  // creates solvers, organizes coupling between them, implements derivative computation
  OptiCaseDual<S,PipeBendOptiSolver,descriptors::POROSITY,PorousBGKdynamics> optiCase(config);

  std::shared_ptr<DistributedObjective<S,PipeBendOptiSolver>> objective;  // implements objective computation
  if constexpr (dissipationObjective) {
    objective = std::make_shared<DissipationObjective<S>>(config);
  } else {
    auto objectiveHelp = std::make_shared<PressureObjectiveGeneric<S>>(config);
    objective = std::make_shared<GenericObjective<S,PipeBendOptiSolver,
      PressureObjectiveGeneric<S>,PorousBGKdynamics, 1>>(objectiveHelp);
  }
  optiCase.setObjective(objective);

  auto optimizer = createOptimizerLBFGS<S>(config, optiCase._dimCtrl);  // implements search algorithm

  // compute start value for optimization
  S startValue;
  config.readOrWarn<S>("Optimization", "StartValue", "", startValue);
  startValue = projection::getInitialControl(startValue, optiCase);
  optimizer->setStartValue(startValue);

  // run optimization
  optimizer->optimize(optiCase);
}

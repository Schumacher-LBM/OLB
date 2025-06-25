#ifndef HELPER_H
#define HELPER_H

namespace olb {

// Temporary method, replace with SuperLatticeFieldReduction in future
template <typename T, typename DESCRIPTOR>
T norm(SuperLattice<T,DESCRIPTOR>& refLattice, auto& converter, auto& objectiveDomain) {
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> ref_u(refLattice, converter);
  int tmp[4]{0}; T refNorm[1]{0};
  SuperL2Norm3D<T>(ref_u,objectiveDomain)(refNorm, tmp);
  return refNorm[0]*refNorm[0];
}

// Temporary method, replace with SuperLatticeFieldReduction in future
template <typename FIELD, typename T, typename DESCRIPTOR>
auto integrate(SuperLattice<T,DESCRIPTOR>& lattice, auto& domain) {
  std::size_t constexpr DIM = DESCRIPTOR::template size<FIELD>();
  SuperLatticeField3D<T,DESCRIPTOR,FIELD> f(lattice);

  int dummy[4]{0}; T F[DIM]{0};
  SuperIntegral<DIM,T> integral(f, domain);
  integral(F, dummy);

  Vector<T,DIM> result;
  std::copy(std::begin(F), std::end(F), result.begin());
  return result;
}

}

#endif

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jeßberger, Jan Kösters
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
 * This example illustrates the usage of forward automatic differentiation in
 * order to evaluate partial derivatives.
 * It is structured as a tutorial that can be read sequentially.
 */

#include <olb.h>


using namespace olb;
using namespace olb::util;  // for ADf data type

using S = FLOATING_POINT_TYPE;

/** ---------------------------------------------------------------------------
 * Part 0: Preliminaries
 * ----------------------------------------------------------------------------
 * The computation of partial derivatives of a function or functor f relies on
 * calling f with the special datatype ADf. Hence, f must not be fixed to
 * accept (and return) a certain type, but instead be templatized or built on
 * automatic type deduction (auto type). Furthermore, explicit typecasts should
 * be omitted in the definition of f.
 *
 * In general, any (polymorphic) function/ functor f: x -> y can be called with
 * arguments x_0, ... x_n of type ADf. If the arguments are initialized s.t.
 * d/dx_j x_i = \delta_{ij}, then f computes the value as well as the partial
 * derivatives of f at x, s.t. we get y.v() = f(x) and
 * y[i].d(j) = d/dx_j f_i(x).
 * Directly using this overloading is the most efficient way to compute f and
 * its derivatives. The function iniAD helps initializing the derivatives (cf.
 * example below).
 * Since one sometimes only needs the derivative (and not the value) of f,
 * helper functions have been defined which directly yield the partial
 * derivatives as an array (cf. examples below).
 *
 * The framework is e.g. applicable to C++ functions, OpenLB functors and also
 * full LBM simulations. For each of these, examples are given below.
 */

/** ---------------------------------------------------------------------------
 * Part 1: Evaluating derivatives of functions
 * ----------------------------------------------------------------------------
 * First, we define some example functions. In general, any (polymorphic)
 * function can be applied to ADf-typed arguments.
 */

// Example 1: One-dimensional lambda function f: S -> T
auto p1 = [](auto x) {
  return x * x;
};

// Example 2: Two-dimensional lambda function f: S* -> T
auto p2 = [](auto x) {
  return x[0] * x[0] + x[1] + 1;
};

// Example 3: General multi-dimensional lambda function f: S* -> T*
// the argument y is modified in place
auto p3 = [](auto y, auto x){
  y[0] = x[0]+x[1];
  y[1] = 3.0 * x[0] * x[1];
  y[2] = x[1] * x[1];
};

// Example 4: One-dimensional free function f: S -> T
template <typename T>
T p4 (T x) {
  return x * x;
}

// Example 5: these extended math functions can be mulit-dimensional
// and can be combined with other functions
auto p5 = [](auto y, auto x){
  y[0] = exp(x[0]) + 2.0*x[1];
  y[1] = 3.0 * x[0] + pow(cos(x[1]), 2);
  y[2] = sqrt(x[1]+1);
};


/// Get derivatives by applying functions to the ADf versions of variables
void functionADfApplication() {
  OstreamManager clout(std::cout,"functionADfApplication");

  // Evaluate p1 at x=2.5
  clout << "Consider p1(x) = x * x" << std::endl;
  const S y0 {2.5};
  S y1 = p1 (y0);
  clout << "p1(2.5) = " << y1 << std::endl;

  // Get ADf version of the argument y0
  ADf<S,1> y0_ad = iniAD(y0);         // auto would be ok
  // Apply p1 to that ADf typed argument
  ADf<S,1> y1_ad = p1 (y0_ad);        // auto would be ok
  clout << "p1(2.5), p1'(2.5) = " << y1_ad << '\n' << std::endl;
  // The value of the result y1_ad can be accessed via "y1_ad.v()", the i-th
  // derivative via "y1_ad.d(i)" (here we only have i=0).
  assert((y1_ad.v() == 6.25));
  assert((y1_ad.d(0) == 5.0));

  // The same procedure works for p2, p3, p4

  // Exemplarily for p3
  clout << "Consider p3(x) = [x0 + x1, 3 * x0 * x1, x1 * x1]" << std::endl;
  const S z0[2] {1.0, 2.0};
  S u1[3];
  p3(u1, z0);
  clout << "p3(1.0, 2.0) = [" << u1[0] << ", "
    << u1[1] << ", " << u1[2] << "]" << std::endl;

  // Get ADf typed argument. We need to tell the array dimensions
  ADf<S,2>* z0_ad = iniAD<2>(z0);     // auto would be ok
  // Initialize output array
  ADf<S,2> u1_ad[3];
  // Apply the function to the ADf argument
  p3 (u1_ad, z0_ad);
  clout << "grad p3(1.0, 2.0) = " << u1_ad[0] << ", "
    << u1_ad[1] << ", " << u1_ad[2] << '\n' << std::endl;
}


/**
 * In the case that only the derivative is needed one can use the wrapper
 * "derivativeFAD", which accepts the "normal" arithmetic argument and directly
 * returns the value/ array of derivatives, depending on the dimensions.
 *
 * Some different interfaces are supported. f can either map value to value (as
 * does p1), or array to value (as does p2), or array to array with in-place
 * modification (p3).
 * For technical reasons, the function must be given as a lambda expression (so
 * the free function p4 needs to be wrapped in a lambda (p6)).
 */

// Example 6: Wrap p4 in a lambda expression
auto p6 = [](auto x) {
  return p4(x);
};


void functionDerivatives() {
  OstreamManager clout(std::cout,"functionDerivatives");

  // Evaluate p1 at x=2.5
  clout << "Consider p1(x) = x * x" << std::endl;
  const S y0 {2.5};
  S y1 = p1 (y0);
  clout << "p1(2.5) = " << y1 << std::endl;

  // Get the derivative of p1 at x= 2.5
  S y2 = derivativeFAD (p1, y0);
  clout << "p1'(2.5) = " << y2 << '\n' << std::endl;

  // This also works for p2 on multi-dimensional space
  clout << "Consider p2(x) = x0 * x0 + x1 + 1" << std::endl;
  const S z0[2] {1.0, 2.0};
  S z1 = p2 (z0);
  clout << "p2(1.0, 2.0) = " << z1 << std::endl;

  // but we need to tell the array size of x.
  // an array is assigned for the result
  auto z2 = derivativeFAD<2> (p2, z0);
  clout << "grad p2(1.0, 2.0) = [" << z2[0] << ", " << z2[1] << "]\n" << std::endl;

  // And for multi-dimensional function p3
  clout << "Consider p3(x) = [x0 + x1, 3 * x0 * x1, x1 * x1]" << std::endl;
  S u1[3];
  p3(u1, z0);
  clout << "p3(1.0, 2.0) = [" << u1[0] << ", "
    << u1[1] << ", " << u1[2] << "]" << std::endl;

  // we need to tell both source and target dimension to get the derivative.
  // the second argument will be modified in-place
  // indices are result[j*dim(source)+i] = d/dx_i (s_j)
  S u2[2*3];
  derivativeFAD<2> (p3, u2, z0, 3);
  clout << "grad p3(1.0, 2.0) = [" << u2[0] << ", "
    << u2[1] << ", " << u2[2] << ", " << u2[3] << ", "
      << u2[4] << ", " << u2[5] << "]\n" << std::endl;

  // For p6, the same results as for p1 can be expected
  // clout << "Consider p6(x) = x * x" << std::endl;
  const S v0 {2.5};
  S v1 = p6 (v0);
  clout << "p6(2.5) = " << v1 << std::endl;

  // Get the derivative of p6 at x= 2.5
  S v2 = derivativeFAD (p6, v0);
  clout << "p6'(2.5) = " << v2 << '\n' << std::endl;
}


/** ---------------------------------------------------------------------------
 * Part 2: Evaluating derivatives of functors
 * ----------------------------------------------------------------------------
 * Similar to functions, we can apply the automatic differentiation concept to
 * functors.
 */

/// If we construct and call a functor with ADf data types, we get values +
/// derivatives.
void functorADfApplication() {
  OstreamManager clout(std::cout,"functorADfApplication");

  // Consider a 1-dimensional analytical linear functor
  // Attention: Contrary to functions, functors have no return value, but the
  // first argument is modified in-place.
  clout << "Consider functor f1(x) = 2.1x + 1.5" << std::endl;
  AnalyticalLinear1D<ADf<S,1>,ADf<S,1>> f1(S(2.1), S(1.5));
  S s1[1] = {0.5};
  // Get ADf typed argument. We need to tell the array dimensions
  ADf<S,1>* s1_ad = iniAD<1>(s1);
  // Assign output array
  ADf<S,1> t1_ad[1];
  // Apply f0 to that ADf type argument
  f1(t1_ad, s1_ad);
  clout << "f1(0.5), f1'(0.5) = " << t1_ad[0] << '\n' << std::endl;

  //Consider analytical constant functor in 2 dimensions
  clout << "Consider functor f2(x) = [3.145193, 2.0]"
    << std::endl;
  AnalyticalConst2D<ADf<S,2>,ADf<S,2>> f2(S(3.145193), S(2.0));
  // Define argument
  S s2[2] = {-1.0, -1.0};
  // Get ADf typed argument. We need to tell the array dimensions
  ADf<S,2>* s2_ad = iniAD<2>(s2);
  // Initialize output array
  ADf<S,2> t2_ad[2];
  // Apply f2 to that ADf type argument
  f2(t2_ad, s2_ad);
  clout << "first component:  f2_1(-1,-1), grad f2_1(-1,-1) = "
    << t2_ad[0] << std::endl;
  clout << "second component: f2_2(-1,-1), grad f2_2(-1,-1) = "
    << t2_ad[1] << '\n' << std::endl;
}

/// If we are only interested in the derivatives, we can make use of a wrapper
/// and therefore avoid the ADf data types
void functorDerivatives() {
  OstreamManager clout(std::cout,"functorDerivatives");

  // Consider a 1-dim analytical functor, defined with the underlying type S
  clout << "Consider functor f(x) = 2.1x + 1.5" << std::endl;
  AnalyticalLinear1D<S,S> f(2.1, 1.5);
  // Construct derivative of f
  AnalyticalDerivativeAD df(f);
  // Define arguments
  S s[1] = {1.3};
  S t[1];
  // Evaluate derivative
  df(t, s);
  // t now contains the derivative of f in s=1.3
  clout << "f'(1.3) = " << t[0] << '\n' << std::endl;

  // Sometimes, one wants to avoid contruction of the original functor f
  // (because it could be expensive). In that case, a different constructor of
  // the derivative functor can be used, where the original functor type and
  // dimensions have to be passed as template arguments (cf. the implementation
  // for the detailed list of arguments):
  clout << "Consider functor f(x) = 2.1x + 1.5 once more" << std::endl;
  AnalyticalDerivativeAD<AnalyticalLinear1D<S,S>,S,S,1,S,S> df2(1, 2.1, 1.5);
  // We evaluate it in the same way
  df2(t, s);
  clout << "f'(1.3) = " << t[0] << '\n' << std::endl;

  // All this works similarly for multidimensional functors
}

/** ---------------------------------------------------------------------------
 * Part 3: Evaluating derivatives of full simulations
 * ----------------------------------------------------------------------------
 * Sometimes, one is interested in the sensitivity of e.g. a drag coefficient
 * w.r.t. a design parameter. This can be done easily via wrapping the
 * simulation in a function which takes the design parameter as argument and
 * returns the drag coefficient. Then, one can apply the methods to derive
 * functions from part 1.
 */


int main(int argc, char **argv)
{
  initialize(&argc, &argv);

  functionADfApplication();
  functionDerivatives();
  functorADfApplication();
  functorDerivatives();

  return 0;
}

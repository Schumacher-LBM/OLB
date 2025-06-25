/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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

#ifndef ANALYTICAL_SOLUTION_TEST_FLOW_3D_H
#define ANALYTICAL_SOLUTION_TEST_FLOW_3D_H

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm

#include "functors/analytical/analyticalF.h"
#include "io/ostreamManager.h"

/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

////////////////////////////////////////////////////////////////////////////////
// 2nd level classes
// note: the number behind the AnalyticalF indicates the SOURCE dimension,
//       the target dimension for AnalyticalF is always 1 to reduce complexity
// note: for LatticeFunctions the number indicates the SOURCE dimension,
//       target dim depends on return variable type, so std::vector<T> is used

/// represents analytical functions without units, 1D-->1D

////////////////////////////////////////////////////////////////////////////////
// 3rd level classes: AnalyticalF



/// represents linear functions through 3 points (planes)
template <typename T, typename S, typename DESCRIPTOR>
class VelocityTestFlow3D : public AnalyticalF3D<T,S> {

protected:
  T u0;
public:
  VelocityTestFlow3D(UnitConverter<S,DESCRIPTOR> const& converter) : AnalyticalF3D<T,S>(3)
  {
    u0 = converter.getCharPhysVelocity();
    this->getName() = "velocityTestFlow3d";
  };

  bool operator()(T output[3], const S input[3])
  {

    T pi = 4.0 * util::atan(1.0);
    T twoPi = 2.0 * pi;
    T x = input[0];
    T y = input[1];
    T z = input[2];

    output[0] = (T) u0 * (util::sin((T) (twoPi * x)) * util::cos((T) (twoPi * z)) - util::cos((T) (twoPi * y)) * util::cos((T) (twoPi * x))) / (T) 4;
    output[1] = (T) u0 * (util::sin((T) (twoPi * y)) * util::sin((T) (twoPi * z)) + util::cos((T) (twoPi * x))) / (T) 4;
    output[2] = (T) -u0 * (util::cos((T) (twoPi * x)) * util::sin((T) (twoPi * z)) + 2 * util::cos((T) (twoPi * y)) * util::sin((T) (twoPi * x)) * (T) z * (T) pi - util::cos((T) (twoPi * y)) * util::cos((T) (twoPi * z))) / (T) 4;

    return true;
  };
};

template <typename T, typename S, typename DESCRIPTOR>
class PressureTestFlow3D : public AnalyticalF3D<T,S> {

protected:
  T charP;
public:
  PressureTestFlow3D(UnitConverter<S,DESCRIPTOR> const& converter) : AnalyticalF3D<T,S>(1)
  {
    charP = converter.getPhysDensity();
    this->getName() = "PressureTestFlow3d";
  };

  bool operator()(T output[3], const S input[3])
  {

    T pi = 4.0 * util::atan(1.0);
    T twoPi = 2.0 * pi;
    T x = input[0];
    T y = input[1];
    T z = input[2];

    output[0] = charP*util::cos((T) (twoPi * x)) * util::sin((T) (twoPi * y)) * z;

    return true;
  };
};

template <typename T, typename S, typename DESCRIPTOR>
class DensityTestFlow3D : public AnalyticalF3D<T,S> {

protected:
  PressureTestFlow3D<T,S,DESCRIPTOR>      pressure;
  UnitConverter<S,DESCRIPTOR> const&      converter;
public:
  DensityTestFlow3D(UnitConverter<S,DESCRIPTOR> const& converter_)
   : AnalyticalF3D<T,S>(1), pressure(converter_), converter(converter_)
  {
    this->getName() = "DensityTestFlow3d";
  };

  bool operator()(T output[3], const S input[3])
  {
    pressure(output, input);
    output[0] = converter.getPhysDensity(converter.getLatticeDensityFromPhysPressure(output[0]));

    return true;
  };
};

template <typename T, typename S, typename DESCRIPTOR>
class StrainRateTestFlow3D : public AnalyticalF3D<T,S> {

protected:
  T MF;
  T nu;
  T u0;
  T p0;
public:
  StrainRateTestFlow3D(UnitConverter<S,DESCRIPTOR> const& converter) : AnalyticalF3D<T,S>(9)
  {
    MF = converter.getConversionFactorForce() / converter.getConversionFactorMass();
    nu = converter.getPhysViscosity();
    u0 = converter.getCharPhysVelocity();
    p0 = converter.getPhysDensity();
    this->getName() = "ForceTestFlow3d";
  };

  bool operator()(T output[9], const S input[3])
  {

    T pi = 4.0 * util::atan(1.0);
    T twoPi = 2.0 * pi;
    T x = input[0];
    T y = input[1];
    T z = input[2];

    T DuDx = u0 * (twoPi*util::cos((T) (twoPi * x)) * util::cos((T) (twoPi * z)) + twoPi *util::cos((T) (twoPi * y)) * util::sin((T) (twoPi * x))) / (T) 4;
    T DuDy = (1./2.)*pi*u0*util::cos(2.*pi*x)*util::sin(2.*pi*y);
    T DuDz = (-(1./2.))*pi*u0*util::sin(2.*pi*x)*util::sin(2.*pi*z);
    T DvDx = (-(1./2.))*pi*u0*util::sin(2.*pi*x);
    T DvDy = u0 * twoPi* (util::cos((T) (twoPi * y)) * util::sin((T) (twoPi * z)) ) / (T) 4;
    T DvDz = (1./2.)*pi*u0*util::cos(2.*pi*z)*util::sin(2.*pi*y);
    T DwDx = (-(1./4.))*u0*(4.*pi*pi*z*util::cos(2.*pi*x)*util::cos(2.*pi*y) - 2.*pi*util::sin(2.*pi*x)*util::sin(2.*pi*z));
    T DwDy = (-(1./4.))*u0*(2.*pi*util::cos(2.*pi*z)*util::sin(2.*pi*y) - 4.*pi*pi*z*util::sin(2.*pi*x)*util::sin(2.*pi*y));
    T DwDz = (-(1./4.))*u0*(2.*pi*util::cos(2.*pi*x)*util::cos(2.*pi*z) + 2.*pi*util::cos(2.*pi*y)*util::sin(2.*pi*x) + 2.*pi*util::cos(2.*pi*y)*util::sin(2.*pi*z));

    output[0] = (DuDx + DuDx)/2.;
    output[1] = (DuDy + DvDx)/2.;
    output[2] = (DuDz + DwDx)/2.;
    output[3] = (DvDx + DuDy)/2.;
    output[4] = (DvDy + DvDy)/2.;
    output[5] = (DvDz + DwDy)/2.;
    output[6] = (DwDx + DuDz)/2.;
    output[7] = (DwDy + DvDz)/2.;
    output[8] = (DwDz + DwDz)/2.;

    return true;
  };
};


template <typename T, typename S, typename DESCRIPTOR>
class DissipationTestFlow3D : public AnalyticalF3D<T,S> {

protected:
  T MF;
  T nu;
  T u0;
  T p0;
public:
  DissipationTestFlow3D(UnitConverter<S,DESCRIPTOR> const& converter) : AnalyticalF3D<T,S>(1)
  {
    MF = converter.getConversionFactorForce() / converter.getConversionFactorMass();
    nu = converter.getPhysViscosity();
    u0 = converter.getCharPhysVelocity();
    p0 = converter.getPhysDensity();
    this->getName() = "ForceTestFlow3d";
  };

  bool operator()(T output[1], const S input[3])
  {

    T pi = 4.0 * util::atan(1.0);
    T twoPi = 2.0 * pi;
    T x = input[0];
    T y = input[1];
    T z = input[2];

    T DuDx = u0 * (twoPi*util::cos((T) (twoPi * x)) * util::cos((T) (twoPi * z)) + twoPi *util::cos((T) (twoPi * y)) * util::sin((T) (twoPi * x))) / (T) 4;
    T DuDy = (1./2.)*pi*u0*util::cos(2.*pi*x)*util::sin(2.*pi*y);
    T DuDz = (-(1./2.))*pi*u0*util::sin(2.*pi*x)*util::sin(2.*pi*z);
    T DvDx = (-(1./2.))*pi*u0*util::sin(2.*pi*x);
    T DvDy = u0 * twoPi* (util::cos((T) (twoPi * y)) * util::sin((T) (twoPi * z)) ) / (T) 4;
    T DvDz = (1./2.)*pi*u0*util::cos(2.*pi*z)*util::sin(2.*pi*y);
    T DwDx = (-(1./4.))*u0*(4.*pi*pi*z*util::cos(2.*pi*x)*util::cos(2.*pi*y) - 2.*pi*util::sin(2.*pi*x)*util::sin(2.*pi*z));
    T DwDy = (-(1./4.))*u0*(2.*pi*util::cos(2.*pi*z)*util::sin(2.*pi*y) - 4.*pi*pi*z*util::sin(2.*pi*x)*util::sin(2.*pi*y));
    T DwDz = (-(1./4.))*u0*(2.*pi*util::cos(2.*pi*x)*util::cos(2.*pi*z) + 2.*pi*util::cos(2.*pi*y)*util::sin(2.*pi*x) + 2.*pi*util::cos(2.*pi*y)*util::sin(2.*pi*z));

    output[0] = 0.5*nu*(util::pow(2*DuDx,2)+util::pow(2*DvDy,2)+util::pow(2*DwDz,2)+2.*util::pow(DuDy+DvDx,2)+2.*util::pow(DuDz+DwDx,2)+2.*util::pow(DvDz+DwDy,2));
    return true;
  };
};

template <typename T, typename S, typename DESCRIPTOR>
class ForceTestFlow3D : public AnalyticalF3D<T,S> {

protected:
  T rho;
  T nu;
  T u0;
  T p0;
public:
  ForceTestFlow3D(UnitConverter<S,DESCRIPTOR> const& converter) : AnalyticalF3D<T,S>(3)
  {
    rho = converter.getPhysDensity();
    nu = converter.getPhysViscosity();
    u0 = converter.getCharPhysVelocity();
    p0 = converter.getPhysDensity();
    this->getName() = "ForceTestFlow3d";
  };

  bool operator()(T output[3], const S input[3])
  {

    T pi = 4.0 * util::atan(1.0);
    T twoPi = 2.0 * pi;
    T x = input[0];
    T y = input[1];
    T z = input[2];

    output[0] = (T) -(T) pi / rho * (-(T) 16 * nu * u0 * (T) pi * rho * util::sin((T) (twoPi * x)) * util::cos((T) (twoPi * z))
                                    + (T) 16 * nu * u0 * (T) pi * rho * util::cos((T) (twoPi * y)) * util::cos((T) (twoPi * x))
                                    - u0 * u0 * rho * util::cos((T) (twoPi * z)) * util::cos((T) (twoPi * y))
                                    + 2 * u0 * u0 * rho * util::cos((T) (twoPi * y)) * util::pow(util::cos((T) (twoPi * x)), 2) * util::cos((T) (twoPi * z))
                                    + u0 * u0 * rho * util::pow(util::cos((T) (twoPi * y)), 2) * util::cos((T) (twoPi * x)) * util::sin((T) (twoPi * x))
                                    - u0 * u0 * util::cos((T) (twoPi * x)) * rho * util::sin((T) (twoPi * z))
                                    + u0 * u0 * util::cos((T) (twoPi * x)) * rho * util::sin((T) (twoPi * z)) * util::pow(util::cos((T) (twoPi * y)), 2)
                                    - u0 * u0 * util::sin((T) (twoPi * y)) * util::pow(util::cos((T) (twoPi * x)), 2) * rho
                                    - u0 * u0 * rho * util::sin((T) (twoPi * x)) * util::cos((T) (twoPi * x))
                                    - 2 * u0 * u0 * util::sin((T) (twoPi * z)) * rho * util::cos((T) (twoPi * y)) * (T) z * (T) pi
                                    + 2 * u0 * u0 * util::sin((T) (twoPi * z)) * rho * util::cos((T) (twoPi * y)) * (T) z * (T) pi * util::pow(util::cos((T) (twoPi * x)), 2)
                                    + u0 * u0 * util::sin((T) (twoPi * x)) * util::sin((T) (twoPi * z)) * rho * util::cos((T) (twoPi * y)) * util::cos((T) (twoPi * z))
                                    + (T) 16 * p0 * util::sin((T) (twoPi * x)) * util::sin((T) (twoPi * y)) * (T) z) / (T) 8;

    output[1] = (T) (T) pi / rho * ((T) 8 * nu * u0 * (T) pi * rho * util::cos((T) (twoPi * x))
                                    + (T) 16 * nu * u0 * (T) pi * rho * util::sin((T) (twoPi * y)) * util::sin((T) (twoPi * z))
                                    - u0 * u0 * rho * util::cos((T) (twoPi * z)) + u0 * u0 * rho * util::cos((T) (twoPi * z)) * util::pow(util::cos((T) (twoPi * x)), 2)
                                    + u0 * u0 * util::sin((T) (twoPi * x)) * rho * util::cos((T) (twoPi * y)) * util::cos((T) (twoPi * x))
                                    + u0 * u0 * util::cos((T) (twoPi * y)) * rho * util::sin((T) (twoPi * y))
                                    + u0 * u0 * util::cos((T) (twoPi * y)) * util::sin((T) (twoPi * z)) * rho * util::cos((T) (twoPi * x))
                                    - u0 * u0 * util::sin((T) (twoPi * y)) * util::cos((T) (twoPi * z)) * rho * util::cos((T) (twoPi * x)) * util::sin((T) (twoPi * z))
                                    - 2 * u0 * u0 * util::sin((T) (twoPi * y)) * util::cos((T) (twoPi * z)) * rho * util::cos((T) (twoPi * y)) * util::sin((T) (twoPi * x)) * (T) z * (T) pi
                                    + (T) 16 * p0 * util::cos((T) (twoPi * x)) * util::cos((T) (twoPi * y)) * (T) z) / (T) 8;

    output[2] = (T) -((T) 16 * nu * u0 * (T) (pi * pi) * rho * util::cos((T) (twoPi * x)) * util::sin((T) (twoPi * z))
                      + (T) 32 * nu * u0 * (T) util::pow((T) pi, (T) 3) * rho * util::cos((T) (twoPi * y)) * util::sin((T) (twoPi * x)) * (T) z
                      - (T) 16 * nu * u0 * (T) (pi * pi) * rho * util::cos((T) (twoPi * y)) * util::cos((T) (twoPi * z))
                      - 2 * (T) (pi * pi) * u0 * u0 * rho * util::sin((T) (twoPi * z)) * util::sin((T) (twoPi * x)) * (T) z
                      - 2 * (T) (pi * pi) * u0 * u0 * util::sin((T) (twoPi * y)) * rho * util::cos((T) (twoPi * x)) * util::sin((T) (twoPi * x)) * (T) z
                      + (T) pi * u0 * u0 * util::sin((T) (twoPi * y)) * rho * util::cos((T) (twoPi * x)) * util::cos((T) (twoPi * z))
                      - u0 * u0 * (T) pi * rho * util::cos((T) (twoPi * x)) * util::cos((T) (twoPi * y))
                      + 2 * u0 * u0 * (T) pi * rho * util::cos((T) (twoPi * y)) * util::pow(util::cos((T) (twoPi * z)), 2) * util::cos((T) (twoPi * x))
                      - 2 * u0 * u0 * (T) (pi * pi) * rho * util::pow(util::cos((T) (twoPi * y)), 2) * (T) z
                      + u0 * u0 * (T) pi * rho * util::pow(util::cos((T) (twoPi * y)), 2) * util::cos((T) (twoPi * z)) * util::sin((T) (twoPi * x))
                      - (T) 8 * p0 * util::cos((T) (twoPi * x)) * util::sin((T) (twoPi * y))) / rho / (T) 8;

    return true;
  };
};

} // end namespace olb

#endif

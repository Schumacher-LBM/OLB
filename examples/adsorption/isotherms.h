/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2022 Florian Raichle
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

#ifndef OLB_APPS_FLORIAN_ADSORPTION3D_PHOSPHATEREACTION_H_
#define OLB_APPS_FLORIAN_ADSORPTION3D_PHOSPHATEREACTION_H_

#include <vector>
#include <cmath>
#include <iostream>

namespace olb {

template<typename T>
class LinearIsotherm: public Isotherm<T> {
 public:
  LinearIsotherm(T isoKonstA, T isoKonstB): Isotherm<T>(isoKonstA, isoKonstB){}

  T getLoading(T c) const override {
    return this->isoKonstA * c;
  }
};

template<typename T>
class LangmuirIsotherm: public Isotherm<T> {
 public:
  LangmuirIsotherm(T isoKonstA, T isoKonstB): Isotherm<T>(isoKonstA, isoKonstB){}

  T getLoading(T c) const override {
    return this->isoKonstA*this->isoKonstB*c/(1+this->isoKonstB*c);
  }
};

template<typename T>
class FreundlichIsotherm: public Isotherm<T> {
 public:
  FreundlichIsotherm(T isoKonstA, T isoKonstB): Isotherm<T>(isoKonstA, isoKonstB){}

  T getLoading(T c) const override {
    return this->isoKonstA * std::pow(c, this->isoKonstB);
  }
};

}
#endif

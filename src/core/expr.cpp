/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#include <iostream>

#include "expr.h"

#include "utilities/oalgorithm.h"

namespace olb {

std::string Expr::Symbol::describe() const {
  return name;
}
std::size_t Expr::Symbol::size() const {
  return 0;
}

Expr::Constant::Constant():
  value{} { }

Expr::Constant::Constant(double x):
  value{x} { }

std::string Expr::Constant::describe() const {
  return std::to_string(value);
}

std::size_t Expr::Constant::size() const {
  return 0;
}

Expr::Binary::Binary(Expr lhs, Op op, Expr rhs):
  _lhs{std::make_shared<Expr>(lhs)},
  op{op},
  _rhs{std::make_shared<Expr>(rhs)} { }

std::string Expr::Binary::describe() const {
  switch (op) {
  case Op::Add: return "(" + _lhs->describe() + ") + (" + _rhs->describe() + ")";
  case Op::Sub: return "(" + _lhs->describe() + ") - (" + _rhs->describe() + ")";
  case Op::Mul: return "(" + _lhs->describe() + ") * (" + _rhs->describe() + ")";
  case Op::Div: return "(" + _lhs->describe() + ") / (" + _rhs->describe() + ")";
  case Op::Pow: return "Pow(" + _lhs->describe() + "," + _rhs->describe() + ")";
  default: throw std::invalid_argument("Unsupported binary operation");
  }
}

std::size_t Expr::Binary::size() const {
  return 1 + _lhs->size() + _rhs->size();
}

Expr::Unary::Unary(Op op, Expr arg):
  op{op},
  _arg{std::make_shared<Expr>(arg)} { }

std::string Expr::Unary::describe() const {
  switch (op) {
  case Op::Sqrt: return "sqrt(" + _arg->describe() + ")";
  case Op::Abs:  return "Abs(" + _arg->describe() + ")";
  case Op::Exp:  return "exp(" + _arg->describe() + ")";
  default: throw std::invalid_argument("Unsupported unary operation");
  }
}

std::size_t Expr::Unary::size() const {
  return 1 + _arg->size();
}

Expr::Expr():
  _payload(Constant{}) { }

Expr::Expr(Expr&& rhs):
  _payload(rhs._payload) { }

Expr::Expr(const Expr& rhs):
  _payload(rhs._payload) { }

Expr::Expr(double v):
  _payload(Constant{v}) { }

Expr::Expr(std::string name):
  _payload(Symbol{name}) { }

Expr::Expr(Expr lhs, Op op, Expr rhs):
  _payload(Binary(lhs, op, rhs)) { }

Expr::Expr(Op op, Expr rhs):
  _payload(Unary(op, rhs)) { }

std::string Expr::describe() const {
  std::string out;
  std::visit([&out](auto& x) {
    out = x.describe();
  }, _payload);
  return out;
};

std::size_t Expr::size() const {
  std::size_t size = 0;
  std::visit([&size](auto& x) {
    size = x.size();
  }, _payload);
  return size;
}

Expr& Expr::operator=(const Expr& rhs) {
  _payload = rhs._payload;

  return *this;
}

Expr& Expr::operator+=(Expr rhs) {
  _payload = Binary(*this, Op::Add, rhs);
  return *this;
}

Expr& Expr::operator-=(Expr rhs) {
  _payload = Binary(*this, Op::Sub, rhs);
  return *this;
}

Expr& Expr::operator*=(Expr rhs) {
  _payload = Binary(*this, Op::Mul, rhs);
  return *this;
}

Expr& Expr::operator/=(Expr rhs) {
  _payload = Binary(*this, Op::Div, rhs);
  return *this;
}

Expr operator+(Expr lhs, Expr rhs) {
  return Expr(lhs, Expr::Op::Add, rhs);
}

Expr operator-(Expr lhs, Expr rhs) {
  return Expr(lhs, Expr::Op::Sub, rhs);
}

Expr operator*(Expr lhs, Expr rhs) {
  return Expr(lhs, Expr::Op::Mul, rhs);
}

Expr operator/(Expr lhs, Expr rhs) {
  return Expr(lhs, Expr::Op::Div, rhs);
}

Expr operator-(Expr rhs) {
  return Expr(-1.) * rhs;
}

bool operator==(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '=='");
}

bool operator!=(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '!='");
}

bool operator>(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '>'");
}

bool operator<(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '<'");
}

bool operator>=(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '>='");
}

bool operator<=(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '<='");
}

namespace util {

Expr sqrt(Expr x) {
  return Expr(Expr::Op::Sqrt, x);
}

Expr fabs(Expr x) {
  return Expr(Expr::Op::Abs, x);
}

Expr pow(Expr base, Expr exp) {
  return Expr(base, Expr::Op::Pow, exp);
}

Expr exp(Expr x) {
  return Expr(Expr::Op::Exp, x);
}

Expr max(Expr a, Expr b) {
  throw std::domain_error("Invalid operator for Expr type: 'max'");
}

Expr min(Expr a, Expr b) {
  throw std::domain_error("Invalid operator for Expr type: 'min'");
}

}

}

/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <cmath>
#include <memory>

#include "ufo/ObsBiasIncrement.h"

#include "oops/util/Logger.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const eckit::Configuration & conf)
  : biasinc_() {
  std::unique_ptr<ObsBiasBase> biasbase(ObsBiasFactory::create(conf));
  if (biasbase) {
    for (std::size_t ii = 0; ii < biasbase->size(); ++ii)
      biasinc_.push_back(0.0);
  }
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other, const bool copy)
  : biasinc_(other.size(), 0.0) {
  if (copy) {
    for (std::size_t ii = 0; ii < other.size(); ++ii)
      biasinc_[ii] = other[ii];
  }
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other,
                                   const eckit::Configuration & conf)
  : biasinc_(other.size(), 0.0) {
  for (std::size_t ii = 0; ii < other.size(); ++ii)
    biasinc_[ii] = other[ii];
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::diff(const ObsBias & b1, const ObsBias & b2) {
  for (std::size_t ii = 0; ii < this->size(); ++ii) {
    biasinc_[ii] = b1[ii] - b2[ii];
  }
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::zero() {
  for (std::size_t ii = 0; ii < this->size(); ++ii) biasinc_[ii] = 0.0;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator=(const ObsBiasIncrement & rhs) {
  for (std::size_t ii = 0; ii < this->size(); ++ii)
    biasinc_[ii] = rhs[ii];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator+=(const ObsBiasIncrement & rhs) {
  for (std::size_t ii = 0; ii < this->size(); ++ii)
    biasinc_[ii] += rhs[ii];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator-=(const ObsBiasIncrement & rhs) {
  for (std::size_t ii = 0; ii < this->size(); ++ii)
    biasinc_[ii] -= rhs[ii];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator*=(const double fact) {
  for (std::size_t ii = 0; ii < this->size(); ++ii)
    biasinc_[ii] *= fact;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::axpy(const double fact, const ObsBiasIncrement & rhs) {
  for (std::size_t ii = 0; ii < this->size(); ++ii)
    biasinc_[ii] += fact * rhs[ii];
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::dot_product_with(const ObsBiasIncrement & rhs) const {
  double zz = 0.0;
  for (std::size_t ii = 0; ii < this->size(); ++ii) {
    zz += biasinc_[ii] * rhs[ii];
  }
  return zz;
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::norm() const {
  double zz = 0.0;
  for (std::size_t ii = 0; ii < this->size(); ++ii) {
    zz += biasinc_[ii]*biasinc_[ii];
  }
  if (this->size() > 0) zz = std::sqrt(zz/this->size());
  return zz;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

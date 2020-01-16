/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <cmath>
#include <memory>

#include "ufo/ObsBiasIncrement.h"

#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ioda::ObsSpace & obs, const eckit::Configuration & conf)
  : biasbase_(LinearObsBiasFactory::create(obs, conf)), conf_(conf) {
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other, const bool copy)
  : biasbase_(), conf_(other.config()) {
  if (other) {
    biasbase_.reset(LinearObsBiasFactory::create(other.obspace(), other.config()));
    if (copy) *biasbase_ = other;
  }
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other,
                                   const eckit::Configuration & conf)
  : biasbase_(LinearObsBiasFactory::create(other.obspace(), conf)), conf_(conf) {
  /*
   * As we don't know the details now, it needs revisit later
   */
  if (biasbase_) *biasbase_ = other;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::diff(const ObsBias & b1, const ObsBias & b2) {
  if (biasbase_) biasbase_->diff(b1, b2);
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::zero() {
  if (biasbase_) biasbase_->zero();
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator=(const ObsBiasIncrement & rhs) {
  if (biasbase_) *biasbase_ = rhs;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator+=(const ObsBiasIncrement & rhs) {
  if (biasbase_) *biasbase_ += rhs;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator-=(const ObsBiasIncrement & rhs) {
  if (biasbase_) *biasbase_ -= rhs;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator*=(const double fact) {
  if (biasbase_) *biasbase_ *= fact;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::axpy(const double fact, const ObsBiasIncrement & rhs) {
  if (biasbase_) biasbase_->axpy(fact, rhs);
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::dot_product_with(const ObsBiasIncrement & rhs) const {
  double zz = 0.0;
  if (biasbase_) zz = biasbase_->dot_product_with(rhs);
  return zz;
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::norm() const {
  double zz = 0.0;
  if (biasbase_) zz = biasbase_->norm();
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::computeObsBiasTL(const GeoVaLs & geovals,
                                        const ioda::ObsDataVector<float> & preds,
                                        ioda::ObsVector & ybiasinc) const {
  if (biasbase_) biasbase_->computeObsBiasTL(geovals, preds, ybiasinc);
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::computeObsBiasAD(GeoVaLs & geovals,
                                        const ioda::ObsDataVector<float> & preds,
                                        const ioda::ObsVector & ybiasinc) {
  if (biasbase_) biasbase_->computeObsBiasAD(geovals, preds, ybiasinc);
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::print(std::ostream & os) const {
  if (biasbase_) os << *biasbase_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/ObsBias.h"

#include "oops/util/Logger.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const eckit::Configuration & conf)
  : biasbase_(ObsBiasFactory::create(conf)), conf_(conf), vars_() {
  if (biasbase_) vars_ = biasbase_->variables();
}

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : biasbase_(ObsBiasFactory::create(other.config())), conf_(other.config()),
    vars_(other.vars_) {
  if (copy && biasbase_) {
    for (std::size_t jj =0; jj < other.size(); ++jj)
      (*biasbase_)[jj] = other[jj];
  }
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator+=(const ObsBiasIncrement & dx) {
  if (biasbase_) *biasbase_+=dx;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBias::read(const eckit::Configuration & conf) {
  if (biasbase_) biasbase_->read(conf);
}

// -----------------------------------------------------------------------------

void ObsBias::write(const eckit::Configuration & conf) const {
  if (biasbase_) biasbase_->write(conf);
}

// -----------------------------------------------------------------------------

void ObsBias::computeObsBias(const GeoVaLs & geovals,
                             ioda::ObsVector & ybias,
                             const ioda::ObsSpace & os) const {
  if (biasbase_) biasbase_->computeObsBias(geovals, ybias, os);
}

// -----------------------------------------------------------------------------

double ObsBias::norm() const {
  double zz = 0.0;
  if (biasbase_) zz = biasbase_->norm();
  return zz;
}

// -----------------------------------------------------------------------------

std::size_t ObsBias::size() const {
  std::size_t zz = 0;
  if (biasbase_) zz = biasbase_->size();
  return zz;
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsBias::variables() const {
  return vars_;
}

// -----------------------------------------------------------------------------

void ObsBias::print(std::ostream & os) const {
  if (biasbase_) os << *biasbase_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

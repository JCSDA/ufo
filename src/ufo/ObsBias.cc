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
  : biasbase_(ObsBiasFactory::create(conf)), conf_(conf), geovars_(), hdiags_() {
  if (biasbase_) {
    geovars_ += biasbase_->requiredGeoVaLs();
    hdiags_  += biasbase_->requiredHdiagnostics();
  }
}

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : biasbase_(ObsBiasFactory::create(other.config())), conf_(other.config()),
    geovars_(), hdiags_() {
  if (copy && biasbase_) *biasbase_ = other;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator+=(const ObsBiasIncrement & dx) {
  if (biasbase_) *biasbase_+=dx;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator=(const ObsBias & rhs) {
/* temporarily comment out for coverage check
  if (biasbase_) {
    *biasbase_ = rhs;
    geovars_ += biasbase_->requiredGeoVaLs();
    hdiags_  += biasbase_->requiredHdiagnostics();
  }
*/
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
                             const ioda::ObsSpace & os,
                             const ObsDiagnostics & ydiags) const {
  if (biasbase_) biasbase_->computeObsBias(geovals, ybias, os, ydiags);
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

void ObsBias::print(std::ostream & os) const {
  if (biasbase_) os << *biasbase_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

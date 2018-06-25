/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsOperator.h"

#include "eckit/config/Configuration.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsOperatorBase.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsOperator::ObsOperator(const ioda::ObsSpace & os, const eckit::Configuration & conf)
  : oper_(ObsOperatorFactory::create(os, conf))
{}

// -----------------------------------------------------------------------------

ObsOperator::~ObsOperator() {}

// -----------------------------------------------------------------------------

void ObsOperator::simulateObs(const GeoVaLs & gvals, ioda::ObsVector & yy, const ObsBias & bias) const {
  oper_->simulateObs(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsOperator::variables() const {
  return oper_->variables();
}

// -----------------------------------------------------------------------------

void ObsOperator::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

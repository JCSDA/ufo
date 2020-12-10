/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "ufo/ObsOperator.h"

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsOperator::ObsOperator(ioda::ObsSpace & os, const eckit::Configuration & conf)
  : oper_(ObsOperatorFactory::create(os, conf)), odb_(os)
{}

// -----------------------------------------------------------------------------

ObsOperator::~ObsOperator() {}

// -----------------------------------------------------------------------------

void ObsOperator::simulateObs(const GeoVaLs & gvals, ioda::ObsVector & yy,
                              const ObsBias & bias, ObsDiagnostics & ydiags) const {
  oper_->simulateObs(gvals, yy, ydiags);
  if (bias) {
    ioda::ObsVector ybias(odb_);
    bias.computeObsBias(ybias, ydiags, bias.computePredictors(gvals, ydiags));
    ybias.save("ObsBias");
  }
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsOperator::requiredVars() const {
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

std::unique_ptr<Locations> ObsOperator::locations() const {
  return oper_->locations();
}

// -----------------------------------------------------------------------------

void ObsOperator::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

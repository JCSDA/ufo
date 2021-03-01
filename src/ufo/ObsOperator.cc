/*
 * (C) Copyright 2021 UCAR
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
#include "ufo/ObsBiasOperator.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsOperator::ObsOperator(ioda::ObsSpace & os, const eckit::Configuration & conf)
  : oper_(ObsOperatorFactory::create(os, conf)), odb_(os)
{
  // We use += rather than = to make sure the Variables objects contain no duplicate entries
  // and the variables are sorted alphabetically.
  oops::Variables operatorVars;
  operatorVars += oper_->simulatedVars();
  oops::Variables obsSpaceVars;
  obsSpaceVars += os.obsvariables();
  if (!(operatorVars == obsSpaceVars))
    throw eckit::UserError("The list of variables simulated by the obs operator differs from "
                           "the list of simulated variables in the obs space",
                           Here());
}

// -----------------------------------------------------------------------------

void ObsOperator::simulateObs(const GeoVaLs & gvals, ioda::ObsVector & yy,
                              const ObsBias & bias, ObsDiagnostics & ydiags) const {
  oper_->simulateObs(gvals, yy, ydiags);
  if (bias) {
    ioda::ObsVector ybias(odb_);
    ObsBiasOperator biasoper(odb_);
    biasoper.computeObsBias(gvals, ybias, bias, ydiags);
    // update H(x) with bias correction
    yy += ybias;
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

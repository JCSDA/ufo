/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsOperator.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Locations.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasOperator.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsOperator::ObsOperator(ioda::ObsSpace & os, const Parameters_ & params)
  : oper_(ObsOperatorFactory::create(os, params.operatorParameters)), odb_(os)
{
  // We use += rather than = to make sure the Variables objects contain no duplicate entries.
  oops::Variables operatorVars;
  operatorVars += oper_->simulatedVars();
  operatorVars.sort();
  oops::Variables obsSpaceVars;
  obsSpaceVars += os.assimvariables();
  obsSpaceVars.sort();
  if (!(operatorVars == obsSpaceVars))
    throw eckit::UserError("The list of variables simulated by the obs operator differs from "
                           "the list of simulated variables in the obs space",
                           Here());
}

// -----------------------------------------------------------------------------

void ObsOperator::simulateObs(const GeoVaLs & gvals, ioda::ObsVector & yy,
                              const ObsBias & biascoeff, ioda::ObsVector & ybias,
                              ObsDiagnostics & ydiags) const {
  oper_->simulateObs(gvals, yy, ydiags);
  if (biascoeff) {
    ObsBiasOperator biasoper(odb_);
    biasoper.computeObsBias(gvals, ybias, biascoeff, ydiags);
    // update H(x) with bias correction
    yy += ybias;
  }
}

// -----------------------------------------------------------------------------

const oops::Variables & ObsOperator::requiredVars() const {
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

ObsOperator::Locations_ ObsOperator::locations() const {
  return oper_->locations();
}

// -----------------------------------------------------------------------------

void ObsOperator::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

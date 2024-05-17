/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsOperator.h"

#include "eckit/config/Configuration.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Locations.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasOperator.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"
#include "ufo/ScopedDefaultGeoVaLFormatChange.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsOperator::ObsOperator(ioda::ObsSpace & os, const eckit::Configuration & config)
  : oper_(), odb_(os)
{
  Parameters_ params;
  params.validateAndDeserialize(config);
  oper_.reset(ObsOperatorFactory::create(os, params.operatorParameters));
  // We use += rather than = to make sure the Variables objects contain no duplicate entries.
  oops::ObsVariables operatorVars;
  operatorVars += oper_->simulatedVars();
  operatorVars.sort();
  oops::ObsVariables obsSpaceVars;
  obsSpaceVars += os.assimvariables();
  obsSpaceVars.sort();
  if (!(operatorVars == obsSpaceVars))
    throw eckit::UserError("The list of variables simulated by the obs operator differs from "
                           "the list of simulated variables in the obs space",
                           Here());
}

// -----------------------------------------------------------------------------
void ObsOperator::simulateObs(const GeoVaLs & gvals, ioda::ObsVector & yy,
                              const ObsBias & biascoeff, const QCFlags_t & qc_flags,
                              ioda::ObsVector & ybias, ObsDiagnostics & ydiags) const
{
  ScopedDefaultGeoVaLFormatChange change(gvals, GeoVaLFormat::SAMPLED);
  oper_->simulateObs(gvals, yy, ydiags, qc_flags);
  if (biascoeff) {
    ObsBiasOperator biasoper(odb_);
    biasoper.computeObsBias(gvals, ybias, biascoeff, ydiags, qc_flags);
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

void ObsOperator::computeReducedVars(const oops::Variables & vars, GeoVaLs & geovals) const {
  return oper_->computeReducedVars(vars, geovals);
}

// -----------------------------------------------------------------------------

void ObsOperator::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

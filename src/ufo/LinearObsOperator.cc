/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/LinearObsOperator.h"

#include <utility>
#include <vector>

#include "ioda/ObsVector.h"
#include "oops/base/Locations.h"
#include "oops/interface/SampledLocations.h"
#include "ufo/LinearObsBiasOperator.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"

namespace ufo {

// -----------------------------------------------------------------------------

LinearObsOperator::LinearObsOperator(ioda::ObsSpace & os, const Parameters_ & params)
  : oper_(LinearObsOperatorFactory::create(os, params.operatorParameters)), odb_(os)
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

void LinearObsOperator::setTrajectory(const GeoVaLs & gvals, const ObsBias & bias) {
  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  oops::Variables vars;
  vars += bias.requiredHdiagnostics();
  std::vector<float> lons(odb_.nlocs());
  std::vector<float> lats(odb_.nlocs());
  std::vector<util::DateTime> times(odb_.nlocs());
  odb_.get_db("MetaData", "latitude", lats);
  odb_.get_db("MetaData", "longitude", lons);
  odb_.get_db("MetaData", "dateTime", times);
  auto locs = std::make_unique<SampledLocations>(lons, lats, times, odb_.distribution());
  ObsDiagnostics ydiags(odb_, SampledLocations_(std::move(locs)), vars);
  oper_->setTrajectory(gvals, ydiags);
  if (bias) {
    biasoper_.reset(new LinearObsBiasOperator(odb_));
    biasoper_->setTrajectory(gvals, bias, ydiags);
  }
}

// -----------------------------------------------------------------------------

void LinearObsOperator::simulateObsTL(const GeoVaLs & gvals, ioda::ObsVector & yy,
                                      const ObsBiasIncrement & bias) const {
  oper_->simulateObsTL(gvals, yy);
  if (bias) {
    ioda::ObsVector ybiasinc(odb_);
    biasoper_->computeObsBiasTL(bias, ybiasinc);
    yy += ybiasinc;
  }
}

// -----------------------------------------------------------------------------

void LinearObsOperator::simulateObsAD(GeoVaLs & gvals, const ioda::ObsVector & yy,
                                      ObsBiasIncrement & bias) const {
  oper_->simulateObsAD(gvals, yy);
  if (bias) {
    ioda::ObsVector ybiasinc(yy);
    biasoper_->computeObsBiasAD(bias, ybiasinc);
  }
}

// -----------------------------------------------------------------------------

const oops::Variables & LinearObsOperator::requiredVars() const {
  return oper_->requiredVars();
}

// -----------------------------------------------------------------------------

void LinearObsOperator::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

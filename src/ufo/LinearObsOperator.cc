/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/LinearObsOperator.h"

#include <vector>

#include "ioda/ObsVector.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

LinearObsOperator::LinearObsOperator(ioda::ObsSpace & os, const eckit::Configuration & conf)
  : oper_(LinearObsOperatorFactory::create(os, conf)), odb_(os)
{}

// -----------------------------------------------------------------------------

void LinearObsOperator::setTrajectory(const GeoVaLs & gvals, const ObsBias & bias) {
  oops::Variables vars;
  vars += bias.requiredHdiagnostics();
  std::vector<float> lons(odb_.nlocs());
  std::vector<float> lats(odb_.nlocs());
  std::vector<util::DateTime> times(odb_.nlocs());
  odb_.get_db("MetaData", "latitude", lats);
  odb_.get_db("MetaData", "longitude", lons);
  odb_.get_db("MetaData", "datetime", times);
  ObsDiagnostics ydiags(odb_, Locations(lons, lats, times, odb_.comm()), vars);
  oper_->setTrajectory(gvals, bias, ydiags);
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
    biasoper_->computeObsBiasTL(gvals, bias, ybiasinc);
    yy += ybiasinc;
  }
}

// -----------------------------------------------------------------------------

void LinearObsOperator::simulateObsAD(GeoVaLs & gvals, const ioda::ObsVector & yy,
                                      ObsBiasIncrement & bias) const {
  oper_->simulateObsAD(gvals, yy);
  if (bias) {
    ioda::ObsVector ybiasinc(yy);
    biasoper_->computeObsBiasAD(gvals, bias, ybiasinc);
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

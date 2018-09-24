/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/LinearObsOperator.h"

#include "ufo/LinearObsOperatorBase.h"

namespace ufo {

// -----------------------------------------------------------------------------

LinearObsOperator::LinearObsOperator(const ioda::ObsSpace & os, const eckit::Configuration & conf)
  : oper_(LinearObsOperatorFactory::create(os, conf))
{}

// -----------------------------------------------------------------------------

LinearObsOperator::~LinearObsOperator() {}

// -----------------------------------------------------------------------------

void LinearObsOperator::setTrajectory(const GeoVaLs & gvals, const ObsBias & bias) {
  oper_->setTrajectory(gvals, bias);
}

// -----------------------------------------------------------------------------

void LinearObsOperator::simulateObsTL(const GeoVaLs & gvals, ioda::ObsVector & yy,
                                      const ObsBiasIncrement & bias) const {
  oper_->simulateObsTL(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

void LinearObsOperator::simulateObsAD(GeoVaLs & gvals, const ioda::ObsVector & yy,
                                      ObsBiasIncrement & bias) const {
  oper_->simulateObsAD(gvals, yy, bias);
}

// -----------------------------------------------------------------------------

const oops::Variables & LinearObsOperator::variables() const {
  return oper_->variables();
}

// -----------------------------------------------------------------------------

void LinearObsOperator::print(std::ostream & os) const {
  os << *oper_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

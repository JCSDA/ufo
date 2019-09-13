/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ioda/ObsVector.h"
#include "ufo/LinearObsOperator.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------

LinearObsOperator::LinearObsOperator(ioda::ObsSpace & os, const eckit::Configuration & conf)
  : oper_(LinearObsOperatorFactory::create(os, conf)), odb_(os)
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
  oper_->simulateObsTL(gvals, yy);
  ioda::ObsVector ybiasinc(odb_);
  bias.computeObsBiasTL(gvals, ybiasinc, odb_);
  yy += ybiasinc;
}

// -----------------------------------------------------------------------------

void LinearObsOperator::simulateObsAD(GeoVaLs & gvals, const ioda::ObsVector & yy,
                                      ObsBiasIncrement & bias) const {
  ioda::ObsVector ybiasinc(yy);
  oper_->simulateObsAD(gvals, yy);
  bias.computeObsBiasAD(gvals, ybiasinc, odb_);
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

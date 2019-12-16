/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ioda/ObsVector.h"
#include "ufo/LinearObsOperator.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

LinearObsOperator::LinearObsOperator(ioda::ObsSpace & os, const eckit::Configuration & conf)
  : oper_(LinearObsOperatorFactory::create(os, conf)), odb_(os), biaspreds_()
{}

// -----------------------------------------------------------------------------

LinearObsOperator::~LinearObsOperator() {}

// -----------------------------------------------------------------------------

void LinearObsOperator::setTrajectory(const GeoVaLs & gvals, const ObsBias & bias) {
  oops::Variables vars;
  if (bias) vars += bias.requiredHdiagnostics();
  ObsDiagnostics ydiags(odb_, Locations(odb_, odb_.windowStart(), odb_.windowEnd()), vars);
  oper_->setTrajectory(gvals, bias, ydiags);
  if (bias) bias.computeObsBiasPredictors(gvals, ydiags, biaspreds_);
}

// -----------------------------------------------------------------------------

void LinearObsOperator::simulateObsTL(const GeoVaLs & gvals, ioda::ObsVector & yy,
                                      const ObsBiasIncrement & bias) const {
  oper_->simulateObsTL(gvals, yy);
  if (bias) {
    ioda::ObsVector ybiasinc(odb_);
    bias.computeObsBiasTL(gvals, *biaspreds_, ybiasinc);
    yy += ybiasinc;
  }
}

// -----------------------------------------------------------------------------

void LinearObsOperator::simulateObsAD(GeoVaLs & gvals, const ioda::ObsVector & yy,
                                      ObsBiasIncrement & bias) const {
  oper_->simulateObsAD(gvals, yy);
  if (bias) {
    ioda::ObsVector ybiasinc(yy);
    bias.computeObsBiasAD(gvals, *biaspreds_, ybiasinc);
  }
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

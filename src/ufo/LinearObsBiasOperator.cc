/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/LinearObsBiasOperator.h"

#include <memory>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------

LinearObsBiasOperator::LinearObsBiasOperator(ioda::ObsSpace & odb)
  : odb_(odb) {
  oops::Log::trace() << "LinearObsBiasOperator::create done." << std::endl;
}

// -----------------------------------------------------------------------------

void LinearObsBiasOperator::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias,
                                          ObsDiagnostics & ydiags) {
  oops::Log::trace() << "LinearObsBiasOperator::setTrajectory starts." << std::endl;
  const std::vector<std::shared_ptr<const PredictorBase>> variablePredictors =
      bias.variablePredictors();
  const std::size_t npreds = variablePredictors.size();

  predData_.resize(npreds, ioda::ObsVector(odb_));
  for (std::size_t p = 0; p < npreds; ++p) {
    variablePredictors[p]->compute(odb_, geovals, ydiags, bias, predData_[p]);
  }

  oops::Log::trace() << "LinearObsBiasOperator::setTrajectory done." << std::endl;
}

// -----------------------------------------------------------------------------

void LinearObsBiasOperator::computeObsBiasTL(const GeoVaLs & geovals,
                                             const ObsBiasIncrement & biascoeffinc,
                                             ioda::ObsVector & ybiasinc) const {
  oops::Log::trace() << "LinearObsBiasOperator::computeObsBiasTL starts." << std::endl;

  const size_t npreds = predData_.size();

  ybiasinc.zero();
  for (size_t jpred = 0; jpred < npreds; ++jpred) {
    ybiasinc.axpy(biascoeffinc.coefficients(jpred), predData_[jpred]);
  }

  oops::Log::trace() << "LinearObsBiasOperator::computeObsBiasTL done." << std::endl;
}

// -----------------------------------------------------------------------------

void LinearObsBiasOperator::computeObsBiasAD(GeoVaLs & geovals,
                                             ObsBiasIncrement & biascoeffinc,
                                             const ioda::ObsVector & ybiasinc) const {
  oops::Log::trace() << "LinearObsBiasOperator::computeObsBiasAD starts." << std::endl;

  const size_t npreds = predData_.size();

  for (std::size_t jpred = 0; jpred < npreds; ++jpred) {
    biascoeffinc.updateCoeff(jpred, predData_[jpred].multivar_dot_product_with(ybiasinc));
  }

  oops::Log::trace() << "LinearObsBiasOperator::computeAD done." << std::endl;
}

// -----------------------------------------------------------------------------

void LinearObsBiasOperator::print(std::ostream & os) const {
  os << "TL/AD bias computed as linear combination of predictors.";
}

// -----------------------------------------------------------------------------

}  // namespace ufo

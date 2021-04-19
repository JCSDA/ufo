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
#include "oops/util/missingValues.h"

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
    variablePredictors[p]->compute(odb_, geovals, ydiags, predData_[p]);
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

  const size_t nlocs  = ybiasinc.nlocs();
  const size_t nvars  = ybiasinc.nvars();
  const size_t npreds = predData_.size();

  // map bias coeffs to eigen matrix (writable)
  std::size_t indx;
  Eigen::VectorXd tmp = Eigen::VectorXd::Zero(nlocs, 1);
  for (std::size_t jch = 0; jch < nvars; ++jch) {
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      indx = jl*nvars+jch;
      if (ybiasinc[indx] != util::missingValue(ybiasinc[indx])) {
        tmp(jl) = ybiasinc[indx];
      }
    }
    // For each channel: ( npreds X 1 ) = ( npreds X nlocs ) * ( nlocs X 1 )
    for (std::size_t jp = 0; jp < npreds; ++jp) {
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        biascoeffinc(jp, jch) += predData_[jp][jl*nvars+jch] * tmp(jl);
      }
    }
    // zero out for next job
    tmp.setConstant(0.0);
  }

  // Sum across the processros
  biascoeffinc.allSumInPlace();

  oops::Log::trace() << "LinearObsBiasOperator::computeAD done." << std::endl;
}

// -----------------------------------------------------------------------------

void LinearObsBiasOperator::print(std::ostream & os) const {
  os << "TL/AD bias computed as linear combination of predictors.";
}

// -----------------------------------------------------------------------------

}  // namespace ufo

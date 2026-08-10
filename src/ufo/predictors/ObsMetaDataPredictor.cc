/*
 * (C) Copyright 2023 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <string>
#include <vector>

#include "ufo/predictors/ObsMetaDataPredictor.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/missingValues.h"

#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<ObsMetaDataPredictor> makerFuncObsMetaDataPredictor_(\
"obsMetadataPredictor");

// -----------------------------------------------------------------------------

ObsMetaDataPredictor::ObsMetaDataPredictor(const Parameters_ & parameters,
const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    order_(parameters.order.value().value_or(1)),
    variable_(parameters.varName) {
  // predictor name is a variable name
  name() = variable_;
  if (parameters.order.value() != boost::none) {
    // override the predictor name to distinguish between predictors of different orders
    name() = name() + "_order_" + std::to_string(order_);
  }
}

// -----------------------------------------------------------------------------

void ObsMetaDataPredictor::compute(const ioda::ObsSpace & odb,
                                   const GeoVaLs &,
                                   const ObsDiagnostics &,
                                   const ObsBias &,
                                   ioda::ObsVector & out) const {
  const size_t nlocs = out.nlocs();
  const size_t nvars = out.nvars();

  std::vector<float> obsMetaDataPred(nlocs, 0.0);

  const float fmiss = util::missingValue<float>();
  const int imiss = util::missingValue<int>();

  // retrieve the predictor

  if (odb.dtype("MetaData", variable_) == ioda::ObsDtype::Integer) {
    std::vector<int> obsMetaDataPred2(nlocs, 0);
    odb.get_db("MetaData", variable_, obsMetaDataPred2);
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      obsMetaDataPred[jloc] = (obsMetaDataPred2[jloc] == imiss)
                            ? fmiss
                            : static_cast<float>(obsMetaDataPred2[jloc])*1.0f;
    }
  } else {
    odb.get_db("MetaData", variable_, obsMetaDataPred);
  }

  const double dmiss = util::missingValue<double>();

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    const bool isMissing = (obsMetaDataPred[jloc] == fmiss);
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      out[jloc*nvars+jvar] = isMissing ? dmiss : std::pow(obsMetaDataPred[jloc], order_);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

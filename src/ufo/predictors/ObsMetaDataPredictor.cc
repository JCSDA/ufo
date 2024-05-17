/*
 * (C) Copyright 2023 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ufo/predictors/ObsMetaDataPredictor.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

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
  if (parameters.order.value() != boost::none) {
    // override the predictor name to distinguish between predictors of different orders
    // since this is a generalized function for any obs MetaData variable, it is prudent
    // to include the variable-name here as well to distinguish between different
    // variables that may be retrieved for the same bias correction task
    name() = name() + "_order_" + std::to_string(order_) + "_" +
             variable_;
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

  // retrieve the predictor

  if (odb.dtype("MetaData", variable_) == ioda::ObsDtype::Integer) {
    std::vector<int> obsMetaDataPred2(nlocs, 0);
    odb.get_db("MetaData", variable_, obsMetaDataPred2);
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      obsMetaDataPred[jloc] = static_cast<float>(obsMetaDataPred2[jloc])*1.0f;
    }
  } else {
    odb.get_db("MetaData", variable_, obsMetaDataPred);
  }

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      out[jloc*nvars+jvar] = pow(obsMetaDataPred[jloc], order_);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

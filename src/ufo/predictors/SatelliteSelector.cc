/*
 * (C) Crown copyright 2021, MetOffice
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ufo/predictors/SatelliteSelector.h"
#include "ufo/utils/Constants.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

namespace ufo {

static PredictorMaker<SatelliteSelector> makerFuncSatelliteSelector_("satelliteSelector");

// -----------------------------------------------------------------------------

SatelliteSelector::SatelliteSelector(const Parameters_ & parameters,
                                     const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars), predictor_(), satid_(parameters.satelliteId),
    metadata_name_(parameters.metadataName) {
  // Setup the predictor that will be run
  std::string predname = parameters.predictor.value().getString("name");
  std::unique_ptr<PredictorParametersBase> params(PredictorFactory::createParameters(predname));
  params->validateAndDeserialize(parameters.predictor.value());
  std::unique_ptr<PredictorBase> pred(PredictorFactory::create(*params, vars));
  predictor_ = std::move(pred);

  // Setup the name and variables for external access
  name() = predictor_->name() + "_satid_" + std::to_string(satid_);
  geovars_ = predictor_->requiredGeovars();
  hdiags_ = predictor_->requiredHdiagnostics();
}

// -----------------------------------------------------------------------------

void SatelliteSelector::compute(const ioda::ObsSpace & odb,
                                const GeoVaLs & gv,
                                const ObsDiagnostics & obsdiags,
                                const ObsBias & biascoeffs,
                                ioda::ObsVector & out) const {
  const size_t nlocs = out.nlocs();
  const size_t nvars = out.nvars();

  predictor_->compute(odb, gv, obsdiags, biascoeffs, out);

  std::vector<int> satid(nlocs);
  odb.get_db("MetaData", metadata_name_, satid);

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    if (satid[jloc] != satid_) {
      for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
        out[jloc*nvars+jvar] = Constants::zero;
      }  // jvar
    }
  }  // jloc
}

// -----------------------------------------------------------------------------

}  // namespace ufo

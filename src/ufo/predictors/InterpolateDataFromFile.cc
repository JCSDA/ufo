/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/predictors/InterpolateDataFromFile.h"

#include <set>
#include <string>
#include <vector>

#include <boost/make_unique.hpp>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/AssociativeContainers.h"
#include "ufo/filters/ObsFilterData.h"

namespace ufo {

namespace {

/// \brief Store all components of `input` in `output`.
void save(const ioda::ObsDataVector<float> &input, ioda::ObsVector &output) {
  const float fmiss = util::missingValue<float>();
  const double dmiss = util::missingValue<double>();
  const size_t numLocs = input.nlocs();
  const size_t numInputVars = input.nvars();
  const size_t numOutputVars = output.nvars();
  for (size_t inputVarIndex = 0; inputVarIndex < numInputVars; ++inputVarIndex) {
    const ioda::ObsDataRow<float> &inputVarValues = input[inputVarIndex];
    const size_t outputVarIndex = output.varnames().find(input.varnames()[inputVarIndex]);
    for (size_t locIndex = 0; locIndex < numLocs; ++locIndex) {
      const float value = inputVarValues[locIndex];
      output[locIndex * numOutputVars + outputVarIndex] =
          (value == fmiss) ? dmiss : static_cast<double>(value);
    }
  }
}

std::set<std::string> getVariableNamesWithoutChannels(const oops::ObsVariables &vars) {
  if (vars.channels().empty()) {
    return std::set<std::string>(vars.variables().begin(), vars.variables().end());
  } else {
    // We assume there's only a single multi-channel variable
    ASSERT(vars.channels().size() == vars.variables().size());
    std::string channellessVariable = vars.variables().front();
    channellessVariable.resize(channellessVariable.find_last_of('_'));
    return std::set<std::string>{channellessVariable};
  }
}

}  // namespace

static PredictorMaker<InterpolateDataFromFile> maker("interpolate_data_from_file");

InterpolateDataFromFile::InterpolateDataFromFile(const Parameters_ & parameters,
                                                 const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars) {
  const std::set<std::string> channellessVariables = getVariableNamesWithoutChannels(vars_);

  for (const VariableCorrectionParameters & varParams : parameters.correctedVariables.value()) {
    if (!oops::contains(channellessVariables, varParams.name))
      throw eckit::UserError("'" + varParams.name.value() +
                             "' is not in the list of bias-corrected variables", Here());
    eckit::LocalConfiguration varConfig = varParams.details.toConfiguration();
    varConfig.set("group", "ObsBias");
    obsFunctions_[varParams.name] = boost::make_unique<DrawValueFromFile<float>>(varConfig);
  }

  for (const auto &varAndObsFunction : obsFunctions_) {
    const DrawValueFromFile<float> &obsFunction = *varAndObsFunction.second;
    const ufo::Variables &requiredVariables = obsFunction.requiredVariables();
    geovars_ += requiredVariables.allFromGroup("GeoVaLs").toOopsVariables();
    hdiags_ += requiredVariables.allFromGroup("ObsDiag").toOopsObsVariables();
  }
}

void InterpolateDataFromFile::compute(const ioda::ObsSpace & /*odb*/,
                                      const GeoVaLs & geovals,
                                      const ObsDiagnostics & obsdiags,
                                      const ObsBias &,
                                      ioda::ObsVector & out) const {
  ObsFilterData obsFilterData(out.space());
  obsFilterData.associate(geovals);
  obsFilterData.associate(obsdiags);

  out.zero();

  for (const auto &varAndObsFunction : obsFunctions_) {
    const std::string &varName = varAndObsFunction.first;
    const DrawValueFromFile<float> &obsFunction = *varAndObsFunction.second;

    oops::ObsVariables currentVars({varName}, vars_.channels());
    ioda::ObsDataVector<float> obsFunctionResult(out.space(), currentVars);
    obsFunction.compute(obsFilterData, obsFunctionResult);
    save(obsFunctionResult, out);
  }
}

}  // namespace ufo

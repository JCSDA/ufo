/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/DrawObsErrorFromFile.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/PropertiesOfNVectors.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

namespace {

eckit::LocalConfiguration makeConfigForDrawValueFromFile(const eckit::LocalConfiguration &config,
                                                         const std::string &group) {
  eckit::LocalConfiguration newConfig(config);
  newConfig.set("group", group);
  return newConfig;
}

}  // namespace

// -----------------------------------------------------------------------------
constexpr char DispersionMeasureParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<DispersionMeasure>
  DispersionMeasureParameterTraitsHelper::namedValues[];

static ObsFunctionMaker<DrawObsErrorFromFile> maker("DrawObsErrorFromFile");


DrawObsErrorFromFile::DrawObsErrorFromFile(const eckit::LocalConfiguration &config)
  : invars_() {
  options_.reset(new DrawObsErrorFromFileParameters());
  options_->deserialize(config);

  // Initialise the DrawValueFromFile object with the `group` option.
  drawValueFromFile_.reset(new DrawValueFromFile<float>
                           (makeConfigForDrawValueFromFile(config, options_->group.value())));

  if (options_->normvariable.value() != boost::none) {
    const boost::optional<Variable> &normvariable = options_->normvariable.value();
    invars_ += *normvariable;
    oops::Log::debug() << "DrawObsErrorFromFile: norm variable = " << *normvariable << std::endl;
  }
}

void DrawObsErrorFromFile::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<float> & out) const {
  const float missing = util::missingValue(missing);

  // Interpolate variances
  drawValueFromFile_->compute(in, out);

  // Transform variances into standard deviations if required
  if (options_->dispersionMeasure.value() == DispersionMeasure::VARIANCE) {
    for (size_t ivar = 0; ivar < out.nvars(); ++ivar)
      for (float &value : out[ivar])
        if (value != missing)
          value = std::sqrt(value);
  }

  // Transform normalized standard deviations to standard deviation if required
  if ((options_->dispersionMeasure.value() == DispersionMeasure::NORMALIZED ||
       options_->dispersionMeasure.value() == DispersionMeasure::FRACTIONAL) &&
      options_->normvariable.value() != boost::none) {
    const boost::optional<Variable> &normvariable = options_->normvariable.value();

    oops::Log::debug() << "Norm variable is: "
                       << (*normvariable).variable() << "  and name: "
                       << (*normvariable).fullName() << "  and group: "
                       << (*normvariable).group() << std::endl;

    ioda::ObsDataVector<float> obvalues(in.obsspace(), (*normvariable).toOopsVariables());
    in.get(*normvariable, obvalues);

    oops::Log::debug() << "Sizes are: obvalues=" << obvalues[0].size() <<
                          "  value size= " << out[0].size() << std::endl;

    // Normalized error is in % so the final obs error value should be divided by 100
    // If fractional is chosen, then no normalization is needed
    float normFactor = 100.;
    if (options_->dispersionMeasure.value() == DispersionMeasure::FRACTIONAL)
      normFactor = 1.;

    for (size_t ivar = 0; ivar < out.nvars(); ++ivar) {
      if (obvalues[ivar].size() == out[ivar].size()) {
        for (size_t iloc = 0; iloc < obvalues[ivar].size(); ++iloc) {
          if (out[ivar][iloc] != missing) {
            out[ivar][iloc] = out[ivar][iloc] * std::abs(obvalues[ivar][iloc]) / normFactor;
            out[ivar][iloc] = std::max(options_->minValue.value(), out[ivar][iloc]);
          }
        }
      }
    }
  }
}

const ufo::Variables & DrawObsErrorFromFile::requiredVariables() const {
  return drawValueFromFile_->requiredVariables();
}

}  // namespace ufo

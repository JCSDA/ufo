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
#include "oops/util/IntSetParser.h"
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
  options_.deserialize(config);

  // Initialise the DrawValueFromFile object with the `group` option.
  drawValueFromFile_.reset(new DrawValueFromFile<float>
                           (makeConfigForDrawValueFromFile(config, options_.group.value())));

  if (options_.normvariable.value() != boost::none) {
    const boost::optional<Variable> &normvariable = options_.normvariable.value();
    invars_ += *normvariable;
    oops::Log::debug() << "DrawObsErrorFromFile: norm variable = " << *normvariable << std::endl;
  }

  // Get channels from options
  if (options_.chlist.value() != boost::none) {
      std::set<int> channels = options_.chlist.value().get();
      channels_ = {std::make_move_iterator(std::begin(channels)),
                   std::make_move_iterator(std::end(channels))};
  }
}

void DrawObsErrorFromFile::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<float> & out) const {
  const float missing = util::missingValue<float>();

  // Interpolate variances
  drawValueFromFile_->compute(in, out);

  // Transform variances into standard deviations if required
  if (options_.dispersionMeasure.value() == DispersionMeasure::VARIANCE) {
    for (size_t ivar = 0; ivar < out.nvars(); ++ivar)
      for (float &value : out[ivar])
        if (value != missing)
          value = std::sqrt(value);
  }

  // DrawValueFromFile can only work on one variable (you can only draw one variable
  // at a time).  Therefore the number of channels must equal the number of variables.
  if (channels_.size() > 0 && out.nvars() != channels_.size()) {
    throw eckit::BadValue("The number of channels must match the number of variables " +
      std::to_string(channels_.size()) + " != " + std::to_string(out.nvars()), Here());
  }

  // Transform normalized standard deviations to standard deviation if required
  if ((options_.dispersionMeasure.value() == DispersionMeasure::NORMALIZED ||
       options_.dispersionMeasure.value() == DispersionMeasure::FRACTIONAL) &&
      options_.normvariable.value() != boost::none) {
    const boost::optional<Variable> &normvariable = options_.normvariable.value();

    oops::Log::debug() << "Norm variable is: "
                       << (*normvariable).variable() << "  and name: "
                       << (*normvariable).fullName() << "  and group: "
                       << (*normvariable).group() << std::endl;

    // Normalized error is in % so the final obs error value should be divided by 100
    // If fractional is chosen, then no normalization is needed
    float normFactor = 100.;
    if (options_.dispersionMeasure.value() == DispersionMeasure::FRACTIONAL)
      normFactor = 1.;

    for (size_t ivar = 0; ivar < out.nvars(); ++ivar) {
      std::vector<float> obvalues;
      in.get(Variable((*normvariable).fullName(), channels_)[ivar], obvalues);
      if (obvalues.size() == out[ivar].size()) {
        for (size_t iloc = 0; iloc < obvalues.size(); ++iloc) {
          if (out[ivar][iloc] != missing && obvalues[iloc] != missing) {
            out[ivar][iloc] = out[ivar][iloc] * std::abs(obvalues[iloc]) / normFactor;
            out[ivar][iloc] = std::max(options_.minValue.value(), out[ivar][iloc]);
          } else {
            out[ivar][iloc] = missing;
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

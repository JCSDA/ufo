/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/DrawObsErrorFromFile.h"

#include "ioda/ObsDataVector.h"

namespace ufo {

namespace {

eckit::LocalConfiguration makeConfigForDrawValueFromFile(const eckit::LocalConfiguration &config) {
  eckit::LocalConfiguration newConfig(config);
  newConfig.set("group", "ErrorVariance");
  return newConfig;
}

}  // namespace

static ObsFunctionMaker<DrawObsErrorFromFile> maker("DrawObsErrorFromFile");


DrawObsErrorFromFile::DrawObsErrorFromFile(const eckit::LocalConfiguration &config)
  : drawValueFromFile_(makeConfigForDrawValueFromFile(config)) {
}

void DrawObsErrorFromFile::compute(const ObsFilterData & in,
                                   ioda::ObsDataVector<float> & out) const {
  const float missing = util::missingValue(missing);

  // Interpolate variances
  drawValueFromFile_.compute(in, out);

  // Transform them into standard deviations
  for (size_t ivar = 0; ivar < out.nvars(); ++ivar)
    for (float &value : out[ivar])
      if (value != missing)
        value = std::sqrt(value);
}

const ufo::Variables & DrawObsErrorFromFile::requiredVariables() const {
  return drawValueFromFile_.requiredVariables();
}

}  // namespace ufo

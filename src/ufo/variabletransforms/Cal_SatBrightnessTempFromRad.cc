/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>

#include "ufo/variabletransforms/Cal_SatBrightnessTempFromRad.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {
constexpr char RadianceUnitsParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<RadianceUnits> RadianceUnitsParameterTraitsHelper::namedValues[];

/************************************************************************************/
//  Cal_SatBrightnessTempFromRad
/************************************************************************************/

static TransformMaker<Cal_SatBrightnessTempFromRad>
    makerCal_SatBrightnessTempFromRad_("SatBrightnessTempFromRad");

Cal_SatBrightnessTempFromRad::Cal_SatBrightnessTempFromRad(
    const Parameters_ &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr), parameters_(options),
      variables_(), channels_(parameters_.transformVariable.value().channels()) {
  variables_ += parameters_.transformVariable.value();
  variables_ += parameters_.spectralVariable.value();
  ASSERT(channels_ == parameters_.spectralVariable.value().channels());
}

/************************************************************************************/

void Cal_SatBrightnessTempFromRad::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> convert input array to brightness temperature"
            << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  // Read in radiance to be corrected
  oops::Variables radianceVar(parameters_.transformVariable.value().toOopsVariables());
  ioda::ObsDataVector<float> radiance(obsdb_, radianceVar);
  data_.get(parameters_.transformVariable.value(), radiance);

  // Read in the spectral variable
  oops::Variables spectralVar(parameters_.spectralVariable.value().toOopsVariables());
  ioda::ObsDataVector<float> spectralVariable(obsdb_, spectralVar);
  data_.get(parameters_.spectralVariable.value(), spectralVariable);

  // Setup Variables
  const size_t nlocs = obsdb_.nlocs();
  const size_t nvars = radiance.nvars();
  float minval = parameters_.minvalue.value().value_or(missingValueFloat);
  float maxval = parameters_.maxvalue.value().value_or(missingValueFloat);

  // Create output array
  std::vector<std::vector<float>> brightnessTemperature(
                         nvars, std::vector<float>(nlocs, missingValueFloat));

  // Sanity check
  ASSERT(radiance.nlocs() == nlocs);

  // Loop over all obs
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (apply[iloc]) {
      for (size_t ichan = 0; ichan < nvars; ++ichan) {
        if ( radiance[ichan][iloc] != missingValueFloat &&
             radiance[ichan][iloc] > 0.0f &&
             spectralVariable[ichan][iloc] != missingValueFloat) {
          double bt = missingValueFloat;
          switch (parameters_.radianceUnits.value()) {
            case RadianceUnits::WAVENUMBER: {
              bt = formulas::inversePlanck(static_cast<double>(radiance[ichan][iloc]),
                                           static_cast<double>(spectralVariable[ichan][iloc]));
              break;
            }
            case RadianceUnits::FREQUENCY: {
              double freq = static_cast<double>(spectralVariable[ichan][iloc]);
              double wvn = freq / Constants::speedOfLight;  // Hz to m-1
              double rad = static_cast<double>(radiance[ichan][iloc]) * freq / wvn;
              bt = formulas::inversePlanck(rad, wvn);
              break;
            }
            case RadianceUnits::WAVELENGTH: {
              double wvl = static_cast<double>(spectralVariable[ichan][iloc]);
              double wvn = 1.0e6 / wvl;  // microns to m-1
              double rad = static_cast<double>(radiance[ichan][iloc]) * wvl / wvn;
              bt = formulas::inversePlanck(rad, wvn);
              break;
            }
          }
          brightnessTemperature[ichan][iloc] = static_cast<float>(bt);
          if (minval != missingValueFloat)
            if (brightnessTemperature[ichan][iloc] < minval)
              brightnessTemperature[ichan][iloc] = missingValueFloat;
          if (maxval != missingValueFloat)
            if (brightnessTemperature[ichan][iloc] > maxval)
              brightnessTemperature[ichan][iloc] = missingValueFloat;
        }
      }  // ichan
    }  // apply
  }  // iloc

  //  Write out the resulting data to Derived group and update qcflags
  for (size_t ichan =0; ichan < nvars; ++ichan) {
    putObservation("brightness_temperature_" + std::to_string(channels_[ichan]),
                   brightnessTemperature[ichan]);
  }
}  // runTransform
}  // namespace ufo

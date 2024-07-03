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
  oops::ObsVariables radianceVar(parameters_.transformVariable.value().toOopsObsVariables());
  ioda::ObsDataVector<float> radiance(obsdb_, radianceVar);
  data_.get(parameters_.transformVariable.value(), radiance);

  // Setup Variables
  const size_t nlocs = obsdb_.nlocs();
  const size_t nvars = radiance.nvars();
  float minval = parameters_.minvalue.value().value_or(missingValueFloat);
  float maxval = parameters_.maxvalue.value().value_or(missingValueFloat);

  // Read in the spectral variable
  const ufo::Variable & var = parameters_.spectralVariable.value();
  std::vector<float> spectralVariable(nvars, missingValueFloat);
  bool spectralVariableCreated = false;
  // If variable is size `Channel`
  if (obsdb_.has(var.group(), var.variable())) {
    std::vector<float> bufr;
    data_.get(var, bufr);
    if (bufr.size() == nvars) {
      for (size_t ichan = 0; ichan < nvars; ++ichan) {
        spectralVariable[ichan] = bufr[ichan];
      }
      spectralVariableCreated = true;
    }
  }
  // If variable includes `Location` dimension
  if (!spectralVariableCreated) {
    ioda::ObsDataVector<float> bufr(obsdb_,
                                    parameters_.spectralVariable.value().toOopsObsVariables());
    data_.get(parameters_.spectralVariable.value(), bufr);
    for (size_t ichan = 0; ichan < nvars; ++ichan) {
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        if (bufr[ichan][iloc] != missingValueFloat) {
          spectralVariable[ichan] = bufr[ichan][iloc];
          break;
        }
      }  // ivar
    }  // ichan
  }

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
             spectralVariable[ichan] != missingValueFloat) {
          double bt = missingValueFloat;
          switch (parameters_.radianceUnits.value()) {
            case RadianceUnits::WAVENUMBER: {
              bt = formulas::inversePlanck(static_cast<double>(radiance[ichan][iloc]),
                                           static_cast<double>(spectralVariable[ichan]),
                                           parameters_.planck1.value(),
                                           parameters_.planck2.value());
              break;
            }
            case RadianceUnits::FREQUENCY: {
              double freq = static_cast<double>(spectralVariable[ichan]);
              double wvn = freq / Constants::speedOfLight;  // Hz to m-1
              double rad = static_cast<double>(radiance[ichan][iloc]) * freq / wvn;
              bt = formulas::inversePlanck(rad, wvn, parameters_.planck1.value(),
                                           parameters_.planck2.value());
              break;
            }
            case RadianceUnits::WAVELENGTH: {
              double wvl = static_cast<double>(spectralVariable[ichan]);
              double wvn = 1.0e6 / wvl;  // microns to m-1
              double rad = static_cast<double>(radiance[ichan][iloc]) * wvl / wvn;
              bt = formulas::inversePlanck(rad, wvn, parameters_.planck1.value(),
                                           parameters_.planck2.value());
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
    putObservation("brightnessTemperature_" + std::to_string(channels_[ichan]),
                   brightnessTemperature[ichan]);
  }
}  // runTransform
}  // namespace ufo

/*
 * (C) Crown copyright 2022, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>

#include "ufo/variabletransforms/Cal_SatRadianceFromScaledRadiance.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace ufo {

/************************************************************************************/
//  Cal_SatRadianceFromScaledRadiance
/************************************************************************************/

static TransformMaker<Cal_SatRadianceFromScaledRadiance>
    makerCal_SatRadianceFromScaledRadiance_("SatRadianceFromScaledRadiance");

Cal_SatRadianceFromScaledRadiance::Cal_SatRadianceFromScaledRadiance(
    const Parameters_ &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr), parameters_(options),
      variables_({parameters_.transformVariable.value()}),
      channels_(parameters_.transformVariable.value().channels()) {
  size_t numScaleFactors = parameters_.numScaleFactors.value();
  for (std::size_t i = 0; i < numScaleFactors; ++i) {
    std::string str = std::to_string(i+1);
    variables_ += Variable(parameters_.scalingVariable.value().fullName() + str);
    variables_ += Variable(parameters_.scalingStartChannel.value().fullName() + str);
    variables_ += Variable(parameters_.scalingEndChannel.value().fullName() + str);
  }
}

/************************************************************************************/

void Cal_SatRadianceFromScaledRadiance::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> convert input array to brightness temperature"
            << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  size_t numScaleFactors = parameters_.numScaleFactors.value();

  // Read in radiance to be corrected
  oops::ObsVariables radianceVar(parameters_.transformVariable.value().toOopsObsVariables());
  ioda::ObsDataVector<float> radiance(obsdb_, radianceVar);
  data_.get(parameters_.transformVariable.value(), radiance);

  // Read in scaling factors
  std::vector<int> channelScaleFactor(numScaleFactors, missingValueInt);
  std::vector<int> startChannelScale(numScaleFactors, missingValueInt);
  std::vector<int> endChannelScale(numScaleFactors, missingValueInt);
  if (parameters_.getFactorsFromMultipleArrays) {
    for (std::size_t i = 0; i < channelScaleFactor.size(); ++i) {
      std::string str = std::to_string(i+1);
      getFirstLocationValue(Variable(parameters_.scalingVariable.value().fullName() + str),
                            channelScaleFactor[i]);
      getFirstLocationValue(Variable(parameters_.scalingStartChannel.value().fullName() + str),
                            startChannelScale[i]);
      getFirstLocationValue(Variable(parameters_.scalingEndChannel.value().fullName() + str),
                            endChannelScale[i]);
    }
  } else {
    Variable var = parameters_.scalingVariable.value();
    obsdb_.get_db(var.group(), var.variable(), channelScaleFactor);
    ASSERT(channelScaleFactor.size() == numScaleFactors);
    var = parameters_.scalingStartChannel.value();
    obsdb_.get_db(var.group(), var.variable(), startChannelScale);
    ASSERT(startChannelScale.size() == numScaleFactors);
    var = parameters_.scalingEndChannel.value();
    obsdb_.get_db(var.group(), var.variable(), endChannelScale);
    ASSERT(endChannelScale.size() == numScaleFactors);
  }

  // Setup Variables
  const size_t nlocs = obsdb_.nlocs();
  const size_t nvars = radiance.nvars();

  // Sanity check
  ASSERT(radiance.nlocs() == nlocs);

  // Loop over all obs
  // takes a scaled radiance and produces a radiance in (W / (m^2.sr.m^-1))
  // radiance is power / (area . solid angle . wavenumber)
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (apply[iloc]) {
      for (size_t ichan = 0; ichan < nvars; ++ichan) {
        if (radiance[ichan][iloc] != missingValueFloat  &&
            radiance[ichan][iloc] > 0.0f) {
          for (size_t iscale = 0; iscale < numScaleFactors; ++iscale) {
            if (channels_[ichan] >= startChannelScale[iscale] &
                channels_[ichan] <= endChannelScale[iscale]) {
              radiance[ichan][iloc] *= std::pow(10, (-1.0f*channelScaleFactor[iscale]));
              break;
            }  // if channels_
          }  // iscale
        }  // if missing
      }  // ichan
    }  // apply
  }  // iloc

  //  Write out the resulting data to Derived group and update qcflags
  for (size_t ichan =0; ichan < nvars; ++ichan) {
    putObservation("radiance_" + std::to_string(channels_[ichan]),
                   radiance[ichan]);
  }
}  // runTransform
}  // namespace ufo

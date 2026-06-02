/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>

#include "oops/base/Variables.h"

#include "ufo/variabletransforms/Cal_SatZenithAngleCorrection.h"

namespace ufo {

/************************************************************************************/
//  Cal_SatZenithAngleCorrection
/************************************************************************************/

static TransformMaker<Cal_SatZenithAngleCorrection>
    makerCal_SatZenithAngleCorrection_("SatZenithAngleCorrection");

Cal_SatZenithAngleCorrection::Cal_SatZenithAngleCorrection(
    const Parameters_ &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr), parameters_(options),
      variables_({parameters_.transformVariable.value()}) {
  ASSERT(parameters_.coeffA.value().size() == parameters_.transformVariable.value().size());
  ASSERT(parameters_.coeffB.value().size() == parameters_.transformVariable.value().size());
  ASSERT(parameters_.coeffC.value().size() == parameters_.transformVariable.value().size());
}

/************************************************************************************/

void Cal_SatZenithAngleCorrection::runTransform(const std::vector<bool> &apply) {
  Variable variableToBeCorrected(parameters_.transformVariable.value());
  oops::Log::trace() << " --> correct input array for the satellite viewing angle"
            << std::endl;
  oops::Log::trace() << "      --> method: " << options_.Method.value() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;
  oops::Log::trace() << "      --> variable being corrected: "
                     << variableToBeCorrected.fullName() << std::endl;

  // Read in variable to be corrected
  oops::ObsVariables varin(variableToBeCorrected.toOopsObsVariables());
  ioda::ObsDataVector<float> varArray(obsdb_, varin);
  data_.get(variableToBeCorrected, varArray);

  // Read in zenith angle
  std::vector<float> zenithAngle;
  getObservation("MetaData", "sensorZenithAngle", zenithAngle, true);

  // Setup Variables
  const size_t nlocs = obsdb_.nlocs();
  const size_t nvars = varArray.nvars();

  // Sanity check
  ASSERT(zenithAngle.size() == nlocs);

  // Loop over all obs
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (apply[iloc]) {
      for (size_t ichan = 0; ichan < nvars; ++ichan) {
        if ( varArray[ichan][iloc] != missingValueFloat &&
             zenithAngle[iloc] != missingValueFloat) {
          varArray[ichan][iloc] += parameters_.coeffA.value()[ichan] *
                                   std::pow(zenithAngle[iloc], parameters_.exponentA.value());
          varArray[ichan][iloc] += parameters_.coeffB.value()[ichan] *
                                   std::pow(zenithAngle[iloc], parameters_.exponentB.value());
          varArray[ichan][iloc] += parameters_.coeffC.value()[ichan] *
                                   std::pow(zenithAngle[iloc], parameters_.exponentC.value());
          if (parameters_.maxvalue.value() != boost::none)
            varArray[ichan][iloc] = std::min(varArray[ichan][iloc],
                                             parameters_.maxvalue.value().get());
          if (parameters_.minvalue.value() != boost::none)
            varArray[ichan][iloc] = std::max(varArray[ichan][iloc],
                                             parameters_.minvalue.value().get());
        } else {
          varArray[ichan][iloc] = missingValueFloat;
        }
      }  // ichan
    }
  }  // iloc

  // Overwrite variable at existing location
  varArray.save(variableToBeCorrected.group());
}
}  // namespace ufo

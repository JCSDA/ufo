/* * (C) Copyright 2024 NASA GMAO
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/AverageObsToGeoValLevels.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"
#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// AverageObsToGeoValLevels: average the observation values on to model levels.

AverageObsToGeoValLevels::AverageObsToGeoValLevels(
        ioda::ObsSpace & obsdb,
        const Parameters_ & parameters,
        std::shared_ptr<ioda::ObsDataVector<int> > flags,
        std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "AverageObsToGeoValLevels constructor" << std::endl;

  // Ensure observations have been grouped into profiles.
  if (obsdb_.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Check the ObsSpace has been extended. If this is not the case
  // then it will not be possible to access profiles in the original and
  // extended sections of the ObsSpace.
  if (!obsdb_.has("MetaData", "extendedObsSpace"))
    throw eckit::UserError("The extended obs space has not been produced", Here());

  // Get parameters from options
  allvars_ += parameters_.observation_vertical_coordinate;
  allvars_ += parameters_.model_vertical_coordinate;
}

// -----------------------------------------------------------------------------

AverageObsToGeoValLevels::~AverageObsToGeoValLevels() {
  oops::Log::trace() << "AverageObsToGeoValLevels destructor" << std::endl;
}

// -----------------------------------------------------------------------------
/// Apply the Average Observations To Model Levels filter.
void AverageObsToGeoValLevels::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "AverageObsToGeoValLevels applyFilter start" << std::endl;
  oops::Log::debug() << "AverageObsToGeoValLevels obserr: " << *obserr_ << std::endl;

  const float missing = util::missingValue<float>();
  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsObsVariables(), "ObsValue");

  // Number of locations (including extended space)
  const size_t nlocs = obsdb_.nlocs();

  oops::Log::debug() << "Number of obervation locations " << nlocs << std::endl;
  // Get model vertical coordinate name
  const oops::Variable model_vcoord_name = oops::Variable
                   (Variable(parameters_.model_vertical_coordinate).variable());

  // Get GeoVals
  const ufo::GeoVaLs * gvals = data_.getGeoVaLs();

  // Get obs vertical coordinate
  std::vector<float> obs_vert_coord_all(nlocs);
  data_.get(parameters_.observation_vertical_coordinate, obs_vert_coord_all);

  // Get the sequence numbers from the observation data.  These will be used to identify
  // which observations belong to which profile.
  std::vector<float> sequenceNumber(nlocs, 0.0);
  std::vector<float> modelLayer(nlocs, missing);
  std::vector<int> actObsAvgQC(nlocs, 0);
  std::vector<float> height(obs_vert_coord_all.begin(), obs_vert_coord_all.end());

  if (obsdb_.has("MetaData", "sequenceNumber")) {
     obsdb_.get_db("MetaData", "sequenceNumber", sequenceNumber);
  }

  if (obsdb_.has("MetaData", "obsLayer")) {
     obsdb_.get_db("MetaData", "obsLayer", modelLayer);
  }

  const std::vector<size_t> &uniqueSeqNums = obsdb_.recidx_all_recnums();
  oops::Log::debug() << "Unique record numbers" << uniqueSeqNums << std::endl;
  oops::Log::debug() << "Number of unique record numbers" << uniqueSeqNums.size() << std::endl;

  // Number of profiles in the original ObsSpace.
  const size_t numUniqSeqNum = uniqueSeqNums.size() / 2;

  for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
    const std::string varname = filtervars.variable(jv).variable();
    oops::Log::debug() << "varname: " << varname << std::endl;

    // Loop over profiles.
    for (size_t iseq = 0; iseq < numUniqSeqNum; ++iseq) {
      // Get locations of profile in the original ObsSpace and
      // the corresponding profile in the extended ObsSpace.
      // Assuming the extended ObsSpace has been configured correctly, which is
      // checked above, the profile in the extended ObsSpace is always located
      // nprofs positions further on than the profile in the original ObsSpace.
      const size_t extSeqNum = uniqueSeqNums[iseq + numUniqSeqNum];
      const std::vector<size_t> &locsOriginal = obsdb_.recidx_vector(uniqueSeqNums[iseq]);
      const std::vector<size_t> &locsExt = obsdb_.recidx_vector(extSeqNum);

      // Get GeoVaLs at the specified location. Note that all the datapoints in the original
      // obsspace are using the same model profile
      std::vector<float> modelZ;
      modelZ.assign(gvals->nlevs(model_vcoord_name), 0.0);
      gvals->getAtLocation(modelZ, model_vcoord_name, locsOriginal[0]);

      const bool nonZeroModelZ = std::any_of(modelZ.begin(), modelZ.end(),
                                     [](float z1) { return std::abs(z1) > 0.0; });
      ASSERT(nonZeroModelZ && "Error: All model levels are zero!");

      // Loop through all the observations
      const size_t nModelLayers = modelZ.size() - 1;
      const size_t nLocsObs = locsOriginal.size();

      // New vector to hold the selected elements
      std::vector<float> obs_vert_coord;
      obs_vert_coord.reserve(nLocsObs);
      for (size_t iloc : locsOriginal) {
          obs_vert_coord.push_back(obs_vert_coord_all[iloc]);
      }

      // Check that model layers and extended obs space are the same
      ASSERT(nModelLayers == locsExt.size() &&
            "Error: number of model layers and extended ObsSpace have different size");

      for (size_t j = 0; j < nModelLayers; ++j) {
          float sumobs = static_cast<float>(ufo::Constants::zero);
          size_t numobs = static_cast<size_t>(ufo::Constants::zero);
          const size_t iLocExt = locsExt[j];
          modelLayer[iLocExt] = j;
          height[iLocExt] = 0.5 * (modelZ[j] + modelZ[j + 1]);
          sequenceNumber[iLocExt] = extSeqNum;
          for (size_t i = 0; i < nLocsObs; ++i) {
           const size_t iloc = locsOriginal[i];
           const float obsZ = obs_vert_coord[i];
            // Check if the observation is between model level j and j+1
            if ((obsZ < modelZ[j]) && (obsZ >= modelZ[j + 1])) {
             if ((obs[jv][iloc] != missing)  && !flagged[jv][iloc]) {
                sumobs += obs[jv][iloc];
                numobs++;
              }
            }
          }
          if (numobs > 0) {
           obs[jv][iLocExt] = sumobs / numobs;
            actObsAvgQC[iLocExt] = 1;
          } else {
            obs[jv][iLocExt] = missing;
          }
      }
    }  // record iseq
    // Save new obs to obsSpace:
    obsdb_.put_db("DerivedObsValue", varname, obs[jv]);
  }  // filter variable jv

  obsdb_.put_db("MetaData", "sequenceNumber", sequenceNumber);
  obsdb_.put_db("MetaData", "modelLayer", modelLayer);
  obsdb_.put_db("MetaData", "height", height);
  obsdb_.put_db("MetaData", "actObsAvgQC", actObsAvgQC);
  oops::Log::debug() << "Total number of locations in  ObsSpace: " << obsdb_.nlocs() << std::endl;
  oops::Log::trace() << "AverageObsToGeoValLevels applyFilter complete" << std::endl;
}

// -----------------------------------------------------------------------------

void AverageObsToGeoValLevels::print(std::ostream & os) const {
  os << "AverageObsToGeoValLevels: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

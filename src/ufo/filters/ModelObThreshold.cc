/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ModelObThreshold.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"
#include "oops/util/PropertiesOfNVectors.h"

#include "ufo/GeoVaLs.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"

namespace ufo {

constexpr char ThresholdTypeParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<ThresholdType>
  ThresholdTypeParameterTraitsHelper::namedValues[];

// -----------------------------------------------------------------------------

ModelObThreshold::ModelObThreshold(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                               std::shared_ptr<ioda::ObsDataVector<int> > flags,
                               std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "ModelObThreshold contructor starting" << std::endl;
  allvars_ += parameters_.model_profile;
  allvars_ += parameters_.model_vcoord;
  allvars_ += parameters_.obs_height;
}

// -----------------------------------------------------------------------------

ModelObThreshold::~ModelObThreshold() {
  oops::Log::trace() << "ModelObThreshold destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief Filter to apply a threshold to a model profile interpolated to the
 * observation height.
 *
 * \details The specified model profile variable is linearly (vertical) interpolated
 * to the observation height using the specified model vertical coordinate variable.
 * This is referred to as the "ModelOb". Note that the ModelOb is not necessarily
 * one of the HofX variables.
 *
 * The observation height must be in the same coordinate system as that specified
 * for the model vertical coordinate, e.g. both pressure.
 *
 * The thresholds to compare the ModelOb against is specified as height-dependent.
 * We supply a vector of threshold values, and a vector of vertical coordinate
 * values corresponding to those thresholds. The coordinate values must be in the same
 * vertical coordinate as the observation, e.g. pressure. The threshold values are
 * then linearly interpolated to the observation height.
 *
 * The observation is flagged for rejection if the ModelOb lies outside the threshold
 * value according to threshold type - min or max. E.g. if the threshold type is min,
 * then the observation is flagged if ModelOb is less than the interpolated threshold
 * value.
 *
 * Example for relative humidity:
 * \code{.unparsed}
 *  obs filters:
 *  - filter: ModelOb Threshold
 *    model profile:
 *      name: GeoVaLs/relative_humidity
 *    model vertical coordinate:
 *      name: GeoVaLs/air_pressure
 *    observation height:
 *      name: MetaData/pressure
 *    thresholds: [50,50,40,30]
 *    coordinate values: [100000,80000,50000,20000]
 *    threshold type: min
 * \endcode
 *
 * \author J.Cotton (Met Office)
 *
 * \date 12/03/2021: Created
 */

void ModelObThreshold::applyFilter(const std::vector<bool> & apply,
                                   const Variables & filtervars,
                                   std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ModelObThreshold priorFilter" << std::endl;
  print(oops::Log::trace());

  const float missing = util::missingValue<float>();
  const size_t nlocs = obsdb_.nlocs();

// Get piece-wise parameters from options.
  const std::vector<double> coord_vals = parameters_.coord_vals.value();
  const std::vector<double> thresholds = parameters_.thresholds.value();
  oops::Log::debug() << "QC coord vals are " << coord_vals << std::endl;
  oops::Log::debug() << "QC thresholds are " << thresholds << std::endl;

// get names of GeoVal variables
  const oops::Variable model_profile_name = oops::Variable
                                                   (Variable(parameters_.model_profile).variable());
  const oops::Variable model_vcoord_name = oops::Variable
                                                    (Variable(parameters_.model_vcoord).variable());

  std::ostringstream errString;
// Ensure same size vectors (coord_vals and threshold); Also ensure more than one value in each.
  if (coord_vals.size() <= 1 || coord_vals.size() != thresholds.size()) {
      errString << "At least one of coord_vals, thresholds is wrong size - either unequal or < 2"
                << oops::listOfVectorSizes(coord_vals, thresholds) << std::endl;
      throw eckit::BadValue(errString.str());
  }

// Get variables from ObsSpace
// Get obs_height, the observation height
  std::vector<float> obs_height;
  data_.get(parameters_.obs_height, obs_height);

// Get GeoVaLs
  const ufo::GeoVaLs * gvals = data_.getGeoVaLs();

// Setup interpolation of height-dependent thresholds
// N.B. inputs to interp must be double precision
  ufo::PiecewiseLinearInterpolation interp_thresholds(coord_vals, thresholds);

// Loop through locations
  for (size_t iloc=0; iloc < nlocs; ++iloc) {
    // interpolate threshold values to observation height
    float bg_threshold = interp_thresholds(obs_height[iloc]);

    // Vectors storing GeoVaL column for current location.
    std::vector <double> model_profile_column;
    std::vector <double> model_vcoord_column;
    model_profile_column.assign(gvals->nlevs(model_profile_name), 0.0);
    model_vcoord_column.assign(gvals->nlevs(model_vcoord_name), 0.0);
    // Get GeoVaLs at the specified location.
    gvals->getAtLocation(model_profile_column, model_profile_name, iloc);
    gvals->getAtLocation(model_vcoord_column, model_vcoord_name, iloc);

    // interpolate model profile values to observation height
    ufo::PiecewiseLinearInterpolation interp_model(model_vcoord_column, model_profile_column);
    float bg_model = interp_model(obs_height[iloc]);

    // Apply threshold
    if (apply[iloc]) {
      // check to see if one of the compared values is missing
      if (bg_model == missing || bg_threshold == missing) {
        for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
          flagged[jv][iloc] = true;
        }
      } else {
      // Check if model value is outside threshold and set flag
        for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
          if (parameters_.threshold_type == ThresholdType::MIN) {
            if (bg_model < (bg_threshold)) {
              flagged[jv][iloc] = true;
            }
          }
          if (parameters_.threshold_type == ThresholdType::MAX) {
            if (bg_model > (bg_threshold)) {
              flagged[jv][iloc] = true;
            }
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void ModelObThreshold::print(std::ostream & os) const {
  os << "ModelObThreshold: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

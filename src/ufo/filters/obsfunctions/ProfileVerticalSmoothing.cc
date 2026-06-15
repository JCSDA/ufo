/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------
 * Apply vertical smoothing to a set of observations in a profile.
 * -----------------------------------------------------------------------------
 */
#include <Eigen/Dense>
#include <math.h>

#include <algorithm>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ProfileVerticalSmoothing.h"
#include "ufo/filters/Variable.h"

#include "eckit/exception/Exceptions.h"

namespace ufo {

static ObsFunctionMaker<ProfileVerticalSmoothing>
    makerProfileVerticalSmoothing_("ProfileVerticalSmoothing");

/* -----------------------------------------------------------------------------
 * Get the required input parameters from the options, and define the input
 * variables as being required for this function.
 * -----------------------------------------------------------------------------
 */
ProfileVerticalSmoothing::ProfileVerticalSmoothing(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "ProfileVerticalSmoothing constructor" << std::endl;
  // Get input variables from options
  parameters_.deserialize(conf);

  if (parameters_.filterWidth.value() == boost::none ||
      parameters_.filterHeight.value() == boost::none) {
    filterWidths_ = {100, 100};
    filterHeights_ = {0, 60000};
  } else {
    // Get the list of filter widths and heights from options
    filterWidths_ = parseCommaSeparatedFloats(parameters_.filterWidth.value().value());
    filterHeights_ = parseCommaSeparatedFloats(parameters_.filterHeight.value().value());

    if (filterWidths_.size() != filterHeights_.size()) {
      throw eckit::BadParameter("ProfileVerticalSmoothing: filterWidth and filterHeight "
                                "must have the same number of entries", Here());
    }
  }

  polyOrder_ = parameters_.polynomialOrder.value();
  if (polyOrder_ < 1) {
    throw eckit::BadParameter("ProfileVerticalSmoothing: polynomialOrder must be positive",
                              Here());
  }

  invars_ += Variable(parameters_.varApply.value());
  invars_ += Variable(parameters_.heightVariable.value());
}

// -----------------------------------------------------------------------------

ProfileVerticalSmoothing::~ProfileVerticalSmoothing() {
  oops::Log::trace() << "ProfileVerticalSmoothing destructor" << std::endl;
}

/* -----------------------------------------------------------------------------
 * Calculate a smoothing of vertical profiles of observations using a local
 * polynomial regression and write the smoothed values back to the observations
 * array.
 * -----------------------------------------------------------------------------
 */
void ProfileVerticalSmoothing::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "ProfileVerticalSmoothing compute start" << std::endl;
  const ioda::ObsSpace & obsdb = in.obsspace();
  const oops::ObsVariables observed = obsdb.obsvariables();

  // Get the record numbers from the observation data.  These will be used to identify
  // which observations belong to which profile.
  const std::vector<size_t> & record_numbers = obsdb.recidx_all_recnums();

  // Get the observed values to smooth
  std::vector<float> obsVar;
  in.get(parameters_.varApply.value(), obsVar);

  // Get the vertical coordinate values
  std::vector<float> obsHeight;
  in.get(parameters_.heightVariable.value(), obsHeight);

  // Loop over the unique profiles
  for (size_t iProfile : record_numbers) {
    const std::vector<size_t> & obs_numbers = obsdb.recidx_vector(iProfile);

    // Construct the profile of observations to be smoothing, whilst masking
    // out any missing values from the profile
    std::vector<float> obsProfile;
    std::vector<float> heightProfile;
    std::vector<size_t> validIndices;
    for (size_t iob : obs_numbers) {
      if (obsVar[iob] != util::missingValue<float>() &&
          obsHeight[iob] != util::missingValue<float>()) {
        obsProfile.push_back(obsVar[iob]);
        heightProfile.push_back(obsHeight[iob]);
        validIndices.push_back(iob);
      }
    }

    if (obsProfile.size() == 0) {
      continue;  // No valid observations in this profile
    }

    // Loop over the observations in the profile
    std::vector<float> smoothedObs(heightProfile.size(), 0.0f);
    for (int iob = 0; iob < heightProfile.size(); ++iob) {
      // Get the filter width for this point in the profile, based on its height and
      // the filter widths and heights specified in the options.
      float filterWidth;
      filterWidth = getFilterWidth(filterHeights_, filterWidths_, heightProfile[iob]);

      // Normalize the heights by the filter width, centring them on the mean height
      // in the profile.
      double meanHeight = 0;
      for (double height : heightProfile) {
        meanHeight += height;
      }
      meanHeight = meanHeight / heightProfile.size();
      std::vector<double> normalisedHeights(heightProfile.size());
      std::transform(heightProfile.begin(), heightProfile.end(), normalisedHeights.begin(),
                      [=](double x) { return (x - meanHeight) / filterWidth; });

      // Calculate weights using Gaussian kernel
      std::vector<double> weights(heightProfile.size());
      std::transform(normalisedHeights.begin(), normalisedHeights.end(), weights.begin(),
                      [=](double x) { return std::sqrt(std::exp(-(x - normalisedHeights[iob]) *
                                    (x - normalisedHeights[iob]) / 2.0)); });

      // Build the matrix of the weighted polynomial terms for the observations in the profile,
      // up to the specified polynomial order.
      Eigen::MatrixXd X(heightProfile.size(), polyOrder_ + 1);
      for (size_t k = 0; k < heightProfile.size(); ++k) {
        for (size_t p = 0; p <= polyOrder_; ++p) {
          X(k, p) = weights[k] * std::pow(normalisedHeights[k], p);
        }
      }

      // Calculate the matrix to be solved, and find its rank to determine how many
      // polynomial coefficients can be solved for given the heights of the
      // observations in the profile and the specified polynomial order.
      Eigen::MatrixXd matrix = X.transpose() * X;
      int rank = matrix.fullPivLu().rank();
      if (rank > 0) {
        // Calculate the right-hand side of the matrix equation
        Eigen::VectorXd coeffs(rank);
        Eigen::VectorXd weightedObsProfile(heightProfile.size());
        for (size_t k = 0; k < heightProfile.size(); ++k) {
          weightedObsProfile(k) = weights[k] * obsProfile[k];
        }
        Eigen::MatrixXd RHS = X.block(0, 0, X.rows(), rank).transpose() * weightedObsProfile;

        // Calculate the polynomial coefficients by solving the matrix equation,
        // to an order determined by the rank of the matrix.
        coeffs = matrix.block(0, 0, rank, rank).fullPivLu().solve(RHS);

        // Evaluate polynomial at normalized height[iob] to get smoothed value for
        // this observation, using the number of coefficients determined by the
        // rank of the matrix.
        for (int p = 0; p < rank; ++p) {
          smoothedObs[iob] += coeffs[p] * std::pow(normalisedHeights[iob], p);
        }
      } else {
        // If the matrix is rank deficient, set the smoothed value to be the original
        // observation value
        smoothedObs[iob] = obsProfile[iob];
      }

      // Print the original and smoothed values for the first profile as useful information
      if (iProfile == 1) {
        oops::Log::trace() << "Profile " << iProfile << " original values: " << obsProfile[iob]
                           << "  smoothed values: " << smoothedObs[iob]
                           << std::endl;
      }
    }

    // Load the smoothed values back into the output array, ensuring that any observations
    // that were missing in the input profile are set to missing in the output profile
    for (size_t iob : obs_numbers) {
      out[0][iob] = util::missingValue<float>();
    }
    for (size_t iob = 0; iob < validIndices.size(); ++iob) {
      out[0][validIndices[iob]] = smoothedObs[iob];
    }
  }

  oops::Log::trace() << "ProfileVerticalSmoothing compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & ProfileVerticalSmoothing::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------
// Helper function to get filter width at a given height
// -----------------------------------------------------------------------------
float ProfileVerticalSmoothing::getFilterWidth(const std::vector<float>& filterHeights,
                                               const std::vector<float>& filterWidths,
                                               const float targetHeight) const {
  // No interpolation needed if only one filter width and height are provided
  if (filterHeights.size() == 1) {
    return filterWidths[0];
  }

  // Deal with out-of-bounds heights
  if (targetHeight < filterHeights.front()) {
    return filterWidths.front();
  }
  if (targetHeight > filterHeights.back()) {
    return filterWidths.back();
  }

  // Locate the first element greater than `height`
  size_t it = std::distance(filterHeights.begin(), std::lower_bound(
      filterHeights.begin(), filterHeights.end(), targetHeight));

  // Calculate the fraction for linear interpolation between the two nearest filter widths
  float fraction = (targetHeight - filterHeights[it-1]) / (filterHeights[it] - filterHeights[it-1]);

  // Return the linearly interpolated filter width
  return filterWidths[it-1] + fraction * (filterWidths[it] - filterWidths[it-1]);
}

// -----------------------------------------------------------------------------
// Parse a comma-separated string into a vector of floats
// -----------------------------------------------------------------------------
std::vector<float> ProfileVerticalSmoothing::parseCommaSeparatedFloats(
  const std::string & inputString) const {
    std::vector<float> result;
    std::stringstream ss(inputString);
    float thisElement{};

    while (ss >> thisElement)
    {
        result.push_back(thisElement);

        if (ss.peek() == ',')
            ss.ignore();
    }

    return result;
}

}  // namespace ufo

/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#include <algorithm>
#include <set>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/GeoCloudCreateCloudColumn.h"
#include "ufo/filters/Variable.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<GeoCloudCreateCloudColumn>
        makerGeoCloudCreateCloudColumn_("GeoCloudCreateCloudColumn");

// -----------------------------------------------------------------------------

GeoCloudCreateCloudColumn::GeoCloudCreateCloudColumn(
        const eckit::LocalConfiguration & conf): invars_(), channels_() {
  oops::Log::trace() << "GeoCloudCreateCloudColumn constructor" << std::endl;
  // Get options from argument
  options_.deserialize(conf);
  // Get channels used
  const std::set<int> chanset = oops::parseIntSet(options_.channels.value());
  channels_.assign(chanset.begin(), chanset.end());
  ASSERT(channels_.size() > 0);
  // Required cloud top pressure
  invars_ += Variable(options_.cloudTopPressure.value());
  // Required cloud amount
  invars_ += Variable(options_.cloudAmount.value());
  // Required cloud optical thickness
  invars_ += Variable(options_.cloudOpticalThickness.value());
  // Required spatial coherence test flag
  invars_ += Variable(options_.spatialCoherenceTest.value());
  // Required model height above mean sea level, air pressure, and surface altitude
  invars_ += Variable("GeoVaLs/height_above_mean_sea_level");
  invars_ += Variable("GeoVaLs/air_pressure");
  invars_ += Variable("GeoVaLs/height_above_mean_sea_level_at_surface");
}

// -----------------------------------------------------------------------------

void GeoCloudCreateCloudColumn::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "GeoCloudCreateCloudColumn compute start" << std::endl;
  const size_t nlocs = in.nlocs();
  const size_t nlevs = in.nlevs(Variable("GeoVaLs/height_above_mean_sea_level"));
  const GeoVaLs * const gv(in.getGeoVaLs());

  // Cloud top pressure
  std::vector<float> cloudTopPressure(nlocs);
  // Cloud amount
  std::vector<float> cloudAmount(nlocs);
  // Cloud optical thickness
  std::vector<float> cloudOpticalThickness(nlocs);
  // Spatial coherence flag
  std::vector<int> spatialCoherenceTest(nlocs);
  // Error Clear
  const float errorClear = options_.errorClear.value();
  // Error Overcast
  const float errorOvercast = options_.errorOvercast.value();
  // Error for clear columns
  const float errorClearClear = options_.errorClearClear.value();
  // Observation error multiplier
  const float obsErrorMultiplier = options_.obsErrorMultiplier.value();
  // Altitude limit for clear obs
  const float altitudeLimitClear = options_.altitudeLimitClear.value();
  // Altitude limit for cloudy obs
  const float altitudeLimitCloudy = options_.altitudeLimitCloudy.value();
  // Tau scaling denominator
  const float tauScaleDenom = options_.tauScaleDenom.value();
  // Tau scaling offset
  const float tauScaleOffset = options_.tauScaleOffset.value();
  // NDepth
  const std::vector<int> nDepth = options_.nDepth.value();
  if (nDepth.size() != nlevs) {
    throw eckit::UserError(
    "Size of nDepth array should be the same as the number of levels.",
    Here());
  }
  // Maximum number of layers for calculating COT
  const int maxZDepthLayers = options_.maxZDepthLayers.value();
  // Flag to use cloud optical thickness for cloud depth calculation
  const bool useCOT = options_.useCOT.value();
  // Initialise optional parameters
  const std::vector<float> tauZDepthHtLow = options_.tauZDepthHtLow.value();
  const std::vector<float> tauZDepthHtUpp = options_.tauZDepthHtUpp.value();
  const std::vector<float> tauZDepthGrad = options_.tauZDepthGrad.value();
  const std::vector<float> tauZDepthConst = options_.tauZDepthConst.value();
  const std::vector<float> tauMaxZDepth = options_.tauMaxZDepth.value();
  if (useCOT) {
    // Tau lower height limit array
    if (tauZDepthHtLow.size() != maxZDepthLayers) {
      throw eckit::UserError(
      "Cloud depth from cloud optical thickness requested but no tauZDepthHtLow provided.",
       Here());
    }
    // Tau upper height limit array
    if (tauZDepthHtUpp.size() != maxZDepthLayers) {
      throw eckit::UserError(
      "Cloud depth from cloud optical thickness requested but no tauZDepthHtUpp provided.",
      Here());
    }
    // Tau gradient array
    if (tauZDepthGrad.size() != maxZDepthLayers) {
      throw eckit::UserError(
      "Cloud depth from cloud optical thickness requested but no tauZDepthGrad provided.",
      Here());
    }
    // Tau constant array
    if (tauZDepthConst.size() != maxZDepthLayers) {
      throw eckit::UserError(
      "Cloud depth from cloud optical thickness requested but no tauZDepthConst provided.",
      Here());
    }
    // Tau maximum depth array
    if (tauMaxZDepth.size() != maxZDepthLayers) {
      throw eckit::UserError(
      "Cloud depth from cloud optical thickness requested but no tauMaxZDepth provided.",
      Here());
    }
  }

  in.get(Variable(options_.cloudTopPressure.value()), cloudTopPressure);
  in.get(Variable(options_.cloudAmount.value()), cloudAmount);
  in.get(Variable(options_.cloudOpticalThickness.value()), cloudOpticalThickness);
  in.get(Variable(options_.spatialCoherenceTest.value()), spatialCoherenceTest);

  const float missing = util::missingValue<float>();
  const std::vector<std::string> obsErrorName = {"cloudAmount"};
  ioda::ObsDataVector<float> obsErrorOut(in.obsspace(),
                                         oops::ObsVariables(obsErrorName, channels_));
  std::vector<float> pressureDiff(nlevs);
  std::vector<int> cloudTopLevel(nlocs, 0);

  // Populate cloud column for cloudy and clear obs locations
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    std::vector<double> height(nlevs);
    gv->getAtLocation(height, oops::Variable{"height_above_mean_sea_level"}, iloc);
    std::vector<double> PLevels(nlevs);
    gv->getAtLocation(PLevels, oops::Variable{"air_pressure"}, iloc);
    std::vector<double> surfaceAltitude(1);
    gv->getAtLocation(surfaceAltitude,
                      oops::Variable{"height_above_mean_sea_level_at_surface"}, iloc);
    size_t ndep = 1;  // Initialise number of cloud depth levels
    // Find model level closest to measured cloud top pressure
    if (spatialCoherenceTest[iloc] < 2) {
      for (size_t ilev = 0; ilev < nlevs; ++ilev) {
        pressureDiff[ilev] = std::abs(PLevels[ilev] - cloudTopPressure[iloc]);
      }
      std::vector<float>::iterator modelCloudTopPressure =
                                   std::min_element(pressureDiff.begin(), pressureDiff.end());
      cloudTopLevel[iloc] = std::distance(pressureDiff.begin(), modelCloudTopPressure);
    }
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      // Clear obs where spatialCoherenceTest >= 2, else it is cloudy
      if (spatialCoherenceTest[iloc] >= 2) {
        // Intially set all levels cloud fraction = 0.0
        out[ilev][iloc] = 0.0f;
        obsErrorOut[ilev][iloc] = errorClearClear * obsErrorMultiplier;
        // Explicitly set levels to missing:
        // - closest to surface
        // - levels above altitude limit for clear obs
        if (ilev == nlevs - 1 || PLevels[ilev] < altitudeLimitClear) {
          out[ilev][iloc] = missing;
          obsErrorOut[ilev][iloc] = missing;
        }
      // Cloudy obs:
      } else {
        // Check cloud top level is set and is within range of max number of levels
        if (cloudTopLevel[iloc] != 0 && cloudTopLevel[iloc] <= nlevs - 1) {
          // Ensure not to populate level nearest surface (ilev == nlev - 1)
          if (cloudTopLevel[iloc] < nlevs - 1) {
            // Populate levels above cloud top with 0 values
            // unless level is above defined altitude limit, then set missing
            if (ilev < cloudTopLevel[iloc]) {
              if (PLevels[ilev] < altitudeLimitCloudy) {
                out[ilev][iloc] = missing;
                obsErrorOut[ilev][iloc] = missing;
              } else {
                out[ilev][iloc] = 0.0f;
                obsErrorOut[ilev][iloc] = errorClear * obsErrorMultiplier;
              }
            // Populate cloud top level with value in cloudAmount
            } else if (ilev == cloudTopLevel[iloc]) {
              // Still possible for cloud top pressure to be above altitude limit,
              // so use same check here as above
              if (PLevels[ilev] > altitudeLimitCloudy) {
                out[ilev][iloc] = cloudAmount[iloc];
                obsErrorOut[ilev][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                          + errorOvercast * out[ilev][iloc]) * obsErrorMultiplier;
              }
              // Deepen cloud if required, if so calculate how many levels to deepen
              if (cloudTopLevel[iloc] < nlevs - 1 && (nDepth[ilev] > 1 ||
                  cloudOpticalThickness[iloc] != missing)) {
                // Use COT as first preference, if available
                if (cloudOpticalThickness[iloc] != missing && useCOT) {
                  // Scale input COT value back to true COT value
                  // Tau would be calculated as:
                  //   std::exp((cloudOpticalThickness[iloc] / tauScaleDenom) + tauScaleOffset)
                  // However, all that is needed for calculating zDepth is the exponent.
                  const float tau = (cloudOpticalThickness[iloc] / tauScaleDenom)
                              + tauScaleOffset;
                  // Calculate height of cloud top level above surface (model orography)
                  float heightAboveOrography = height[ilev] - surfaceAltitude[0];
                  // Find which height "range" heightAboveOrography falls in
                  // to apply the correct coefficients
                  size_t iRange = 0;
                  while (tauZDepthHtLow[iRange] >= 0.0 &&
                         iRange <= maxZDepthLayers &&
                         heightAboveOrography > tauZDepthHtUpp[iRange]) {
                    iRange++;
                    if (iRange > maxZDepthLayers) {
                      break;
                    }
                  }
                  // Geometrical thickness calculated by linear relationship
                  float zDepth = tauZDepthGrad[iRange] * tau + tauZDepthConst[iRange];
                  // Ensure within allowed range from zero to maximum
                  zDepth = std::max(std::min(zDepth, tauMaxZDepth[iRange]), 0.0f);
                  // Only consider larger values of ndep if zDepth is at least as
                  // far down as next-but-one model level
                  float zDelta = height[ilev] - height[ilev + 2];
                  while (zDelta < zDepth && ndep < (nlevs - cloudTopLevel[iloc])) {
                    ndep++;
                    if (ndep < (nlevs - cloudTopLevel[iloc] - 1)) {
                      zDelta = height[ilev] - height[ilev + ndep + 1];
                    } else {
                      zDelta = height[ilev];
                    }
                  }
                  // Check not above highest height range layer
                  if (tauZDepthHtLow[iRange] < 0.0f) {
                    ndep = 1;
                  }
                // Otherwise use input NDepth array as default
                } else {
                  ndep = nDepth[ilev];
                }
              }
            } else if (ilev > cloudTopLevel[iloc]) {
              // Deepen cloud if ndep value found earlier is greater than 1
              // but only if ilev is between limits set by cloud top and ndep
              if (ndep > 1 && ilev < (cloudTopLevel[iloc] + ndep)) {
                if (PLevels[ilev] > altitudeLimitCloudy) {
                  out[ilev][iloc] = cloudAmount[iloc];
                  obsErrorOut[ilev][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                            + errorOvercast * out[ilev][iloc]) * obsErrorMultiplier;
                }
              // Otherwise, for levels below cloud top, populate with missing
              } else {
                out[ilev][iloc] = missing;
                obsErrorOut[ilev][iloc] = missing;
              }
            }
          // Explicitly deal with situations where cloud is set on level nearest surface
          } else if (cloudTopLevel[iloc] == nlevs - 1) {
            out[ilev][iloc] = missing;
            obsErrorOut[ilev][iloc] = missing;
          }
          // Explicitly set level closest to surface = missing
          if (ilev == nlevs - 1) {
            out[ilev][iloc] = missing;
            obsErrorOut[ilev][iloc] = missing;
          }
          // Explicitly set levels with pressure above altitude limit = missing
          // Mostly to catch obs where cloud top is put on model level nearest surface
          if (PLevels[ilev] < altitudeLimitCloudy) {
            out[ilev][iloc] = missing;
            obsErrorOut[ilev][iloc] = missing;
          }
        }
      }
    }
  }

  obsErrorOut.save(options_.outputGroup.value());
  oops::Log::trace() << "GeoCloudCreateCloudColumn compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & GeoCloudCreateCloudColumn::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo


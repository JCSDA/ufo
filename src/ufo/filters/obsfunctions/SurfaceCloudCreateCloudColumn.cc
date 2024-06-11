/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#include <set>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/SurfaceCloudCreateCloudColumn.h"
#include "ufo/filters/Variable.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<SurfaceCloudCreateCloudColumn>
        makerSurfaceCloudCreateCloudColumn_("SurfaceCloudCreateCloudColumn");

// -----------------------------------------------------------------------------

SurfaceCloudCreateCloudColumn::SurfaceCloudCreateCloudColumn(
        const eckit::LocalConfiguration & conf): invars_(), channels_() {
  // Get options from argument
  options_.deserialize(conf);
  // Get channels used
  std::set<int> chanset = oops::parseIntSet(options_.channels.value());
  channels_.assign(chanset.begin(), chanset.end());
  ASSERT(channels_.size() > 0);
  // Required cloud fraction at base
  invars_ += Variable(options_.cloudFractionAtBase.value());
  // Required cloud fraction at base for METAR layers 2 and 3
  invars_ += Variable(options_.cloudFractionAtBaseMetar2.value());
  invars_ += Variable(options_.cloudFractionAtBaseMetar3.value());
  // Required cloud fraction at base for deep convective cloud (3 layers)
  invars_ += Variable(options_.cloudFractionAtBaseDCC1.value());
  invars_ += Variable(options_.cloudFractionAtBaseDCC2.value());
  invars_ += Variable(options_.cloudFractionAtBaseDCC3.value());
  // Required model level cloud base height
  invars_ += Variable(options_.modelLevelCloudBaseHeight.value());
  // Required model level cloud base height for METAR layers 2 and 3
  invars_ += Variable(options_.modelLevelCloudBaseHeightMetar2.value());
  invars_ += Variable(options_.modelLevelCloudBaseHeightMetar3.value());
  // Required model level cloud base height for deep convective cloud (3 layers)
  invars_ += Variable(options_.modelLevelCloudBaseHeightDCC1.value());
  invars_ += Variable(options_.modelLevelCloudBaseHeightDCC2.value());
  invars_ += Variable(options_.modelLevelCloudBaseHeightDCC3.value());
  // Required model geopotenial height
  invars_ += Variable("GeoVaLs/height");
}

// -----------------------------------------------------------------------------

void SurfaceCloudCreateCloudColumn::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const size_t nlevs = in.nlevs(Variable("GeoVaLs/height"));
  const GeoVaLs * const gv(in.getGeoVaLs());

  // Cloud fraction at base
  std::vector<float> cloudFractionAtBase(nlocs);
  std::vector<float> cloudFractionAtBaseMetar2(nlocs);
  std::vector<float> cloudFractionAtBaseMetar3(nlocs);
  std::vector<float> cloudFractionAtBaseDCC1(nlocs);
  std::vector<float> cloudFractionAtBaseDCC2(nlocs);
  std::vector<float> cloudFractionAtBaseDCC3(nlocs);
  // Model level cloud base height
  std::vector<float> modelLevelCBH(nlocs);
  std::vector<float> modelLevelCBHMetar2(nlocs);
  std::vector<float> modelLevelCBHMetar3(nlocs);
  std::vector<float> modelLevelCBHDCC1(nlocs);
  std::vector<float> modelLevelCBHDCC2(nlocs);
  std::vector<float> modelLevelCBHDCC3(nlocs);
  // Error Clear
  const float errorClear = options_.errorClear.value();
  // Error Overcast
  const float errorOvercast = options_.errorOvercast.value();
  // Observation error multiplier
  const float obsErrorMultiplier = options_.obsErrorMultiplier.value();
  // Altitude limit
  const float altitudeLimit = options_.altitudeLimit.value();
  // numNDepthLevels
  const int numNDepthLevels1 = options_.numNDepthLevels1.value();
  const int numNDepthLevels2 = options_.numNDepthLevels2.value();
  const int numNDepthLevels3 = options_.numNDepthLevels3.value();
  const int numNDepthLevels4 = options_.numNDepthLevels4.value();

  // Set up NDepth array for adding cloud depth
  // - These values come from an analysis of the model
  //   "climatology" of cloud thickness.
  // - A value of 0 is given to the level nearest the surface
  //   so that it is ignored.
  // - Levels above that are assigned with thicknesses in terms of
  //   model level, which can be set as input parameters.
  // - Defaults: 33 * 1, 23 * 2, 5 * 3, 8 * 4 (plus 1 * 0 for surface level)
  // - If model levels are added above current model top (70 levels)
  //   they should be assigned a thickness of 1.
  std::vector<int> NDepth;
  NDepth.insert(NDepth.end(), numNDepthLevels1, 1);
  NDepth.insert(NDepth.end(), numNDepthLevels2, 2);
  NDepth.insert(NDepth.end(), numNDepthLevels3, 3);
  NDepth.insert(NDepth.end(), numNDepthLevels4, 4);
  NDepth.insert(NDepth.end(), 1, 0);

  in.get(Variable(options_.cloudFractionAtBase.value()), cloudFractionAtBase);
  in.get(Variable(options_.cloudFractionAtBaseMetar2.value()), cloudFractionAtBaseMetar2);
  in.get(Variable(options_.cloudFractionAtBaseMetar3.value()), cloudFractionAtBaseMetar3);
  in.get(Variable(options_.cloudFractionAtBaseDCC1.value()), cloudFractionAtBaseDCC1);
  in.get(Variable(options_.cloudFractionAtBaseDCC2.value()), cloudFractionAtBaseDCC2);
  in.get(Variable(options_.cloudFractionAtBaseDCC3.value()), cloudFractionAtBaseDCC3);
  in.get(Variable(options_.modelLevelCloudBaseHeight.value()), modelLevelCBH);
  in.get(Variable(options_.modelLevelCloudBaseHeightMetar2.value()), modelLevelCBHMetar2);
  in.get(Variable(options_.modelLevelCloudBaseHeightMetar3.value()), modelLevelCBHMetar3);
  in.get(Variable(options_.modelLevelCloudBaseHeightDCC1.value()), modelLevelCBHDCC1);
  in.get(Variable(options_.modelLevelCloudBaseHeightDCC2.value()), modelLevelCBHDCC2);
  in.get(Variable(options_.modelLevelCloudBaseHeightDCC3.value()), modelLevelCBHDCC3);

  const float missing = util::missingValue<float>();
  const std::vector<std::string> obsErrorName = {"cloudAmount"};
  ioda::ObsDataVector<float> obsErrorOut(in.obsspace(),
                                         oops::ObsVariables(obsErrorName, channels_));
  const int nlayers = 3;
  std::vector<float> cloudBaseHeight(nlayers);
  float NextCloudLayerUp = missing;

  // Populate cloud column for cloudy and clear obs locations
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    std::fill(cloudBaseHeight.begin(), cloudBaseHeight.end(), missing);
    std::vector<float> Height(nlevs);
    gv->getAtLocation(Height, oops::Variable{"height"}, iloc);
    cloudBaseHeight[0] = modelLevelCBH[iloc];
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      // Now check whether the Height of each model level is above or below
      // the measured cloud base height:
      //  If it's below, populate with 0.0 (no cloud)
      //  If it's above, populate with missing
      //  If it's matching, place the cloud fraction
      if ((Height[ilev] - modelLevelCBH[iloc]) < -1.0e-1) {
        if ((modelLevelCBHDCC1[iloc] != missing && modelLevelCBH[iloc] < modelLevelCBHDCC1[iloc])
            || modelLevelCBHDCC1[iloc] == missing) {
          out[ilev][iloc] = 0.0f;
          obsErrorOut[ilev][iloc] = errorClear * obsErrorMultiplier;
        } else {
          out[ilev][iloc] = missing;
          obsErrorOut[ilev][iloc] = missing;
        }
      } else if (std::abs(Height[ilev] - modelLevelCBH[iloc]) < 1.0e-2) {
        out[ilev][iloc] = cloudFractionAtBase[iloc] / 8.0f;  // convert from oktas
        obsErrorOut[ilev][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                  + errorOvercast * out[ilev][iloc]) * obsErrorMultiplier;
      } else if ((Height[ilev] - modelLevelCBH[iloc]) > 1.0e-1) {
        out[ilev][iloc] = missing;
        obsErrorOut[ilev][iloc] = missing;
      }
    }

    // Deepen the cloud where appropriate
    // An input value of NDepth = 0 is ignored, and so does not suppress cloud at that level
    // Only deepen the cloud if ob is greater than or equal to one okta
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      if (std::abs(modelLevelCBH[iloc] - Height[ilev]) < 1.0e-2 && NDepth[ilev] > 1 &&
          cloudFractionAtBase[iloc] >= 1.0f) {
        const size_t ndep = NDepth[ilev]-1;
        for (size_t znumber = 1; znumber <= ndep; znumber++) {
          if (Height[ilev-znumber] <= altitudeLimit) {
            out[ilev-znumber][iloc] = cloudFractionAtBase[iloc] / 8.0f;
            obsErrorOut[ilev-znumber][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                              + errorOvercast * out[ilev][iloc])
                                              * obsErrorMultiplier;
          }
        }
      }
    }

    // Next add in METAR extra cloud layers, if present
    // N.b. This may overwrite cloud already set, but that
    // is desired behaviour as this is the best information.
    // Second METAR layer
    if (modelLevelCBHMetar2[iloc] != missing && cloudFractionAtBaseMetar2[iloc] != missing) {
      cloudBaseHeight[1] = modelLevelCBHMetar2[iloc];
      for (size_t ilev = 0; ilev < nlevs; ++ilev) {
        // set cloud at this level
        if (std::abs(Height[ilev] - modelLevelCBHMetar2[iloc]) < 1.0e-2 &&
            Height[ilev] <= altitudeLimit) {
          out[ilev][iloc] = cloudFractionAtBaseMetar2[iloc] / 8.0f;  // convert from oktas
          obsErrorOut[ilev][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                    + errorOvercast * out[ilev][iloc]) * obsErrorMultiplier;
        }
        if (std::abs(modelLevelCBHMetar2[iloc] - Height[ilev]) < 1.0e-2 && NDepth[ilev] > 1 &&
           cloudFractionAtBaseMetar2[iloc] >= 1.0f) {  // deepen
          const size_t ndep = NDepth[ilev]-1;
          for (size_t znumber = 1; znumber <= ndep; znumber++) {
            if (Height[ilev-znumber] <= altitudeLimit) {
              out[ilev-znumber][iloc] = cloudFractionAtBaseMetar2[iloc] / 8.0f;
              obsErrorOut[ilev-znumber][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                                + errorOvercast * out[ilev][iloc])
                                                * obsErrorMultiplier;
            }
          }
        }
      }
    }
    // Third METAR layer
    if (modelLevelCBHMetar3[iloc] != missing && cloudFractionAtBaseMetar3[iloc] != missing) {
      cloudBaseHeight[2] = modelLevelCBHMetar3[iloc];
      for (size_t ilev = 0; ilev < nlevs; ++ilev) {
        // set cloud at this level
        if (std::abs(Height[ilev] - modelLevelCBHMetar3[iloc]) < 1.0e-2 &&
            Height[ilev] <= altitudeLimit) {
          out[ilev][iloc] = cloudFractionAtBaseMetar3[iloc] / 8.0f;  // convert from oktas
          obsErrorOut[ilev][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                    + errorOvercast * out[ilev][iloc]) * obsErrorMultiplier;
        }
        if (std::abs(modelLevelCBHMetar3[iloc] - Height[ilev]) < 1.0e-2 && NDepth[ilev] > 1 &&
           cloudFractionAtBaseMetar3[iloc] >= 1.0f) {  // deepen
          const size_t ndep = NDepth[ilev]-1;
          for (size_t znumber = 1; znumber <= ndep; znumber++) {
            if (Height[ilev-znumber] <= altitudeLimit) {
              out[ilev-znumber][iloc] = cloudFractionAtBaseMetar3[iloc] / 8.0f;
              obsErrorOut[ilev-znumber][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                                + errorOvercast * out[ilev][iloc])
                                                * obsErrorMultiplier;
            }
          }
        }
      }
    }
    // Now add the deep convective cloud (DCC) if any exists
    // and only deepen if it will not interfere with a higher normal cloud layer
    // First DCC layer
    if (modelLevelCBHDCC1[iloc] != missing && cloudFractionAtBaseDCC1[iloc] != missing) {
      for (size_t ilev = 0; ilev < nlevs; ++ilev) {
        if ((Height[ilev] - modelLevelCBHDCC1[iloc]) < -1.0e-1 && Height[ilev] <= altitudeLimit) {
          if ((modelLevelCBH[iloc] != missing && modelLevelCBH[iloc] >= modelLevelCBHDCC1[iloc])
             || modelLevelCBH[iloc] == missing) {
            out[ilev][iloc] = 0.0f;
            obsErrorOut[ilev][iloc] = errorClear * obsErrorMultiplier;
          }
        // set cloud at this level
        } else if (std::abs(Height[ilev] - modelLevelCBHDCC1[iloc]) < 1.0e-2 &&
                   Height[ilev] <= altitudeLimit) {
          out[ilev][iloc] = cloudFractionAtBaseDCC1[iloc] / 8.0f;  // convert from oktas
          obsErrorOut[ilev][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                    + errorOvercast * out[ilev][iloc]) * obsErrorMultiplier;
        }
        if (std::abs(modelLevelCBHDCC1[iloc] - Height[ilev]) < 1.0e-2 && NDepth[ilev] > 1 &&
           cloudFractionAtBaseDCC1[iloc] >= 1.0f) {  // deepen
          for (size_t lnumber = 0; lnumber <= nlayers; lnumber++) {
            if (cloudBaseHeight[lnumber] == missing
               || (lnumber == 3 && Height[ilev] > cloudBaseHeight[lnumber])) {
              NextCloudLayerUp = altitudeLimit + 1.0f;
            } else if (Height[ilev] < cloudBaseHeight[lnumber]) {
              NextCloudLayerUp = cloudBaseHeight[lnumber];
              break;
            }
          }
          const size_t ndep = NDepth[ilev]-1;
          for (size_t znumber = 1; znumber <= ndep; znumber++) {
            if (Height[ilev-znumber] <= altitudeLimit &&
               Height[ilev-znumber] < NextCloudLayerUp) {
              out[ilev-znumber][iloc] = cloudFractionAtBaseDCC1[iloc] / 8.0f;
              obsErrorOut[ilev-znumber][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                                + errorOvercast * out[ilev][iloc])
                                                * obsErrorMultiplier;
            }
          }
        }
      }
    }
    // Second DCC layer
    if (modelLevelCBHDCC2[iloc] != missing && cloudFractionAtBaseDCC2[iloc] != missing) {
      NextCloudLayerUp = missing;
      for (size_t ilev = 0; ilev < nlevs; ++ilev) {
        // set cloud at this level
        if (std::abs(Height[ilev] - modelLevelCBHDCC2[iloc]) < 1.0e-2
           && Height[ilev] <= altitudeLimit) {
          out[ilev][iloc] = cloudFractionAtBaseDCC2[iloc] / 8.0f;  // convert from oktas
          obsErrorOut[ilev][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                    + errorOvercast * out[ilev][iloc]) * obsErrorMultiplier;
        }
        if (std::abs(modelLevelCBHDCC2[iloc] - Height[ilev]) < 1.0e-2 && NDepth[ilev] > 1 &&
           cloudFractionAtBaseDCC2[iloc] >= 1.0f) {  // deepen
          for (size_t lnumber = 0; lnumber <= nlayers; lnumber++) {
            if (cloudBaseHeight[lnumber] == missing
               || (lnumber == 3 && Height[ilev] > cloudBaseHeight[lnumber])) {
              NextCloudLayerUp = altitudeLimit + 1.0f;
            } else if (Height[ilev] < cloudBaseHeight[lnumber]) {
              NextCloudLayerUp = cloudBaseHeight[lnumber];
              break;
            }
          }
          const size_t ndep = NDepth[ilev]-1;
          for (size_t znumber = 1; znumber <= ndep; znumber++) {
            if (Height[ilev-znumber] <= altitudeLimit &&
               Height[ilev-znumber] < NextCloudLayerUp) {
              out[ilev-znumber][iloc] = cloudFractionAtBaseDCC2[iloc] / 8.0f;
              obsErrorOut[ilev-znumber][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                                + errorOvercast * out[ilev][iloc])
                                                * obsErrorMultiplier;
            }
          }
        }
      }
    }
    // Third DCC layer
    if (modelLevelCBHDCC3[iloc] != missing && cloudFractionAtBaseDCC3[iloc] != missing) {
      NextCloudLayerUp = missing;
      for (size_t ilev = 0; ilev < nlevs; ++ilev) {
        // set cloud at this level
        if (std::abs(Height[ilev] - modelLevelCBHDCC3[iloc]) < 1.0e-2 &&
           Height[ilev] <= altitudeLimit) {
          out[ilev][iloc] = cloudFractionAtBaseDCC3[iloc] / 8.0f;  // convert from oktas
          obsErrorOut[ilev][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                    + errorOvercast * out[ilev][iloc]) * obsErrorMultiplier;
        }
        if (std::abs(modelLevelCBHDCC3[iloc] - Height[ilev]) < 1.0e-2 && NDepth[ilev] > 1 &&
           cloudFractionAtBaseDCC3[iloc] >= 1.0f) {  // deepen
          for (size_t lnumber = 0; lnumber <= nlayers; lnumber++) {
            if (cloudBaseHeight[lnumber] == missing
               || (lnumber == 3 && Height[ilev] > cloudBaseHeight[lnumber])) {
              NextCloudLayerUp = altitudeLimit + 1.0;
            } else if (Height[ilev] < cloudBaseHeight[lnumber]) {
              NextCloudLayerUp = cloudBaseHeight[lnumber];
              break;
            }
          }
          const size_t ndep = NDepth[ilev]-1;
          for (size_t znumber = 1; znumber <= ndep; znumber++) {
            if (Height[ilev-znumber] <= altitudeLimit &&
               Height[ilev-znumber] < NextCloudLayerUp) {
              out[ilev-znumber][iloc] = cloudFractionAtBaseDCC3[iloc] / 8.0f;
              obsErrorOut[ilev-znumber][iloc] = (errorClear * (1.0f - out[ilev][iloc])
                                                + errorOvercast * out[ilev][iloc])
                                                * obsErrorMultiplier;
            }
          }
        }
      }
    }
  }

  obsErrorOut.save(options_.outputGroup.value());
}

// -----------------------------------------------------------------------------

const ufo::Variables & SurfaceCloudCreateCloudColumn::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

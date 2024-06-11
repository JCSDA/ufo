/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/SurfaceCloudModelLevelCBH.h"
#include "ufo/filters/Variable.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<SurfaceCloudModelLevelCBH>
        makerSurfaceCloudModelLevelCBH_("SurfaceCloudModelLevelCBH");

// -----------------------------------------------------------------------------

SurfaceCloudModelLevelCBH::SurfaceCloudModelLevelCBH(
        const eckit::LocalConfiguration & conf): invars_() {
  // Get options from argument
  options_.deserialize(conf);
  // Required cloud base height above mean sea level
  invars_ += Variable(options_.cloud_base_height.value());
  // Required model geopotenial height
  invars_ += Variable("GeoVaLs/height");
}

// -----------------------------------------------------------------------------

void SurfaceCloudModelLevelCBH::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const size_t nlevs = in.nlevs(Variable("GeoVaLs/height"));
  const GeoVaLs * const gv(in.getGeoVaLs());

  // Cloud Base Height
  std::vector<float> CBH(nlocs);
  // Number of clear levels
  std::vector<int> ZClear(nlocs);

  in.get(Variable(options_.cloud_base_height.value()), CBH);

  const float missing = util::missingValue<float>();

  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    std::vector<float> Height(nlevs);
    gv->getAtLocation(Height, oops::Variable{"height"}, iloc);
    ZClear[iloc] = -1;
    for (size_t ilev = 0; ilev < nlevs; ilev++) {
      if (ilev == nlevs-1 && Height[ilev] >= CBH[iloc]) {
        if (CBH[iloc] == missing) {
          ZClear[iloc] = -1;
        } else {
          ZClear[iloc] = nlevs;
        }
      } else if (Height[ilev] < CBH[iloc]) {
        ZClear[iloc] = ilev;
        break;
      }
    }
    // ZClear is the last level with no cloud
    // GeoVaLs are backwards so ZClear - 1 is first cloudy level
    // (as opposed to ZClear + 1 in OPS)
    if (ZClear[iloc] == -1) {
      out[0][iloc] = missing;
    } else {
      out[0][iloc] = Height[ZClear[iloc] - 1];
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & SurfaceCloudModelLevelCBH::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

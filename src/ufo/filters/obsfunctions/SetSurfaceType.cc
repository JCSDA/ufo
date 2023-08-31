/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SetSurfaceType.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/SurfaceReportConstants.h"

namespace ufo {

  static ObsFunctionMaker<SetSurfaceType> makerSetSurfaceType_("SetSurfaceType");

  SetSurfaceType::SetSurfaceType(const eckit::LocalConfiguration & conf)
    : invars_() {
    // Initialize options
    options_.validateAndDeserialize(conf);

    // Include list of required data from GeoVaLs
    invars_ += Variable("GeoVaLs/ice_area_fraction");
    invars_ += Variable("GeoVaLs/surface_altitude");

    // Include list of required data from ObsSpace
    invars_ += Variable("MetaData/latitude");

    if (options_.UseReportSurface.value()) {
      invars_ += Variable(options_.SurfaceMetaDataName.value());
    }

    if (options_.UseReportElevation.value()) {
      invars_ += Variable("MetaData/heightOfSurface");
    }

    if (options_.UseAAPPSurfaceClass.value()) {
      invars_ += Variable("MetaData/surfaceClassAAPP");
    }

    if (options_.UseSurfaceWaterFraction.value()) {
      invars_ += Variable("MetaData/waterAreaFraction");
    }
  }

  // -----------------------------------------------------------------------------

  SetSurfaceType::~SetSurfaceType() {}

  // -----------------------------------------------------------------------------

  void SetSurfaceType::compute(const ObsFilterData & in,
                               ioda::ObsDataVector<float> & out) const {
    // Get dimension
    const size_t nlocs = in.nlocs();

    std::vector<float> ice_area_frac(nlocs), model_height(nlocs), latitude(nlocs), elevation(nlocs);
    std::vector<int> land_sea(nlocs);
    std::vector<int> surftype(nlocs, options_.SurfaceTypeDefault.value());

    // mandatory variables
    in.get(Variable("GeoVaLs/ice_area_fraction"), ice_area_frac);
    in.get(Variable("MetaData/latitude"), latitude);
    in.get(Variable("GeoVaLs/surface_altitude"), model_height);

    int surftype_land_   = options_.SurfaceTypeLand.value();
    int surftype_sea_    = options_.SurfaceTypeSea.value();
    int surftype_seaice_ = options_.SurfaceTypeSeaIce.value();

    float heighttolerance = options_.HeightTolerance.value();

    // if available and requested, set elevation to ob surface height
    if (options_.UseReportElevation.value()) {
      if (in.has(Variable("MetaData/heightOfSurface"))) {
        in.get(Variable("MetaData/heightOfSurface"), elevation);
      } else {
        oops::Log::warning() << "UseReportElevation is true but MetaData/heightOfSurface"
                             << " not present. Using model data for elevation data\n";
          elevation = model_height;
      }
    } else {  // otherwise
      elevation = model_height;
    }

    // set land_sea mask according to elevation data from ob (if available and requested) and model
    // note that the elevation > 0 comes directly from OPS and so negative surface report elevations
    // require the model height to be non-zero to be recognised as land which is likely but this is
    // probably not logically correct.
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      land_sea[iloc] = (elevation[iloc] > 0.0f ||
                        model_height[iloc] > 0.0f + heighttolerance ||
                        model_height[iloc] < 0.0f - heighttolerance) ?
        surftype_land_ : surftype_sea_;
    }

    // if available and requested, set closest appropriate surface type using reported surface
    if (options_.UseReportSurface.value()) {
      if (in.has(Variable(options_.SurfaceMetaDataName.value()))) {
        std::vector<int> reported_surftype(nlocs);

        // load reported surf type taking into account variable name from the yaml config file
        in.get(Variable(options_.SurfaceMetaDataName.value()), reported_surftype);

        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          if ( reported_surftype[iloc] == AAPP_surftype::sea ||
               reported_surftype[iloc] == BUFR_surftype::ocean ||
               reported_surftype[iloc] == BUFR_surftype::posice ) {
            surftype[iloc] = surftype_sea_;
          } else if ( reported_surftype[iloc] == BUFR_surftype::land ) {
            surftype[iloc] = surftype_land_;
          } else if ( reported_surftype[iloc] == BUFR_surftype::ice ) {
            surftype[iloc] = surftype_seaice_;
          } else if ( reported_surftype[iloc] == BUFR_surftype::coast ||
                      reported_surftype[iloc] == BUFR_surftype::nrcoast) {
            surftype[iloc] = land_sea[iloc];
          }
        }
      } else {
        oops::Log::warning() << "UseReportSurface is true but "
                             << options_.SurfaceMetaDataName.value() << " not present. "
                             << "Using elevation data to determine initial surface type\n";
        surftype = land_sea;
      }
    } else {  // set surface type using land_sea directly
        surftype = land_sea;
    }

    // if available and requested, set closest appropriate surface type using water_fraction
    if (options_.UseSurfaceWaterFraction.value()) {
      if (in.has(Variable("MetaData/waterAreaFraction"))) {
        std::vector<float> water_fraction(nlocs);
        in.get(Variable("MetaData/waterAreaFraction"), water_fraction);

        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          if (water_fraction[iloc] > options_.MinWaterFrac.value()) {
            surftype[iloc] = surftype_sea_;
          } else {
            surftype[iloc] = surftype_land_;
          }
        }
      } else {
        oops::Log::warning() << "UseSurfaceWaterFraction is true but MetaData/waterAreaFraction "
                             << "not present. Ignoring\n";
          }
    }

    // Set sea ice surfaces
    // Only sea spots can be reclassified as ice (land points that may be covered
    // with ice are left as land as we don't have a suitable method of determining
    // where they are). The presence of seaice is determined by enforcing a
    // threshold on the model seaice fraction.
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (surftype[iloc] == surftype_sea_ &&
          ice_area_frac[iloc] >= options_.MinIceFrac.value()) {
        surftype[iloc] = surftype_seaice_;
      }
    }

    // AAPP can provide additional surface type information derived from radiances.
    // This can be used to help identify sea and seaice surfaces correctly,
    // although it can give odd results at very high latitudes

    // if available and requested, set closest appropriate surface type using AAPP surface class
    if (options_.UseAAPPSurfaceClass.value()) {
      if (in.has(Variable("MetaData/surfaceClassAAPP"))) {
        std::vector<int> AAPP_surface_class(nlocs);
        in.get(Variable("MetaData/surfaceClassAAPP"), AAPP_surface_class);

        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          if (AAPP_surface_class[iloc] == AAPP_surfclass::sea) {
          // reclassify surface as sea if prior classification hasn't as long as not 'highland'
            if (surftype[iloc] != surftype_sea_ &&
                elevation[iloc] < options_.HighlandHeight.value()) {
              surftype[iloc] = surftype_sea_;
            }
          } else if (AAPP_surface_class[iloc] >= AAPP_surfclass::newice &&  // AAPP_surface_class
                     AAPP_surface_class[iloc] <= AAPP_surfclass::desert) {  // must be valid
            if (surftype[iloc] == surftype_sea_ &&
                std::abs(latitude[iloc]) >= options_.IceLimitSoft.value() &&
                ice_area_frac[iloc] >= 0.0f) {
              surftype[iloc] = surftype_seaice_;
            }
          }
        }
      } else {
        oops::Log::warning() << "UseAAPPSurfaceClass is true but MetaData/surfaceClassAAPP not "
                             << "present. Ignoring\n";
          }
    }

    // Any sea point south of IceShelfLimit is assumed to be ice
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if (surftype[iloc] == surftype_sea_) {
        if (latitude[iloc] < -1.0f * options_.IceLimitHard.value()) {
          surftype[iloc] = surftype_seaice_;
        }
      }
    }

    // Finally assign surftype to obsfunction output

    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      out[0][iloc] = static_cast <float> (surftype[iloc]);
      }
  }

  // -----------------------------------------------------------------------------

  const ufo::Variables & SetSurfaceType::requiredVariables() const {
    return invars_;
  }

  // -----------------------------------------------------------------------------

}  // namespace ufo

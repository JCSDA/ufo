/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_SURFACEREPORTCONSTANTS_H_
#define UFO_UTILS_SURFACEREPORTCONSTANTS_H_

//-------------------------------------------------------------------------------------------------

namespace ufo {

//-------------------------------------------------------------------------------------------------

/// AAPP surface classification
  // these appear to be shifted by 1 (i.e. 0-8) compared to the original definition (1-9)
  // in NWPSAF-MF-UD-003 (AAPP documentation - data formats)

struct AAPP_surfclass {
  static constexpr int newice    = 0;  // Bare young ice (new ice, no snow)
  static constexpr int dryland   = 1;  // Dry land (dry w/ or w/o sig. vegetation)
  static constexpr int drysnow   = 2;  // Dry snow (snow with water < 2%, over land)
  static constexpr int multiice  = 3;  // Multi-year ice (old ice w/ snow cover)
  static constexpr int sea       = 4;  // Sea (open water, no islands, ice-free)
  static constexpr int wetforest = 5;  // Wet forest (est. forest w/ wet canopy)
  static constexpr int wetland   = 6;  // Wet land (non-forested land w/ wet surface)
  static constexpr int wetsnow   = 7;  // Wet snow (w/ water > 2%, over land/ice)
  static constexpr int desert    = 8;  // Desert
};

struct BUFR_surftype {
  // BUFR code surface definitons (013040 0009)
  static constexpr int land         = 0;
  static constexpr int nrcoast      = 2;
  static constexpr int ice          = 3;
  static constexpr int posice       = 4;  // possible ice
  static constexpr int ocean        = 5;
  static constexpr int coast        = 6;
  static constexpr int inland_water = 7;
  static constexpr int snow_cover   = 8;
  static constexpr int sea_ice      = 9;
  static constexpr int std_water    = 10;
  static constexpr int snow         = 11;
  static constexpr int missing      = 15;
};

struct AAPP_surftype {
  static constexpr int land         = 0;
  static constexpr int sea          = 1;
  static constexpr int coast        = 2;
};

struct RTTOV_surftype {
  // RTTOV surface type constants from rttov_const_mod
  static constexpr int land    = 0;
  static constexpr int sea     = 1;
  static constexpr int seaice  = 2;
};

struct CRTM_surftype {
  static constexpr int invalid = 0;
  static constexpr int land    = 1;
  static constexpr int water   = 2;
  static constexpr int snow    = 4;
  static constexpr int ice     = 8;
};

//--------------------------------------------------------------------------------------------------
}  // namespace ufo

#endif  // UFO_UTILS_SURFACEREPORTCONSTANTS_H_

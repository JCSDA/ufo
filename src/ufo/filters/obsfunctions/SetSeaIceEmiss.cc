/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SetSeaIceEmiss.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

  static ObsFunctionMaker<SetSeaIceEmiss> makerSetSeaIceEmiss_("SetSeaIceEmiss");

  SetSeaIceEmiss::SetSeaIceEmiss(const eckit::LocalConfiguration & conf)
    : invars_() {
    // Initialize options
    options_.validateAndDeserialize(conf);

    // Get channels from options
    std::set<int> channelset = oops::parseIntSet(options_.channelList);
    std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
    ASSERT(channels_.size() > 0);

    // Include list of required data from GeoVaLs/Obs
    // Model ice area fraction [0-1]
    invars_ += Variable("GeoVaLs/ice_area_fraction");

    // Sensor Zenith Angle (degrees)
    invars_ += Variable("MetaData/sensorZenithAngle");
    // AAPP Surface Classification [0-9]
    invars_ += Variable("MetaData/surfaceClassAAPP");
    // RTTOV Surface type [0-2]
    invars_ += Variable("MetaData/surfaceQualifier");
  }

  // -----------------------------------------------------------------------------

  SetSeaIceEmiss::~SetSeaIceEmiss() {}

  // -----------------------------------------------------------------------------

  /// \details From Fast Models for Land Surface Emissivity report [Hewison and English (1999)]:
  /// This is a semi-empirical model that uses Fresnel’s formulae to calculate the specular
  /// reflectivity of a dielectric surface, whose permittivity can be described by a single Debye
  /// relaxation, neglecting the ionic conductivity term, as this is negligible for frequencies
  /// above 20GHz.
  void SetSeaIceEmiss::compute(const ObsFilterData & in,
                               ioda::ObsDataVector<float> & out) const {
    // Get dimensions and initialise constants
    const size_t nlocs = in.nlocs();
    const size_t nchans = channels_.size();
    const int missingValueInt = util::missingValue<int>();

    // constant emissivities defined for different ice types
    static constexpr float emissivity_multiice = 0.84f;
    static constexpr float emissivity_newice = 0.95f;
    static constexpr float emissivity_ice_default = 0.92f;

    // grease ice FASTEM emissivity parameters
    static constexpr float perm_static = 23.7f;   // relative static permittivity
    static constexpr float perm_infinite = 7.7f;  // relative permittivity at infinite frequency
    static constexpr float freqr = 17.3f;         // relaxation frequency

    // define fraction of V and H pol for different channel polarisations
    // index 0: [0.5(V+H)]
    // index 1: [V]
    // index 2: [H]
    std::array<std::array<float, 3>, 3> constexpr zvpol
    {{{0.5f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f}, {0.0f, 1.0f, 0.0f}}};
    std::array<std::array<float, 3>, 3> constexpr zhpol
    {{{0.5f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}}};

    // inputs
    std::vector<float> ice_area_frac(nlocs);
    in.get(Variable("GeoVaLs/ice_area_fraction"), ice_area_frac);

    std::vector<float> satzenith(nlocs);
    in.get(Variable("MetaData/sensorZenithAngle"), satzenith);

    std::vector<int> surftype(nlocs);
    in.get(Variable("MetaData/surfaceQualifier"), surftype);

    std::vector<int> AAPP_surface_class(nlocs);
    in.get(Variable("MetaData/surfaceClassAAPP"), AAPP_surface_class);

    // Get orbit Height (km) from options
    float orbit_height = options_.orbitHeight.value();

    // Get polarisation indices from options
    std::vector<size_t> polindex = options_.polList.value();

    // Get frequencies (GHz) from options
    std::vector<float> freqghz = options_.freqList.value();

    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      if ( surftype[iloc] == RTTOV_surftype::seaice ) {
        // set default ice emissivity
        for (size_t ichan = 0; ichan < nchans; ++ichan) {
          out[ichan][iloc] = emissivity_ice_default;
        }

        if ( ice_area_frac[iloc] >= 0.0f &&
           AAPP_surface_class[iloc] != AAPP_surfclass::sea &&
           AAPP_surface_class[iloc] != missingValueInt) {
          // set ice emissivity to constant depending on AAPP surface classification
          float ice_emiss = (AAPP_surface_class[iloc] == AAPP_surfclass::multiice) ?
            emissivity_multiice : emissivity_newice;

          // calculate geometry
          // satellite zenith angle
          float satzenith_rad = satzenith[iloc] * Constants::deg2rad;

          float czen = std::cos(satzenith_rad);
          float czen2 = czen * czen;
          float szen2 = 1.0f - czen2;

          // view angle to satellite accounting for Earth curvature
          float orbit_height_frac = Constants::mean_earth_rad /
                                   (Constants::mean_earth_rad + orbit_height);
          float sinsatview = sin(satzenith_rad) * orbit_height_frac;
          float snad2 = sinsatview * sinsatview;
          float cnad2 = 1.0f - snad2;

          for (size_t ichan = 0; ichan < nchans; ++ichan) {
            // get index for zhpol and zvpol arrays to determine fraction of horizontal and vertical
            // polatisation for channel
            size_t ipol = polindex[ichan];

            // calculate permittivity using single Debye formula following Hewison & English (1999)
            float fen = freqghz[ichan] / freqr;
            float fen2 = fen * fen;
            float den1r = 1.0f / (1.0f + fen2);

            float perm_real = (perm_static + perm_infinite * fen2) * den1r;
            float perm_imag = fen * (perm_static - perm_infinite) * den1r;
            std::complex<float> xperm(perm_real, perm_imag);

            // calculate complex reflection coefficients (and corrections)
            std::complex<float> perm1 = sqrt(xperm - szen2);
            std::complex<float> perm2 = xperm * czen;

            std::complex<float> rhth = (czen - perm1) / (czen + perm1);
            std::complex<float> rvth = (perm2 - perm1) / (perm2 + perm1);

            float rverts = std::norm(rvth);
            float rhorzs = std::norm(rhth);

            // although the Hewison and English model has scope to account for surface (Bragg)
            // scattering from small-scale features ('roughness') using an empirical term, the
            // parameter is set to 0 for grease ice and is not calculated here.

            // Small scale roughness correction
            // float delta  = 4.0f * M_PI * freqr * 0.1f * small_rough
            // float delta2 = delta * delta
            // float small_rough_cor = exp(-delta2 * czen2)
            float small_rough_cor = 1.0f;

            // specular reflectivities calculated from Fresnel’s formula considerably overestimate
            // the difference between polarisations so a 'large scale roughness correction'
            // is implemented to determine the depolarisation due to this using the expression
            // float qdepol = 0.35f - 0.35f * exp(-0.60f * freq_ghz * large_rough * large_rough)

            // This parameter can take any value from 0.0 for a perfectly specular reflector
            //                                     to 0.5 for a perfectly Lambertian surface.
            // however in the version of the code used operationally this is set to 0.15
            float const qdepol = 0.15f;

            // calculate emissivity from corrected specular reflectivity in h and v
            float evertr = 1.0f - rverts * small_rough_cor;
            float ehorzr = 1.0f - rhorzs * small_rough_cor;

            // calculate emissivity in v and h polarisation accounting for depolarisation
            float evert = evertr * (1.0f - qdepol) + ehorzr * qdepol;
            float ehorz = ehorzr * (1.0f - qdepol) + evertr * qdepol;

            // calculate emissivity in channel polarisation
            float pem = evert * (zvpol[0][ipol] + zvpol[1][ipol] * snad2 + zvpol[2][ipol] * cnad2) +
                        ehorz * (zhpol[0][ipol] + zhpol[1][ipol] * snad2 + zhpol[2][ipol] * cnad2);

            // Finally assign emissivity to obsfunction output assuming that remaining surface
            // fraction is covered with 'grease ice'
            out[ichan][iloc] =
              ice_area_frac[iloc] * ice_emiss + (1.0f - ice_area_frac[iloc]) * pem;
          }
        }
      }
    }
  }
  // -----------------------------------------------------------------------------

  const ufo::Variables & SetSeaIceEmiss::requiredVariables() const {
    return invars_;
  }

  // -----------------------------------------------------------------------------

}  // namespace ufo

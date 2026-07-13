/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/SymmCldImpactUtils.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "oops/util/missingValues.h"

namespace ufo {

// -----------------------------------------------------------------------------

void computeSCIHarnish(const std::vector<float> & bak,
                       const std::vector<float> & obs,
                       float btlim,
                       size_t nlocs,
                       size_t jvar,
                       size_t nvars,
                       int order,
                       std::vector<float> & out) {
  const float fmiss = util::missingValue<float>();
  float Cmod, Cobs, sci;

  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    // Valid IR BT lies within (0, 400) K; this rejects missing values (large
    // negative, in both their float and double encodings), NaN/Inf, and fill
    // values (large positive).
    if (bak[jloc] > 0.0f && bak[jloc] < 400.0f &&
        obs[jloc] > 0.0f && obs[jloc] < 400.0f) {
      Cmod = std::max(0.0f, btlim - bak[jloc]);
      Cobs = std::max(0.0f, btlim - obs[jloc]);
      sci = 0.5f * (Cmod + Cobs);
      out[jloc * nvars + jvar] = std::pow(sci, order);
    } else {
      out[jloc * nvars + jvar] = fmiss;
    }
  }
}

// -----------------------------------------------------------------------------

void computeSCIOkamoto(const std::vector<float> & clr,
                       const std::vector<float> & bak,
                       const std::vector<float> & obs,
                       bool scale_by_omb,
                       float sigmoid_c1,
                       float sigmoid_c2,
                       size_t nlocs,
                       size_t jvar,
                       size_t nvars,
                       int order,
                       std::vector<float> & out) {
  const float fmiss = util::missingValue<float>();
  float Cmod, Cobs, Comb, dx, frac, sci;

  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    // Valid IR BT lies within (0, 400) K; this rejects missing values (large
    // negative, in both their float and double encodings), NaN/Inf, and fill
    // values (large positive).
    if (bak[jloc] > 0.0f && bak[jloc] < 400.0f &&
        obs[jloc] > 0.0f && obs[jloc] < 400.0f) {
      // Fall back to all-sky BT when clear-sky BT is anomalous, e.g. missing
      // or near-zero (CRTM sometimes outputs zero for clear-sky BT when
      // Cloud_Seeding = False).
      float clrval = (clr[jloc] > 0.0f && clr[jloc] < 400.0f) ? clr[jloc] : bak[jloc];
      Cmod = std::abs(clrval - bak[jloc]);
      Cobs = std::abs(clrval - obs[jloc]);
      sci = 0.5f * (Cmod + Cobs);
      if (scale_by_omb) {
        Comb = std::min(std::abs(obs[jloc] - bak[jloc]), 100.0f);
        dx = (Comb - sigmoid_c2) * 0.01f;
        frac = std::max(0.1f, 1.0f / (1.0f + std::exp(-sigmoid_c1 * dx)));
        sci = frac * sci;
      }
      out[jloc * nvars + jvar] = std::pow(sci, order);
    } else {
      out[jloc * nvars + jvar] = fmiss;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

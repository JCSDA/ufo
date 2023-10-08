/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/WindsSPDBCheck.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<WindsSPDBCheck> makerObsFuncWindsSPDBCheck_("WindsSPDBCheck");

// -----------------------------------------------------------------------------

WindsSPDBCheck::WindsSPDBCheck(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // Initialize error_min, and max from options. Make sure they are sane.
  const std::vector<float> &error_min = options_.error_min.value();
  const std::vector<float> &error_max = options_.error_max.value();
  const std::vector<float> &cgross = options_.cgross.value();
  const std::vector<int> &wndtype = options_.wndtype.value();
  const std::string checkvars = options_.checkvars.value();

  ASSERT(error_min < error_max);

  // We need to retrieve the observed wind components.
  invars_ += Variable("ObsValue/windEastward");
  invars_ += Variable("ObsValue/windNorthward");

  // We need to retrieve the observed wind components.
  invars_ += Variable("ObsType/" + checkvars);
  invars_ += Variable("MetaData/pressure");

  // Typical use would be HofX group, but during testing, we include option for GsiHofX
  std::string test_hofx = options_.test_hofx.value();
  invars_ += Variable(test_hofx + "/windEastward");
  invars_ += Variable(test_hofx + "/windNorthward");

  // The starting (un-inflated) value of obserror. If running in sequence of filters,
  // then it is probably found in ObsErrorData, otherwise, it is probably ObsError.
  std::string errgrp = options_.testObserr.value();
  std::string flggrp = options_.testQCflag.value();
  invars_ += Variable(flggrp + "/" + checkvars);
  invars_ += Variable(errgrp + "/" + checkvars);

  // TODO(gthompsn): Need to include a check that whatever HofX/obserr group name used exists.
}

// -----------------------------------------------------------------------------

WindsSPDBCheck::~WindsSPDBCheck() {}

// -----------------------------------------------------------------------------

void WindsSPDBCheck::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();

  // Ensure that only one output variable is expected.
  ASSERT(out.nvars() == 1);

  // Get min, max, gross error values
  const std::vector<float> &error_min = options_.error_min.value();
  const std::vector<float> &error_max = options_.error_max.value();
  const std::vector<float> &cgross = options_.cgross.value();
  const std::vector<int> &wndtype = options_.wndtype.value();
  const std::string checkvars = options_.checkvars.value();

  int nwndtype = wndtype.size();
  ASSERT(nwndtype == error_max.size());
  ASSERT(nwndtype == error_min.size());
  ASSERT(nwndtype == cgross.size());

  std::vector<int> itype(nlocs);
  in.get(Variable("ObsType/" + checkvars), itype);

  std::vector<float> obs_pressure(nlocs);
  in.get(Variable("MetaData/pressure"), obs_pressure);

  // Retrieve Winds observations of wind components
  std::vector<float> u(nlocs), v(nlocs);
  in.get(Variable("ObsValue/windEastward"), u);
  in.get(Variable("ObsValue/windNorthward"), v);

  // Retrieve Model HofX wind components
  const std::string test_hofx = options_.test_hofx.value();
  std::vector<float> um(nlocs), vm(nlocs);
  in.get(Variable(test_hofx + "/windEastward"), um);
  in.get(Variable(test_hofx + "/windNorthward"), vm);

  // Get original ObsError of eastward_wind (would make little sense if diff from northward)
  std::vector<float> currentObserr(nlocs);
  std::vector<int> qcflagdata(nlocs);    //!< effective qcflag
  std::string errgrp = options_.testObserr.value();
  std::string flggrp = options_.testQCflag.value();
  in.get(Variable(errgrp + "/" + checkvars), currentObserr);
  in.get(Variable(flggrp + "/" + checkvars), qcflagdata);

  float uf, vf, presw, qcgross, errmax, errmin;
  float spdb = 0.0f;
  float obserr = 1.0f;
  int ii0 = -1;
  float udiff, vdiff, residual;

  for (size_t jj = 0; jj < nlocs; ++jj) {
    obserr = currentObserr[jj];
    // huge bounds given initially, i.e., no gross-check
    out[0][jj] = 1.e20;

    ii0 = -1;
    for (size_t ii = 0; ii < nwndtype; ii++) {
      if (itype[jj] == wndtype[ii]) {
        ii0 = ii;
        break;
      }
    }
    if (ii0 != -1) {
      qcgross = cgross[ii0];
      errmax = error_max[ii0];
      errmin = error_min[ii0];
      udiff = u[jj]-um[jj];
      vdiff = v[jj]-vm[jj];
      residual = sqrt(udiff*udiff + vdiff*vdiff);
      if (itype[jj] == 244 || itype[jj] == 245 || itype[jj] == 246 ||
          itype[jj] == 253 || itype[jj] == 254 || itype[jj] == 257 ||
          itype[jj] == 258 || itype[jj] == 259) {
         //  Tightening the gross-check bounds for AVHRR/POES IR(244),
         //  NESDIS IR LONG-WAVE(245), NESDIS IMAGER WATER VAPOR(246),
         //  EUMETSAT IR LONG-WAVE AND VISIBLE(253), EUMETSAT IMAGER WATER VAPOR(254)
         //  MODIS/POES IR LONG-WAVE(257),MODIS/POES IMAGER WATER VAPOR(258)
         //  and MODIS/POES IMAGER WATER VAPOR(259) when SPDB < 0.0
        if (u[jj] != missing && v[jj] != missing) {
          uf = um[jj];
          vf = vm[jj];
          presw = 0.01 * obs_pressure[jj];
          spdb = sqrt(u[jj]*u[jj]+v[jj]*v[jj]) - sqrt(uf*uf+vf*vf);
          if (spdb < 0.0f) {
            if (itype[jj] == 244) {
              qcgross = 0.7*qcgross;
            } else if (itype[jj] == 245 || itype[jj] == 246) {
              if (presw > 300.0 && presw < 400.0) qcgross = 0.7 * qcgross;
            } else if (itype[jj] == 253 || itype[jj] == 254) {
              if (presw > 200.0 && presw < 400.0) qcgross = 0.7 * qcgross;
            } else if (itype[jj] >= 257 && itype[jj] <= 259) {
              qcgross = 0.7*qcgross;
            }
          }
        }
      }
      if (qcflagdata[jj] == 0) {
        obserr = std::max(errmin, std::min(obserr, errmax));
        if (checkvars == "windEastward")
               out[0][jj] = qcgross*obserr*(std::abs(udiff)/residual);
        if (checkvars == "windNorthward")
               out[0][jj] = qcgross*obserr*(std::abs(vdiff)/residual);
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & WindsSPDBCheck::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

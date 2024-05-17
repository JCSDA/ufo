/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorBoundMW.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsErrorFactorLatRad.h"
#include "ufo/filters/obsfunctions/ObsErrorFactorTransmitTopRad.h"
#include "ufo/filters/obsfunctions/ObsErrorModelRamp.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorBoundMW> makerObsErrorBoundMW_("ObsErrorBoundMW");

// -----------------------------------------------------------------------------

ObsErrorBoundMW::ObsErrorBoundMW(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get sensor information from options
  const std::string &sensor = options_.sensor.value();

  // Get instrument and satellite from sensor
  std::string inst, sat;
  splitInstSat(sensor, inst, sat);
  ASSERT(inst == "amsua" || inst == "atms" || inst == "mhs");

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Get test groups from options
  const std::string &errgrp = options_.testObserr.value();
  const std::string &flaggrp = options_.testQCflag.value();

  // Include list of required data from ObsSpace
  invars_ += Variable(flaggrp+"/brightnessTemperature", channels_);
  invars_ += Variable(errgrp+"/brightnessTemperature", channels_);

  // Include list of required data from GeoVaLs
  invars_ += Variable("GeoVaLs/water_area_fraction");

  const Variable &obserrlat = options_.obserrBoundLat.value();
  invars_ += obserrlat;

  const Variable &obserrtaotop = options_.obserrBoundTransmittop.value();
  invars_ += obserrtaotop;

  const Variable &obserrtopo = options_.obserrBoundTopo.value();
  invars_ += obserrtopo;

  // const Variable &obserr = options_.obserrFunction.value();
  // invars_ += obserr;
  if (options_.obserrFunction.value() != boost::none) {
    const boost::optional<Variable> &obserrvar = options_.obserrFunction.value();
    invars_ += *obserrvar;
  }

  if (options_.obserrOriginal.value() != boost::none) {
    const std::vector<float> &obserr0 = options_.obserrOriginal.value().get();
  }
}

// -----------------------------------------------------------------------------

ObsErrorBoundMW::~ObsErrorBoundMW() {}

// -----------------------------------------------------------------------------

void ObsErrorBoundMW::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get sensor information from options
  const std::string &sensor = options_.sensor.value();

  // Get instrument and satellite from sensor
  std::string inst, sat;
  splitInstSat(sensor, inst, sat);
  ASSERT(inst == "amsua" || inst == "atms" || inst == "mhs");

  // Get dimensions
  size_t nlocs = in.nlocs();
  size_t nchans = channels_.size();

  // Get observation error bounds from options
  const std::vector<float> &obserr_bound_max = options_.obserrBoundMax.value();

  // Get area fraction of each surface type
  std::vector<float> water_frac(nlocs);
  in.get(Variable("GeoVaLs/water_area_fraction"), water_frac);

  // Get error factor from ObsFunction
  const Variable &obserrlat = options_.obserrBoundLat.value();
  ioda::ObsDataVector<float> errflat(in.obsspace(), obserrlat.toOopsObsVariables());
  in.get(obserrlat, errflat);

  // Get error factor from ObsFunction
  const Variable &obserrtaotop = options_.obserrBoundTransmittop.value();
  ioda::ObsDataVector<float> errftaotop(in.obsspace(), obserrtaotop.toOopsObsVariables());
  in.get(obserrtaotop, errftaotop);

  // Get error factor from ObsFunction
  const Variable &obserrtopo = options_.obserrBoundTopo.value();
  ioda::ObsDataVector<float> errftopo(in.obsspace(), obserrtopo.toOopsObsVariables());
  in.get(obserrtopo, errftopo);

  // Get all-sky observation error from ObsFunction
  std::unique_ptr<ioda::ObsDataVector<float>> obserr;
  if (options_.obserrFunction.value() != boost::none) {
    const boost::optional<Variable> &obserrvar = options_.obserrFunction.value();
    obserr.reset(new ioda::ObsDataVector<float>(in.obsspace(),
                 (*obserrvar).toOopsObsVariables()));
    in.get(*obserrvar, *obserr);
  }

  // Set channel numbers
  int ich238, ich314, ich503, ich528, ich536, ich544, ich549, ich890;
  if (inst == "amsua") {
    ich238 = 1, ich314 = 2, ich503 = 3, ich528 = 4, ich536 = 5;
    ich544 = 6, ich549 = 7, ich890 = 15;
  } else if (inst == "atms") {
    ich238 = 1, ich314 = 2, ich503 = 3, ich528 = 5, ich536 = 6;
    ich544 = 7, ich549 = 8, ich890 = 16;
  }
  float thresholdfactor = 3.0;
  if (options_.thresholdfactor.value() != boost::none) {
    thresholdfactor = options_.thresholdfactor.value().get();
  }

  // Output integrated error bound for gross check
  std::vector<float> obserrdata(nlocs);
  std::vector<int> qcflagdata(nlocs);
  const std::string &errgrp = options_.testObserr.value();
  const std::string &flaggrp = options_.testQCflag.value();
  const float missing = util::missingValue<float>();
  float varinv = 0.0;
  if (options_.obserrFunction.value() != boost::none) {
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      int channel = ichan + 1;
      in.get(Variable(flaggrp+"/brightnessTemperature", channels_)[ichan], qcflagdata);
      in.get(Variable(errgrp+"/brightnessTemperature", channels_)[ichan], obserrdata);
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        if (flaggrp == "PreQC") obserrdata[iloc] == missing ? qcflagdata[iloc] = 100
                                                            : qcflagdata[iloc] = 0;
        (qcflagdata[iloc] == 0) ? (varinv = 1.0 / pow(obserrdata[iloc], 2)) : (varinv = 0.0);
        out[ichan][iloc] = (*obserr)[ichan][iloc];
        if (varinv > 0.0) {
          if (water_frac[iloc] > 0.99) {
            if (inst == "amsua") {
              if (channel <= ich536  || channel == ich890) {
                out[ichan][iloc] = thresholdfactor * (*obserr)[ichan][iloc]
                                       * (1.0 / pow(errflat[0][iloc], 2))
                                       * (1.0 / pow(errftaotop[ichan][iloc], 2))
                                       * (1.0 / pow(errftopo[ichan][iloc], 2));
              } else {
                out[ichan][iloc] = std::fmin((thresholdfactor * (*obserr)[ichan][iloc]
                                       * (1.0 / pow(errflat[0][iloc], 2))
                                       * (1.0 / pow(errftaotop[ichan][iloc], 2))
                                       * (1.0 / pow(errftopo[ichan][iloc], 2))),
                                          obserr_bound_max[ichan]);
              }
            }
            if (inst == "atms") {
              if (channel <= ich536  || channel >= ich890) {
                out[ichan][iloc] = std::fmin((thresholdfactor * (*obserr)[ichan][iloc]
                                       * (1.0 / pow(errflat[0][iloc], 2))
                                       * (1.0 / pow(errftaotop[ichan][iloc], 2))
                                       * (1.0 / pow(errftopo[ichan][iloc], 2))), 10.0);
              } else {
                out[ichan][iloc] = std::fmin((thresholdfactor * (*obserr)[ichan][iloc]
                                       * (1.0 / pow(errflat[0][iloc], 2))
                                       * (1.0 / pow(errftaotop[ichan][iloc], 2))
                                       * (1.0 / pow(errftopo[ichan][iloc], 2))),
                                          obserr_bound_max[ichan]);
              }
            }
            if (inst == "mhs") {
//            GEOS all-sky only.
              out[ichan][iloc] = std::fmin((thresholdfactor * (*obserr)[ichan][iloc]
                                     * (1.0 / pow(errflat[0][iloc], 2))
                                     * (1.0 / pow(errftaotop[ichan][iloc], 2))
                                     * (1.0 / pow(errftopo[ichan][iloc], 2))),
                                        obserr_bound_max[ichan]);
            }
          } else {
            out[ichan][iloc] = std::fmin((thresholdfactor * (*obserr)[ichan][iloc]
                                   * (1.0 / pow(errflat[0][iloc], 2))
                                   * (1.0 / pow(errftaotop[ichan][iloc], 2))
                                   * (1.0 / pow(errftopo[ichan][iloc], 2))),
                                      obserr_bound_max[ichan]);
          }
        }
      }
    }
  }

  std::vector<float> obserr0(nchans);
  if (options_.obserrOriginal.value() != boost::none) {
    const std::vector<float> obserr0 = options_.obserrOriginal.value().get();

    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      int channel = ichan + 1;
      in.get(Variable(flaggrp+"/brightnessTemperature", channels_)[ichan], qcflagdata);
      in.get(Variable(errgrp+"/brightnessTemperature", channels_)[ichan], obserrdata);
      for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        if (flaggrp == "PreQC") obserrdata[iloc] == missing ? qcflagdata[iloc] = 100
                                                            : qcflagdata[iloc] = 0;
        (qcflagdata[iloc] == 0) ? (varinv = 1.0 / pow(obserrdata[iloc], 2)) : (varinv = 0.0);
        out[ichan][iloc] = obserr0[ichan];
        if (varinv > 0.0) {
            out[ichan][iloc] = std::fmin((thresholdfactor * obserr0[ichan]
                                   * (1.0 / pow(errflat[0][iloc], 2))
                                   * (1.0 / pow(errftaotop[ichan][iloc], 2))
                                   * (1.0 / pow(errftopo[ichan][iloc], 2))),
                                      obserr_bound_max[ichan]);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorBoundMW::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

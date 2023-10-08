/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#include "ufo/filters/obsfunctions/MetOfficeAllSkyErrorModel.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/SurfaceReportConstants.h"

namespace ufo {

static ObsFunctionMaker<MetOfficeAllSkyErrorModel>
  makerMetOfficeAllSkyErrorModel_("MetOfficeAllSkyErrorModel");

// -----------------------------------------------------------------------------

MetOfficeAllSkyErrorModel::MetOfficeAllSkyErrorModel(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.deserialize(conf);

  // Get channels from options
  std::set<int> channelset = oops::parseIntSet(options_.channelList);
  std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  ASSERT(channels_.size() > 0);

  // Include list of required data from ObsSpace
  invars_ += Variable("OneDVar/liquidWaterPath");
  invars_ += Variable("OneDVar/iceWaterPath");
  invars_ += Variable("OneDVar/transmittance", channels_);
  invars_ += Variable("MetaData/surfaceQualifier");
  invars_ += Variable("MetaDataError/instrumentNoise", channels_);
  invars_ += Variable("SurfEmiss/emissivityError", channels_);
  invars_ += Variable("QCflagsData/brightnessTemperature", channels_);
}

// -----------------------------------------------------------------------------

MetOfficeAllSkyErrorModel::~MetOfficeAllSkyErrorModel() {}

// -----------------------------------------------------------------------------

void MetOfficeAllSkyErrorModel::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  // Get parameters from options
  const std::vector<float> &fixland = options_.fixland.value().get();
  const std::vector<float> &fixsea = options_.fixsea.value().get();
  const std::vector<float> &fixice = options_.fixice.value().get();
  const std::vector<float> &taulinsea = options_.taulinsea.value().get();
  const std::vector<float> &taulinland = options_.taulinland.value().get();
  const std::vector<float> &taulinice = options_.taulinice.value().get();
  const std::vector<float> &tausqsea = options_.tausqsea.value().get();
  const std::vector<float> &tausqland = options_.tausqland.value().get();
  const std::vector<float> &tausqice = options_.tausqice.value().get();
  const std::vector<float> &lwpcoef = options_.lwpcoef.value().get();
  const std::vector<float> &iwpcoef = options_.iwpcoef.value().get();
  const float ScaleVarErrorAboveEmissError = options_.ScaleVarErrorAboveEmissError.value();
  const bool UseEmissivityAtlas = options_.UseEmissivityAtlas.value();
  const float maxlwp = options_.maxlwp.value();
  const float maxiwp = options_.maxiwp.value();

  // Get dimensions
  const size_t nlocs = out.nlocs();
  const size_t nchans = out.nvars();

  std::vector<float> ob_lwp(nlocs);
  std::vector<float> ob_iwp(nlocs);
  std::vector<int> ob_surf_type(nlocs);

  const float missing = util::missingValue<float>();
  const int missingValueInt = util::missingValue<int>();

  if (in.has(Variable("OneDVar/liquidWaterPath"))) {
    in.get(Variable("OneDVar/liquidWaterPath"), ob_lwp);
  } else {
    std::fill(ob_lwp.begin(), ob_lwp.end(), missing);
  }

  if (in.has(Variable("OneDVar/iceWaterPath"))) {
    in.get(Variable("OneDVar/iceWaterPath"), ob_iwp);
  } else {
    std::fill(ob_iwp.begin(), ob_iwp.end(), missing);
  }

  if (in.has(Variable("MetaData/surfaceQualifier"))) {
    in.get(Variable("MetaData/surfaceQualifier"), ob_surf_type);
  } else {
    std::fill(ob_surf_type.begin(), ob_surf_type.end(), missingValueInt);
    // surface type must be defined
    throw eckit::UserError("MetaData/surfaceQualifier not found");
  }

  std::vector<float> ob_instr_noise(nlocs);
  std::vector<float> surf_emiss_error(nlocs);  // should be error on surface emissivity
  std::vector<float> tausurf(nlocs);
  std::vector<int> flagsQC(nlocs);
  float VarErrorScaling;

  RTTOV_surftype surftype;

  for (size_t ich = 0; ich < nchans; ++ich) {
    // mandatory NEDT for all channels
    if (in.has(Variable("MetaDataError/instrumentNoise", channels_)[ich])) {
      in.get(Variable("MetaDataError/instrumentNoise", channels_)[ich], ob_instr_noise);
    } else {
      throw eckit::UserError("MetaDataError/instrumentNoise not found");
    }

    // transmittance
    if (in.has(Variable("OneDVar/transmittance", channels_)[ich]) &&
      (options_.taulinsea.value() != boost::none || options_.taulinice.value() != boost::none ||
       options_.taulinland.value() != boost::none || options_.tausqsea.value() != boost::none ||
       options_.tausqice.value() != boost::none || options_.tausqland.value() != boost::none) ) {
      in.get(Variable("OneDVar/transmittance", channels_)[ich], tausurf);
      oops::Log::debug() << "Transmittance: " << tausurf << std::endl;
    } else {
      std::fill(tausurf.begin(), tausurf.end(), missing);
      oops::Log::warning() << "OneDVar/transmittance not found or not needed" << std::endl;
    }

    // surface MW emissivity error
    if (in.has(Variable("SurfEmiss/emissivityError", channels_)[ich]) &&
        (options_.taulinland.value() != boost::none
         || options_.tausqland.value() != boost::none) ) {
      in.get(Variable("SurfEmiss/emissivityError", channels_)[ich], surf_emiss_error);
    } else {
      std::fill(surf_emiss_error.begin(), surf_emiss_error.end(), 0.0f);
      oops::Log::warning() << "SurfEmiss/emissivityError not found or not needed"
                           << std::endl;
    }

    if (in.has(Variable("QCflagsData/brightnessTemperature", channels_)[ich])) {
      in.get(Variable("QCflagsData/brightnessTemperature", channels_)[ich], flagsQC);
    } else {
      std::fill(flagsQC.begin(), flagsQC.end(), missingValueInt);
    }

    for (size_t iloc=0; iloc < nlocs; ++iloc) {
        VarErrorScaling = missing;

        // skip adding to obs error when the the ob is rejected at location loc and channel ich
        if (flagsQC[iloc] != 0) {
          out[ich][iloc] = missing;
          continue;
        } else {
          out[ich][iloc] = 0.0f;
        }

        // mandatory NEDT for all channels
        if (ob_instr_noise[iloc] != missing) {
          out[ich][iloc] += ob_instr_noise[iloc]*ob_instr_noise[iloc];
        } else {
          throw eckit::BadValue("MetaDataError/instrumentNoise has invalid value");
        }

        // over sea
        if (options_.fixsea.value() != boost::none && ob_surf_type[iloc] == surftype.sea) {
          out[ich][iloc] += fixsea[ich]*fixsea[ich];
        }

        // over sea ice
        if (options_.fixice.value() != boost::none && ob_surf_type[iloc] == surftype.seaice) {
          out[ich][iloc] += fixice[ich]*fixice[ich];
        }  else if (options_.fixland.value() != boost::none
                    && ob_surf_type[iloc] != surftype.sea) {
          // over land, also over sea ice if fixice option is not provided
          out[ich][iloc] += fixland[ich]*fixland[ich];
        }

        // default value for trasmittance outside valid range is 1.0
        float predictor = (tausurf[iloc] < 0.0f) ? 1.0f : tausurf[iloc];

        // TauLinSea applied over sea only. Multiplied by surface to space trans
        if (options_.taulinsea.value() != boost::none && ob_surf_type[iloc] == surftype.sea) {
          out[ich][iloc] += predictor*taulinsea[ich]*predictor*taulinsea[ich];
        }

        // TauLinIce applied over sea ice only. Multiplied by surface to space trans
        if (options_.taulinice.value() != boost::none && ob_surf_type[iloc] == surftype.seaice) {
          out[ich][iloc] += predictor*taulinice[ich]*predictor*taulinice[ich];
        } else if (options_.taulinland.value() != boost::none
                   && ob_surf_type[iloc] != surftype.sea) {
          // TauLinLand applied over land. Note that land coeffs are
          // applied over sea ice if the ice coefficients have not been supplied
          if (ScaleVarErrorAboveEmissError < 1.0f && UseEmissivityAtlas
              && ob_surf_type[iloc] == surftype.land) {
              VarErrorScaling = std::min(std::max(surf_emiss_error[iloc],
                ScaleVarErrorAboveEmissError) / ScaleVarErrorAboveEmissError, 100.0f);
          } else {
              VarErrorScaling = 1.0f;
          }
          out[ich][iloc] += predictor*taulinland[ich]*VarErrorScaling
                            *predictor*taulinland[ich]*VarErrorScaling;
        }

        // TauSqSea applied over sea only. Multiplied by surface to (space trans)**2
        if (options_.tausqsea.value() != boost::none && ob_surf_type[iloc] == surftype.sea) {
          out[ich][iloc] += predictor*predictor*tausqsea[ich]*predictor*predictor*tausqsea[ich];
        }

        // TauSqIce applied over sea ice only. Multiplied by surface to (space trans)**2
        if (options_.tausqice.value() != boost::none && ob_surf_type[iloc] == surftype.seaice) {
          out[ich][iloc] += predictor*predictor*tausqice[ich]*predictor*predictor*tausqice[ich];
        } else if (options_.tausqland.value() != boost::none
                   && ob_surf_type[iloc] != surftype.sea) {
          // TauSqLand Land, or sea ice if ice coeffs not set
          if (VarErrorScaling == missing && ScaleVarErrorAboveEmissError < 1.0f
              && UseEmissivityAtlas && ob_surf_type[iloc] == surftype.land) {
              VarErrorScaling = std::min(std::max(surf_emiss_error[iloc],
                ScaleVarErrorAboveEmissError) / ScaleVarErrorAboveEmissError, 100.0f);
          } else {
              VarErrorScaling = 1.0f;
          }
          out[ich][iloc] += predictor*predictor*tausqland[ich]*VarErrorScaling
            *predictor*predictor*tausqland[ich]*VarErrorScaling;
        }

        // 1D-Var liquidWaterPath and iceWaterPath. Note both of these two predictors
        // must be available or prescribed max values for each of them are used instead.
        if (options_.lwpcoef.value() != boost::none && options_.iwpcoef.value() != boost::none
            && lwpcoef[ich] > 0.0f && iwpcoef[ich] > 0.0f) {
          if (ob_lwp[iloc] != missing &&  ob_iwp[iloc] != missing) {
            out[ich][iloc] = sqrt(out[ich][iloc]) + lwpcoef[ich]*ob_lwp[iloc] +
                    iwpcoef[ich]*ob_iwp[iloc];
            out[ich][iloc] *= out[ich][iloc];
          } else {
            out[ich][iloc] = sqrt(out[ich][iloc]) + lwpcoef[ich]*maxlwp +
                      iwpcoef[ich]*maxiwp;
            out[ich][iloc] *= out[ich][iloc];
          }
        }

        // final step: convert error to standard deviation
        out[ich][iloc] = sqrt(out[ich][iloc]);
     }
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & MetOfficeAllSkyErrorModel::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

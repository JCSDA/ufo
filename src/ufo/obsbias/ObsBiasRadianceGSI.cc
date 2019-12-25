/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */
#include "ufo/obsbias/ObsBiasRadianceGSI.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <utility>

#include "ObsBiasRadianceGSI.interface.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/Constants.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

namespace ufo {

static ObsBiasMaker<ObsBiasRadianceGSI> makerBiasRadianceGSI_("GSI");

// -----------------------------------------------------------------------------

ObsBiasRadianceGSI::ObsBiasRadianceGSI(const ioda::ObsSpace & odb,
                                       const eckit::Configuration & conf)
  : ObsBiasBase(), odb_(odb), geovars_(), hdiags_(), tlapmean_(),
    newpc4pred_(false), adp_anglebc_(false), emiss_bc_(false), predictors_(), predNames_() {
// Default predictor names from GSI
  predictors_ = {"constant",
                 "scan_angle",
                 "cloud_liquid_water",
                 "lapse_rate_squared",
                 "lapse_rate",
                 "cosine_of_latitude_times_orbit_node",
                 "sine_of_latitude",
                 "emissivity",
                 "scan_angle_4th_order",
                 "scan_angle_3rd_order",
                 "scan_angle_2nd_order",
                 "scan_angle_1st_order"
                };
// Retrive the channels
  channels_ = odb_.obsvariables().channels();

// Parse predictors from the conf
  if (conf.has("ObsBias.predictors")) {
    predictors_.clear();
    predictors_ = conf.getStringVector("ObsBias.predictors");
    predNames_.reset(new oops::Variables(predictors_, channels_));
  }

// GeoVals needed from model
  const std::vector<std::string> vv0{"air_temperature",
                                     "air_pressure",
                                     "water_area_fraction",
                                     "average_surface_temperature_within_field_of_view"
                                    };
  geovars_.reset(new oops::Variables(vv0));

// Hdiags needed from H diagnostics
  std::vector<std::string> vv{"brightness_temperature_jacobian_surface_emissivity_CH",
                              "transmittances_of_atmosphere_layer_CH",
                              "brightness_temperature_CH"
                              };
  hdiags_.reset(new oops::Variables(vv));

// Parse Sensor_ID from the conf
  const eckit::LocalConfiguration obsoprconf(conf, "ObsOperator");
  sensor_id_ = obsoprconf.getString("ObsOptions.Sensor_ID");

// Replace "_CH" in hdiags_ with digitial Channel ID
  std::vector<std::string> vvtmp;
  for (std::size_t jv = 0; jv < vv.size(); ++jv) {
    std::size_t found = vv[jv].find("_CH");
    if (found != std::string::npos) {
      for (std::size_t jc = 0; jc < channels_.size(); ++jc)
        vvtmp.push_back(vv[jv].replace(found, 3, '_'+std::to_string(channels_[jc])));
    } else {
      vvtmp.push_back(vv[jv]);
    }
  }
  if (vvtmp.size() > 0) hdiags_.reset(new oops::Variables(vvtmp));

// Read ObsBias parameters first guess if available
  const eckit::LocalConfiguration biasConf(conf, "ObsBias");
  if (biasConf.has("abias_in")) {
    read(biasConf);
  } else {
    biascoeffs_.resize(channels_.size() * predictors_.size(), 0.0);
  }

  newpc4pred_ = biasConf.getBool("newpc4pred", false);

  adp_anglebc_ = biasConf.getBool("adp_anglebc", false);

  emiss_bc_ = biasConf.getBool("emiss_bc", false);

  oops::Log::trace() << "ObsBiasRadianceGSI created." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSI::read(const eckit::Configuration & biasconf) {
  oops::Log::trace() << "ObsBiasRadianceGSI::read from file: "
                     << biasconf.getString("abias_in") << " Starting "<< std::endl;
  const std::string filename = biasconf.getString("abias_in");
  std::ifstream infile(filename);

  std::size_t ich;     //  sequential number
  std::string nusis;   //  sensor/instrument/satellite
  std::size_t nuchan;  //  channel number
  float tlap, tsum;
  std::size_t ntlapupdate;

  if ( infile.is_open() )
  {
    float par;
    while ( !infile.eof() )
    {
      infile >> ich;
      infile >> nusis;
      infile >> nuchan;
      infile >> tlap;
      infile >> tsum;
      infile >> ntlapupdate;
      if ( nusis == sensor_id_ &&
           std::find(channels_.begin(), channels_.end(), nuchan) != channels_.end() ) {
        tlapmean_.push_back(tlap);
        for (std::size_t j=0; j < predictors_.size(); ++j) {
          infile >> par;
          biascoeffs_.push_back(static_cast<double>(par));
        }
      } else {
        for (std::size_t j=0; j < predictors_.size(); ++j) {
          infile >> par;
        }
      }
    }
    infile.close();
    oops::Log::trace() << "ObsBiasRadianceGSI::read from file: "
                      << biasconf.getString("abias_in") << " Done " << std::endl;
  } else {
    oops::Log::error() << "Unable to open file : " << filename << std::endl;
    ABORT("Unable to open bias correction parameters file ");
  }
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSI::write(const eckit::Configuration & conf) const {
  oops::Log::trace() << "ObsBiasRadianceGSI::write to file not implmented: "
                     << conf.getString("abias_out") << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSI::computeObsBias(ioda::ObsVector & ybias,
                                        std::unique_ptr<ioda::ObsDataVector<float>> & predTerms)
                                        const {
  const std::size_t npred = predictors_.size();
  const std::size_t nchanl = channels_.size();
  const std::size_t nlocs = ybias.nlocs();
  ASSERT(ybias.nlocs() == odb_.nlocs());

  for (std::size_t n = 0; n < npred; ++n) {
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        (*predTerms)[n*nchanl+jc][jl] = biascoeffs_[jc*npred+n] * (*predTerms)[n*nchanl+jc][jl];
      }
    }
  }

  ybias.zero();
  // Loop through each location
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    // Loop through each channel
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      // Linear combination
      for (std::size_t n = 0; n < npred; ++n) {
        ybias[jl*nchanl+jc] += (*predTerms)[n*nchanl+jc][jl];
      }
    }
  }

  oops::Log::trace() << "ObsBiasRadianceGSI::computeObsBias done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSI::computeObsBiasPredictors(
                                    const GeoVaLs & geovals,
                                    const ObsDiagnostics & ydiags,
                                    std::unique_ptr<ioda::ObsDataVector<float>> & preds) const {
  const std::size_t npred = predictors_.size();
  const std::size_t nlocs = odb_.nlocs();
  const std::size_t nchanl = channels_.size();

  ASSERT(preds && preds->nvars() == npred*nchanl);

  // Following variables should be moved to yaml file ?
  const float ssmis_precond = 0.01;  //  default preconditioner for ssmis bias terms

  // Obstype
  bool no85GHz {false};  // revisit please
  bool amsre  {sensor_id_.find("amsre")  != std::string::npos};
  bool ssmi   {sensor_id_.find("ssmi")   != std::string::npos &&
               sensor_id_.find("ssmis")  == std::string::npos};
  bool ssmis  {sensor_id_.find("ssmis")  != std::string::npos};
  bool amsua  {sensor_id_.find("amsua")  != std::string::npos};
  bool atms   {sensor_id_.find("atms")   != std::string::npos};
  bool amsr2  {sensor_id_.find("amsr2")  != std::string::npos};
  bool gmi    {sensor_id_.find("gmi")    != std::string::npos};
  bool saphir {sensor_id_.find("saphir") != std::string::npos};

  // common vectors storage
  std::vector <float> pred(nlocs, 0.0);

  std::vector<float> profile;

  std::vector<float> zasat(nlocs, 0.0);
  if (std::find(predictors_.begin(), predictors_.end(),
              "scan_angle") != predictors_.end() ||
      std::find(predictors_.begin(), predictors_.end(),
              "cloud_liquid_water") != predictors_.end()) {
    odb_.get_db("MetaData", "sensor_zenith_angle", zasat);
  }

  std::vector<float> cenlat(nlocs, 0.0);
  std::vector<float> node(nlocs, 0.0);
  if (std::find(predictors_.begin(), predictors_.end(),
              "sine_of_latitude") != predictors_.end() ||
      std::find(predictors_.begin(), predictors_.end(),
              "cosine_of_latitude_times_orbit_node") != predictors_.end()) {
    odb_.get_db("MetaData", "latitude", cenlat);
    odb_.get_db("MetaData", "sensor_azimuth_angle", node);
  }

  // retrieve the average surface temperature
  std::vector<float> tsavg5(nlocs, 0.0);
  if (std::find(predictors_.begin(), predictors_.end(),
              "lapse_rate_squared") != predictors_.end() ||
      std::find(predictors_.begin(), predictors_.end(),
              "lapse_rate") != predictors_.end()         ||
      std::find(predictors_.begin(), predictors_.end(),
              "cloud_liquid_water") != predictors_.end()) {
    geovals.get(tsavg5, "average_surface_temperature_within_field_of_view");
  }

  // calculate the lapse rate
  std::vector<std::vector<float>> tlap(nchanl , std::vector<float>(nlocs, 0.0));
  if (std::find(predictors_.begin(), predictors_.end(),
              "lapse_rate_squared") != predictors_.end() ||
      std::find(predictors_.begin(), predictors_.end(),
              "lapse_rate") != predictors_.end()) {
    // Retrieve the transmittances_of_atmosphere_layer from Hdiag
    std::string varname;
    std::vector<std::vector<std::vector<float>>> ptau5;
    std::vector<std::vector<float>> tmpvar;
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      varname = "transmittances_of_atmosphere_layer_" +
                std::to_string(channels_[jc]);
      tmpvar.clear();
      for (std::size_t js = 0; js < ydiags.nlevs(varname); ++js) {
        ydiags.get(pred, varname, js+1);
        tmpvar.push_back(pred);
      }
      ptau5.push_back(tmpvar);
    }
    // Retrieve the temperature
    std::vector<std::vector<float>> tvp;
    std::size_t nlevs = geovals.nlevs("air_temperature");
    for (std::size_t js = 0; js < nlevs; ++js) {
      geovals.get(pred, "air_temperature", js+1);
      tvp.push_back(pred);
    }
    nlevs = geovals.nlevs("air_pressure");
    float tlapchn;
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
          tlapchn = (ptau5[jc][nlevs-2][jl]-ptau5[jc][nlevs-1][jl])*(tsavg5[jl]-tvp[nlevs-2][jl]);
          for (std::size_t k = 1; k < nlevs-1; ++k) {
            tlapchn = tlapchn+(ptau5[jc][nlevs-k-2][jl]-ptau5[jc][nlevs-k-1][jl])*
                      (tvp[nlevs-k][jl]-tvp[nlevs-k-2][jl]);
          }
          if (!newpc4pred_) {
            tlapchn = 0.01*tlapchn;
          }
          tlap[jc][jl] = tlapchn-tlapmean_[jc];
      }
    }
  }

  // retrieve the sensor view angle
  std::vector<float> view_angle(nlocs, 0.0);
  if (std::find(predictors_.begin(), predictors_.end(),
              "scan_angle_4th_order") != predictors_.end() ||
      std::find(predictors_.begin(), predictors_.end(),
              "scan_angle_3rd_order") != predictors_.end() ||
      std::find(predictors_.begin(), predictors_.end(),
              "scan_angle_2nd_order") != predictors_.end() ||
      std::find(predictors_.begin(), predictors_.end(),
              "scan_angle_1st_order") != predictors_.end()) {
    odb_.get_db("MetaData", "sensor_view_angle", view_angle);
  }

  // Compute predictors one-by-one
  std::size_t indx = 0;
  for (std::size_t n = 0; n < npred; ++n) {
    if (predictors_[n] == "constant") {
      if (!newpc4pred_) {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = 0.01;
        }
      } else {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = 1.0;
        }
      }
      indx += nchanl;
    } else if (predictors_[n] == "scan_angle") {
      if (!newpc4pred_) {
        if (ssmi || ssmis || amsre || gmi || amsr2) {
          for (std::size_t jc = 0; jc < nchanl; ++jc) {
            for (std::size_t jl = 0; jl < nlocs; ++jl)
              (*preds)[indx+jc][jl] = 0.0;
          }
        } else {
          for (std::size_t jc = 0; jc < nchanl; ++jc) {
            for (std::size_t jl = 0; jl < nlocs; ++jl)
              (*preds)[indx+jc][jl] = 0.10*pow(1.0/cos(zasat[jl]) - 1.0, 2) - .015;
          }
        }
      } else {
        if (adp_anglebc_) {
          for (std::size_t jc = 0; jc < nchanl; ++jc) {
            for (std::size_t jl = 0; jl < nlocs; ++jl)
              (*preds)[indx+jc][jl] = 0.0;
          }
        } else {
          for (std::size_t jc = 0; jc < nchanl; ++jc) {
            for (std::size_t jl = 0; jl < nlocs; ++jl)
              (*preds)[indx+jc][jl] = pow(1.0/cos(zasat[jl]) - 1.0, 2);
          }
        }
      }
      indx += nchanl;
    } else if (predictors_[n] == "cloud_liquid_water") {
      // Retrieve the brightness temperature observation
      std::vector<std::vector<float>> tb;
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        odb_.get_db("ObsValue", "brightness_temperature_" +
                   std::to_string(channels_[jc]), pred);
        tb.emplace_back(pred);
      }
      // Transpose the array
      std::vector<std::vector<float>> tb_obs;
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        profile.clear();
        for (std::size_t jc = 0; jc < nchanl; ++jc)
          profile.emplace_back(tb[jc][jl]);
        tb_obs.emplace_back(profile);
      }
      // Retrieve the simulated brightness temperature from Hdiag
      tb.clear();
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        ydiags.get(pred, "brightness_temperature_" + std::to_string(channels_[jc]));
        tb.emplace_back(pred);
      }
      // Transpose the array
      std::vector<std::vector<float>> tsim;
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        profile.clear();
        for (std::size_t jc = 0; jc < nchanl; ++jc)
          profile.emplace_back(tb[jc][jl]);
        tsim.emplace_back(profile);
      }
      // Retrieve the surface wind spped from geovals
      std::vector<float> sfc_speed(nlocs, 0.0);
      geovals.get(sfc_speed, "surface_wind_speed");
      // Retrieve the scan_position from ObsSpace
      std::vector<float> nadir(nlocs);
      odb_.get_db("MetaData", "scan_position", nadir);
      std::vector<float> clw(nlocs);
      // Calculate the cloud liquid water
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        calc_clw_f90(static_cast<int>(nadir[jl]), tb_obs[jl][0], tsim[jl][0], channels_[0], nchanl,
                     no85GHz, amsua, ssmi, ssmis, amsre, atms, amsr2, gmi, saphir,
                     tsavg5[jl], sfc_speed[jl], zasat[jl], clw[jl]);
      }
      if (amsre) {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = clw[jl];
        }
      } else {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = clw[jl]*cos(zasat[jl]*cos(zasat[jl]));
        }
      }
      indx += nchanl;
    } else if (predictors_[n] == "lapse_rate_squared") {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl)
          (*preds)[indx+jc][jl] = tlap[jc][jl]*tlap[jc][jl];
      }
      indx += nchanl;
    } else if (predictors_[n] == "lapse_rate") {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl)
          (*preds)[indx+jc][jl] = tlap[jc][jl];
      }
      indx += nchanl;
    } else if (predictors_[n] == "cosine_of_latitude_times_orbit_node") {
      if (ssmis) {
        if (!newpc4pred_) {
          for (std::size_t jc = 0; jc < nchanl; ++jc) {
            for (std::size_t jl = 0; jl < nlocs; ++jl) {
              if (node[jl] < 1000) {
                (*preds)[indx+jc][jl] = ssmis_precond*node[jl]*cos(cenlat[jl]*Constants::deg2rad);
              } else {
                (*preds)[indx+jc][jl] = 0.0;
              }
            }
          }
        } else {
          for (std::size_t jc = 0; jc < nchanl; ++jc) {
            for (std::size_t jl = 0; jl < nlocs; ++jl)
              (*preds)[indx+jc][jl] = node[jl]*cos(cenlat[jl]*Constants::deg2rad);
          }
        }
      } else {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
            for (std::size_t jl = 0; jl < nlocs; ++jl)
              (*preds)[indx+jc][jl] = 0.0;
        }
      }
      indx += nchanl;
    } else if (predictors_[n] == "sine_of_latitude") {
      if (ssmis) {
        if (!newpc4pred_) {
          for (std::size_t jc = 0; jc < nchanl; ++jc) {
            for (std::size_t jl = 0; jl < nlocs; ++jl) {
              if (node[jl] < 1000) {
                (*preds)[indx+jc][jl] = ssmis_precond*sin(cenlat[jl]*Constants::deg2rad);
              } else {
                (*preds)[indx+jc][jl] = 0.0;
              }
            }
          }
        } else {
          for (std::size_t jc = 0; jc < nchanl; ++jc) {
            for (std::size_t jl = 0; jl < nlocs; ++jl)
              (*preds)[indx+jc][jl] = sin(cenlat[jl]*Constants::deg2rad);
          }
        }
      } else {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = 0.0;
        }
      }
      indx += nchanl;
    } else if (predictors_[n] == "emissivity") {
      if (adp_anglebc_ && emiss_bc_) {
        std::vector<float> h2o_frac(nlocs, 0.0);
        geovals.get(h2o_frac, "water_area_fraction");
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          ydiags.get(pred, "brightness_temperature_jacobian_surface_emissivity_" +
                           std::to_string(channels_[jc]));
          for (std::size_t jl = 0; jl < nlocs; ++jl) {
            if (h2o_frac[jl] < 0.99 && std::fabs(pred[jl]) > 0.001) {
              (*preds)[indx+jc][jl] = pred[jl];
            } else {
              (*preds)[indx+jc][jl] = 0.0;
            }
          }
        }
      } else {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = 0.0;
        }
      }
      indx += nchanl;
    } else if (predictors_[n] == "scan_angle_4th_order") {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl)
          (*preds)[indx+jc][jl] = pow(view_angle[jl]*Constants::deg2rad, 4);
      }
      indx += nchanl;
    } else if (predictors_[n] == "scan_angle_3rd_order") {
      if (adp_anglebc_) {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = pow(view_angle[jl]*Constants::deg2rad, 3);
        }
      } else {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = 0.0;
        }
      }
      indx += nchanl;
    } else if (predictors_[n] == "scan_angle_2nd_order") {
      if (adp_anglebc_) {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = pow(view_angle[jl]*Constants::deg2rad, 2);
        }
      } else {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = 0.0;
        }
      }
      indx += nchanl;
    } else if (predictors_[n] == "scan_angle_1st_order") {
      if (adp_anglebc_) {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = view_angle[jl]*Constants::deg2rad;
        }
      } else {
        for (std::size_t jc = 0; jc < nchanl; ++jc) {
          for (std::size_t jl = 0; jl < nlocs; ++jl)
            (*preds)[indx+jc][jl] = 0.0;
        }
      }
      indx += nchanl;
    } else {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl) {
          (*preds)[indx+jc][jl] = 0.0;
        }
      }
      indx += nchanl;
      oops::Log::info() << "predictor: " << predictors_[n] << " is not implemented, "
                        << " ZERO will be filled in" << std::endl;
    }
  }

  oops::Log::trace() << "ObsBiasRadianceGSI::computeObsBiasPredictors done." << std::endl;
}

// -----------------------------------------------------------------------------

double ObsBiasRadianceGSI::norm() const {
  double zz = 0.0;
  for (unsigned int jj = 0; jj < biascoeffs_.size(); ++jj) {
    zz += biascoeffs_[jj] * biascoeffs_[jj];
  }
  if (biascoeffs_.size() > 0) zz = std::sqrt(zz/this->size());
  return zz;
}

// -----------------------------------------------------------------------------

ObsBiasRadianceGSI & ObsBiasRadianceGSI::operator+=(const ObsBiasIncrement & dx) {
  for (unsigned int jj = 0; jj < biascoeffs_.size(); ++jj)
    biascoeffs_[jj] += dx[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasRadianceGSI & ObsBiasRadianceGSI::operator=(const ObsBias & rhs) {
  for (unsigned int jj = 0; jj < biascoeffs_.size(); ++jj)
    biascoeffs_[jj] = rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSI::print(std::ostream & os) const {
  os << "ObsBiasRadianceGSI::print " << sensor_id_ << std::endl;
  std::size_t pred_size = predictors_.size();
  std::size_t jc;
  for (jc = 0; jc < channels_.size(); ++jc) {
    os << "Channel : " << channels_[jc] << std::endl;
    for (std::size_t n = 0; n < pred_size; ++n)
      os << biascoeffs_[jc*pred_size+n] << " ";
    os << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

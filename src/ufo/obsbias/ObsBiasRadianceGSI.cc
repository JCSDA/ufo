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

ObsBiasRadianceGSI::ObsBiasRadianceGSI(const eckit::Configuration & conf)
  : ObsBiasBase(), geovars_(), hdiags_(), tlapmean_(),
    newpc4pred_(false), adp_anglebc_(false), emiss_bc_(false),
    predictors_() {
// Default predictor names
  predictors_ = {"BCPred_Constant_",
                 "BCPred_Scan_Angle_",
                 "BCPred_Cloud_Liquid_Water_",
                 "BCPred_Lapse_Rate_Squared_",
                 "BCPred_Lapse_Rate_",
                 "BCPred_Cosine_Latitude_times_Node_",
                 "BCPred_Sine_Latitude_",
                 "BCPred_Emissivity_",
                 "BCPres_Fourth_Order_View_Angle_",
                 "BCPres_Third_Order_View_Angle_",
                 "BCPres_Second_Order_View_Angle_",
                 "BCPres_First_Order_View_Angle_"
                };
// Parse predictors from the conf
  if (conf.has("ObsBias.predictors")) {
    predictors_.clear();
    predictors_ = conf.getStringVector("ObsBias.predictors");
  }
// GeoVals needed from model
  const std::vector<std::string> vv0{"air_temperature",
                                     "air_pressure",
                                     "air_pressure_levels",
                                     "water_area_fraction",
                                     "land_area_fraction",
                                     "ice_area_fraction",
                                     "surface_snow_area_fraction",
                                     "surface_temperature_where_sea",
                                     "surface_temperature_where_land",
                                     "surface_temperature_where_ice",
                                     "surface_temperature_where_snow"
                                    };
  geovars_.reset(new oops::Variables(vv0));

// Hdiags needed from H diagnostics
  std::vector<std::string> vv{"brightness_temperature_jacobian_surface_emissivity_CH",
                              "optical_thickness_of_atmosphere_layer_CH",
                              "brightness_temperature_CH"
                              };
  hdiags_.reset(new oops::Variables(vv));

// Parse Sensor_ID from the conf
  const eckit::LocalConfiguration obsoprconf(conf, "ObsOperator");
  sensor_id_ = obsoprconf.getString("ObsOptions.Sensor_ID");

// Parse channels from the conf
  const eckit::LocalConfiguration simconf(conf, "ObsSpace.simulate");
  const oops::Variables observed(simconf);
  channels_ = observed.channels();

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
      tlapmean_.push_back(tlap);
      infile >> tsum;
      infile >> ntlapupdate;
      for (std::size_t j=0; j < predictors_.size(); ++j) {
        infile >> par;
        if ( nusis == sensor_id_ &&
            std::find(channels_.begin(), channels_.end(), nuchan) != channels_.end() )
          biascoeffs_.push_back(static_cast<double>(par));
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

void ObsBiasRadianceGSI::computeObsBias(const GeoVaLs & geovals,
                                        ioda::ObsVector & ybias,
                                        const ioda::ObsSpace & odb,
                                        const ObsDiagnostics & ydiags) const {
  const std::size_t npred = predictors_.size();
  const std::size_t nchanl = channels_.size();
  const std::size_t nlocs = ybias.nlocs();
  ASSERT(ybias.nlocs() == odb.nlocs());

  // Allocate predictors
  std::vector<float> preds;

  // Compute the predictors
  this->computeObsBiasPredictors(geovals, odb, ydiags, preds);

  std::size_t index = 0;
  // Loop through each location
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    // Loop through each channel
    std::size_t idx_coeffs = 0;
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      // Linear combination
      ybias[index] = 0.0;
      for (std::size_t n = 0; n < npred; ++n) {
        ybias[index] += biascoeffs_[idx_coeffs] * preds.at(n*nchanl*nlocs+jc*nlocs+jl);
        ++idx_coeffs;
      }
      ++index;
    }
  }
  oops::Log::trace() << "ObsBiasRadianceGSI::computeObsBias done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSI::computeObsBiasPredictors(const GeoVaLs & geovals,
                                                  const ioda::ObsSpace & odb,
                                                  const ObsDiagnostics & ydiags,
                                                  std::vector<float> & preds) const {
  const std::size_t npred = predictors_.size();
  const std::size_t nlocs = odb.nlocs();
  const std::size_t nchanl = channels_.size();

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

  // Temporary storage for one predictor vector (size of nlocs)
  std::vector <float> pred(nlocs, 0.0);

  /*
   * pred(1,:)  = global offset
   */

  if (!newpc4pred_) {
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      for (std::size_t jl = 0; jl < nlocs; ++jl)
        preds.emplace_back(0.01);
    }
  } else {
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      for (std::size_t jl = 0; jl < nlocs; ++jl)
        preds.emplace_back(1.0);
    }
  }

  // Retrieve the sensor_zenith_angle from ObsSpace
  std::vector<float> zasat(nlocs);
  odb.get_db("MetaData", "sensor_zenith_angle", zasat);

  /*
   * pred(2,:)  = zenith angle predictor, is not used and set to zero now
   */

  if (!newpc4pred_) {
    if (ssmi || ssmis || amsre || gmi || amsr2) {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl)
        preds.emplace_back(0.0);
      }
    } else {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl)
          preds.emplace_back(0.10*pow(1.0/cos(zasat[jl]) - 1.0, 2) - .015);
      }
    }
  } else {
    if (adp_anglebc_) {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl)
          preds.emplace_back(0.0);
      }
    } else {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl)
          preds.emplace_back(pow(1.0/cos(zasat[jl]) - 1.0, 2));
      }
    }
  }

  // Comnpute the surface temperature (revisit needed)
  std::vector<float> h2o_frac(nlocs);
  std::vector<float> land_frac(nlocs);
  std::vector<float> ice_frac(nlocs);
  std::vector<float> snow_frac(nlocs);
  std::vector<float> h2o_t(nlocs);
  std::vector<float> land_t(nlocs);
  std::vector<float> ice_t(nlocs);
  std::vector<float> snow_t(nlocs);

  geovals.get(h2o_frac, "water_area_fraction");
  geovals.get(land_frac, "land_area_fraction");
  geovals.get(ice_frac, "ice_area_fraction");
  geovals.get(snow_frac, "surface_snow_area_fraction");

  geovals.get(h2o_t, "surface_temperature_where_sea");
  geovals.get(land_t, "surface_temperature_where_land");
  geovals.get(ice_t, "surface_temperature_where_ice");
  geovals.get(snow_t, "surface_temperature_where_snow");

  std::vector<float> tsavg5(nlocs);
  for (std::size_t jl = 0; jl < nlocs; ++jl)
    tsavg5[jl] = (h2o_frac[jl]*h2o_t[jl] + land_frac[jl]*land_t[jl] +
                  ice_frac[jl]*ice_t[jl] + snow_frac[jl]*snow_t[jl]) /
                 (h2o_frac[jl] + land_frac[jl] + ice_frac[jl] +  snow_frac[jl]);

  // Retrieve the brightness temperature observation
  std::vector<std::vector<float>> tb;
  for (std::size_t jc = 0; jc < nchanl; ++jc) {
    odb.get_db("ObsValue", "brightness_temperature_" +
               std::to_string(channels_[jc]), pred);
    tb.emplace_back(pred);
  }
  // Transpose the array
  std::vector<std::vector<float>> tb_obs;
  std::vector<float> profile;
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
  odb.get_db("MetaData", "scan_position", nadir);

  /*
   * pred(3,:)  = cloud liquid water predictor for clear-sky microwave radiance assimilation
   */

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
        preds.emplace_back(clw[jl]);
    }
  } else {
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      for (std::size_t jl = 0; jl < nlocs; ++jl)
        preds.emplace_back(clw[jl]*cos(zasat[jl]*cos(zasat[jl])));
    }
  }

  // Retrieve the optical_thickness_of_atmosphere_layer from Hdiag
  std::string varname;
  std::vector<std::vector<float>> data;
  for (std::size_t jc = 0; jc < nchanl; ++jc) {
    varname = "optical_thickness_of_atmosphere_layer_" +
              std::to_string(channels_[jc]);
    for (std::size_t js = 0; js < ydiags.nlevs(varname); ++js) {
      ydiags.get(pred, "optical_thickness_of_atmosphere_layer_" +
                 std::to_string(channels_[jc]), js+1);
      data.emplace_back(pred);
    }
  }
  // transpose
  std::vector<std::vector<float>> ptau5;
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    profile.clear();
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      varname = "optical_thickness_of_atmosphere_layer_" +
                std::to_string(channels_[jc]);
      std::size_t nlevs = ydiags.nlevs(varname);
      for (std::size_t js = 0; js < nlevs; ++js)
        profile.emplace_back(data[jc*nlevs+js][jl]);
    }
    ptau5.emplace_back(profile);
  }

  // Retrieve the temperature (revisit the unit, in GSI, it is sensible T ?)
  data.clear();
  std::size_t nsig = geovals.nlevs("air_temperature");
  for (std::size_t js = 0; js < nsig; ++js) {
    geovals.get(pred, "air_temperature", js+1);
    data.emplace_back(pred);
  }
  // transpose
  std::vector<std::vector<float>> tvp;
  for (std::size_t jc = 0; jc < nlocs; ++jc) {
    profile.clear();
    for (std::size_t js = 0; js < nsig; ++js)
      profile.emplace_back(data[js][jc]);
    tvp.emplace_back(profile);
  }

  nsig = geovals.nlevs("air_pressure");

  // pick up tlapmean
  std::vector<float> tlapmean;
  for (std::size_t jc = 0; jc < nchanl; ++jc) {
    tlapmean.emplace_back(tlapmean_[channels_[jc]]);
  }

  std::vector<std::vector<float>> tlap;
  std::vector<float> tlapp(nchanl);

  // Compute tlap
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    calc_tlap_f90(newpc4pred_, nsig, nchanl, ptau5[jl][0], tsavg5[jl],
                  tvp[jl][0], tlapmean[0], tlapp[0]);
    tlap.emplace_back(tlapp);
  }

  /*
   * pred(4,:)  = square of temperature laps rate predictor
   */

  for (std::size_t jc = 0; jc < nchanl; ++jc) {
    for (std::size_t jl = 0; jl < nlocs; ++jl)
      preds.emplace_back(tlap[jl][jc]*tlap[jl][jc]);
  }

  /*
   * pred(5,:)  = temperature laps rate predictor
   */

  for (std::size_t jc = 0; jc < nchanl; ++jc) {
    for (std::size_t jl = 0; jl < nlocs; ++jl)
      preds.emplace_back(tlap[jl][jc]);
  }

  // Retrieve the sensor_azimuth_angle and latitude from ObsSpace
  std::vector<float> cenlat(nlocs);
  std::vector<float> node(nlocs);
  odb.get_db("MetaData", "latitude", cenlat);
  odb.get_db("MetaData", "sensor_azimuth_angle", node);

  /*
   * pred(6,:)  = cosinusoidal predictor for SSMI/S ascending/descending bias
   */

  if (ssmis) {
    if (!newpc4pred_) {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl) {
          if (node[jl] < 1000) {
            preds.emplace_back(ssmis_precond*node[jl]*cos(cenlat[jl]*Constants::deg2rad));
          } else {
            preds.emplace_back(0.0);
          }
        }
      }
    } else {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl)
          preds.emplace_back(node[jl]*cos(cenlat[jl]*Constants::deg2rad));
      }
    }
  } else {
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl)
          preds.emplace_back(0.0);
    }
  }

  /*
   * pred(7,:)  = sinusoidal predictor for SSMI/S
   */

  if (ssmis) {
    if (!newpc4pred_) {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl) {
          if (node[jl] < 1000) {
            preds.emplace_back(ssmis_precond*sin(cenlat[jl]*Constants::deg2rad));
          } else {
            preds.emplace_back(0.0);
          }
        }
      }
    } else {
      for (std::size_t jc = 0; jc < nchanl; ++jc) {
        for (std::size_t jl = 0; jl < nlocs; ++jl)
          preds.emplace_back(sin(cenlat[jl]*Constants::deg2rad));
      }
    }
  } else {
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      for (std::size_t jl = 0; jl < nlocs; ++jl)
        preds.emplace_back(0.0);
    }
  }

  /*
   * pred(8,:)  = emissivity sensitivity predictor for land/sea differences
   */

  if (adp_anglebc_ && emiss_bc_) {
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      ydiags.get(pred, "brightness_temperature_jacobian_surface_emissivity_" +
                       std::to_string(channels_[jc]));
      for (std::size_t jl = 0; jl < nlocs; ++jl) {
        if (h2o_frac[jl] < 0.99 && abs(pred[jl]) > 0.01) {
          preds.emplace_back(pred[jl]);
        } else {
          preds.emplace_back(0.0);
        }
      }
    }
  } else {
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      for (std::size_t jl = 0; jl < nlocs; ++jl)
        preds.emplace_back(0.0);
    }
  }

  // Retrieve the sensor_view_angle from ObsSpace
  odb.get_db("MetaData", "sensor_view_angle", pred);

  /*
   * pred(9,:)  = fourth order polynomial of angle bias correction
   */

  for (std::size_t jc = 0; jc < nchanl; ++jc) {
    for (std::size_t jl = 0; jl < nlocs; ++jl)
      preds.emplace_back(pow(pred[jl]*Constants::deg2rad, 4));
  }

  if (adp_anglebc_) {
    /*
     * pred(10,:) = third order polynomial of angle bias correction
     */

    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      for (std::size_t jl = 0; jl < nlocs; ++jl)
        preds.emplace_back(pow(pred[jl]*Constants::deg2rad, 3));
    }

    /*
     * pred(11,:) = second order polynomial of angle bias correction
     */

    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      for (std::size_t jl = 0; jl < nlocs; ++jl)
        preds.emplace_back(pow(pred[jl]*Constants::deg2rad, 2));
    }

    /*
     * pred(12,:) = first order polynomial of angle bias correction
     */

    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      for (std::size_t jl = 0; jl < nlocs; ++jl)
        preds.emplace_back(pred[jl]*Constants::deg2rad);
    }
  }

  oops::Log::trace() << "ObsBia>sRadianceGSI::computeObsBiasPredictors done." << std::endl;
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

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
#include <set>

#include "ufo/utils/Constants.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ioda/ObsVector.h"

namespace ufo {

static ObsBiasMaker<ObsBiasRadianceGSI> makerBiasRadianceGSI_("GSI");

const std::vector<std::string> ObsBiasRadianceGSI::predictors_
                                = {"BCPred_Constant_",
                                   "BCPred_Scan_Angle_",
                                   "BCPred_Cloud_Liquid_Water_",
                                   "BCPred_Lapse_Rate_Squared_",
                                   "BCPred_Lapse_Rate_",
                                   "BCPred_Cosine_Latitude_times_Node_",
                                   "BCPred_Sine_Latitude_",
                                   "BCPred_Emissivity_"
                                  };
// -----------------------------------------------------------------------------

ObsBiasRadianceGSI::ObsBiasRadianceGSI(const eckit::Configuration & conf)
  : ObsBiasBase(conf), varin_() {
// GeoVals needed from model
  const std::vector<std::string> vv{"air_temperature"};
  varin_.reset(new oops::Variables(vv));

// Parse Sensor_ID from the conf
  const eckit::LocalConfiguration obsoprconf(conf, "ObsOperator");
  sensor_id_ = obsoprconf.getString("ObsOptions.Sensor_ID");

// Parse channels from the conf
  const eckit::LocalConfiguration simconf(conf, "ObsSpace.simulate");
  const oops::Variables observed(simconf);
  channels_ = observed.channels();

  // Read ObsBias parameters first guess if available
  const eckit::LocalConfiguration biasConf(conf, "ObsBias");
  if (biasConf.has("abias_in")) {
    read(biasConf);
  } else {
    biascoeffs_.clear();
    for (std::size_t jc = 0; jc < channels_.size(); ++jc)
      for (std::size_t n = 0; n < predictors_.size() + 4; ++n)
        biascoeffs_.push_back(0.0);
  }

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
  float tslp, tslpm;
  std::size_t level;

  if ( infile.is_open() )
  {
    float par;
    while ( !infile.eof() )
    {
      infile >> ich;
      infile >> nusis;
      infile >> nuchan;
      infile >> tslp;
      infile >> tslpm;
      infile >> level;
      for (std::size_t j=0; j < predictors_.size() + 4; ++j) {
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
                                        const ioda::ObsSpace & odb) const {
  std::size_t npred = predictors_.size();
  std::size_t nchanl = channels_.size();
  std::size_t nlocs = ybias.nlocs();
  ASSERT(ybias.nlocs() == odb.nlocs());

  // retrieve the bias predictors from obs file, it will be computed online
  // From GSI
  // radiance bias correction terms are as follows:
  // pred(1,:)  = global offset
  // pred(2,:)  = zenith angle predictor, is not used and set to zero now
  // pred(3,:)  = cloud liquid water predictor for clear-sky microwave radiance assimilation
  // pred(4,:)  = square of temperature laps rate predictor
  // pred(5,:)  = temperature laps rate predictor
  // pred(6,:)  = cosinusoidal predictor for SSMI/S ascending/descending bias
  // pred(7,:)  = sinusoidal predictor for SSMI/S
  // pred(8,:)  = emissivity sensitivity predictor for land/sea differences
  // pred(9,:)  = fourth order polynomial of angle bias correction
  // pred(10,:) = third order polynomial of angle bias correction
  // pred(11,:) = second order polynomial of angle bias correction
  // pred(12,:) = first order polynomial of angle bias correction

  std::vector< std::vector<float> > preds;
  std::vector <float> pred(nlocs, 0.0);
  for (std::size_t v = 0; v < npred; ++v) {
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      odb.get_db("BiasPred", predictors_[v]+std::to_string(jc+1), nlocs, pred.data());
      preds.push_back(pred);
    }
  }
  std::vector<float> viewing_angle(nlocs);
  odb.get_db("MetaData", "sensor_view_angle", nlocs, viewing_angle.data());
  for (std::size_t jc = 0; jc < nchanl; ++jc) {
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      pred[jl] = pow(viewing_angle[jl]*Constants::deg2rad, 4);
    }
    preds.push_back(pred);
  }
  for (std::size_t jc = 0; jc < nchanl; ++jc) {
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      pred[jl] = pow(viewing_angle[jl]*Constants::deg2rad, 3);
    }
    preds.push_back(pred);
  }
  for (std::size_t jc = 0; jc < nchanl; ++jc) {
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      pred[jl] = pow(viewing_angle[jl]*Constants::deg2rad, 2);
    }
    preds.push_back(pred);
  }

  for (std::size_t jc = 0; jc < nchanl; ++jc) {
    for (std::size_t jl = 0; jl < nlocs; ++jl) {
      pred[jl] = viewing_angle[jl]*Constants::deg2rad;
    }
    preds.push_back(pred);
  }

  std::size_t index = 0;
  // Loop through each location
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    // Loop through each channel
    std::size_t idx_coeffs = 0;
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      // Linear combination
      ybias[index] = 0.0;
      for (std::size_t n = 0; n < npred + 4; ++n) {
        ybias[index] += biascoeffs_[idx_coeffs] * preds[n*nchanl+jc][jl];
        ++idx_coeffs;
      }
      ++index;
    }
  }
  ybias.save("ObsBias");
  oops::Log::trace() << "ObsBiasRadianceGSI::computeObsBias done." << std::endl;
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

void ObsBiasRadianceGSI::print(std::ostream & os) const {
  os << "ObsBiasRadianceGSI::print " << sensor_id_ << std::endl;
  std::size_t pred_size = predictors_.size() + 4;
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

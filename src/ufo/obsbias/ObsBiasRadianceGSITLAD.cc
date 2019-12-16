/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "ufo/obsbias/ObsBiasRadianceGSITLAD.h"

#include "eckit/mpi/Comm.h"

#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ioda/ObsVector.h"

namespace ufo {

static LinearObsBiasMaker<ObsBiasRadianceGSITLAD> makerBiasRadianceGSITLAD_("GSI");

// -----------------------------------------------------------------------------

ObsBiasRadianceGSITLAD::ObsBiasRadianceGSITLAD(const ioda::ObsSpace & odb,
                                               const eckit::Configuration & conf)
  : odb_(odb), predictors_() {
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

// Parse Sensor_ID from the conf
  const eckit::LocalConfiguration obsoprconf(conf, "ObsOperator");
  sensor_id_ = obsoprconf.getString("ObsOptions.Sensor_ID");

// Parse channels from the conf
  const eckit::LocalConfiguration simconf(conf, "ObsSpace.simulate");
  const oops::Variables observed(simconf);
  channels_ = observed.channels();

  // Initialize ObsBias parameters
  biascoeffsinc_.clear();
  for (std::size_t jc = 0; jc < channels_.size(); ++jc)
    for (std::size_t n = 0; n < predictors_.size(); ++n)
      biascoeffsinc_.push_back(0.0);

  oops::Log::trace() << "ObsBiasRadianceGSITLAD created." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSITLAD::read(const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsBiasRadianceGSITLAD::read from not implemented " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSITLAD::write(const eckit::Configuration & conf) const {
  oops::Log::trace() << "ObsBiasRadianceGSITLAD::write to file not implmented " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSITLAD::computeObsBiasTL(const GeoVaLs & geovals,
                                              const ioda::ObsDataVector<float> & preds,
                                              ioda::ObsVector & ybiasinc) const {
  std::size_t npred = predictors_.size();
  std::size_t nchanl = channels_.size();
  std::size_t nlocs = ybiasinc.nlocs();
  ASSERT(ybiasinc.nlocs() == odb_.nlocs());

  std::size_t index = 0;
  // Loop through each locations
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    std::size_t idx_coeffs = 0;
    // Loop through each channel
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      ybiasinc[index] = 0.0;
      // Linear combination
      for (std::size_t n = 0; n < npred; ++n) {
        ybiasinc[index] += biascoeffsinc_[idx_coeffs] * preds[n*nchanl+jc][jl];
        ++idx_coeffs;
      }
      ++index;
    }
  }
  oops::Log::trace() << "ObsBiasRadianceGSITLAD::computeObsBiasTL done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSITLAD::computeObsBiasAD(GeoVaLs & geovals,
                                              const ioda::ObsDataVector<float> & preds,
                                              const ioda::ObsVector & ybiasinc) {
  std::size_t npred = predictors_.size();
  std::size_t nchanl = channels_.size();
  std::size_t nlocs = ybiasinc.nlocs();
  ASSERT(ybiasinc.nlocs() == odb_.nlocs());

  std::size_t index = 0;
  // Loop through each locations
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    std::size_t idx_coeffs = 0;
    // Loop through each channel
    for (std::size_t jc = 0; jc < nchanl; ++jc) {
      // Linear combination
      if (ybiasinc[index] != util::missingValue(ybiasinc[index])) {
        for (std::size_t n = 0; n < npred; ++n) {
          biascoeffsinc_[idx_coeffs] += ybiasinc[index] * preds[n*nchanl+jc][jl];
          ++idx_coeffs;
        }
      } else {
        idx_coeffs += npred;
      }
      ++index;
    }
  }
  // Sum across the processros
  if (odb_.isDistributed())
    odb_.comm().allReduceInPlace(biascoeffsinc_.begin(), biascoeffsinc_.end(), eckit::mpi::sum());
  oops::Log::trace() << "ObsBiasRadianceGSITLAD::computeObsBiasAD done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSITLAD::diff(const ObsBias & b1, const ObsBias & b2) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj]= b1[jj] - b2[jj];
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSITLAD::zero() {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj]= 0.0;
}

// -----------------------------------------------------------------------------

ObsBiasRadianceGSITLAD & ObsBiasRadianceGSITLAD::operator=(const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] = rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasRadianceGSITLAD & ObsBiasRadianceGSITLAD::operator+=(const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] += rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasRadianceGSITLAD & ObsBiasRadianceGSITLAD::operator-=(const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] -= rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasRadianceGSITLAD & ObsBiasRadianceGSITLAD::operator*=(const double fact) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] *= fact;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSITLAD::axpy(const double fact, const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] += fact * rhs[jj];
}

// -----------------------------------------------------------------------------

double ObsBiasRadianceGSITLAD::dot_product_with(const ObsBiasIncrement & rhs) const {
  double zz = 0.0;
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj) {
    zz += biascoeffsinc_[jj] * rhs[jj];
  }
  return zz;
}

// -----------------------------------------------------------------------------

double ObsBiasRadianceGSITLAD::norm() const {
  double zz = 0.0;
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj) {
    zz += biascoeffsinc_[jj] * biascoeffsinc_[jj];
  }
  if (biascoeffsinc_.size() > 0) zz = std::sqrt(zz/biascoeffsinc_.size());
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBiasRadianceGSITLAD::print(std::ostream & os) const {
  os << "ObsBiasRadianceGSITLAD::print " << sensor_id_ << std::endl;
  std::size_t pred_size = predictors_.size();
  std::size_t jc;
  for (jc = 0; jc < channels_.size(); ++jc) {
    os << "Channel : " << channels_[jc] << std::endl;
    for (std::size_t n = 0; n < pred_size; ++n)
      os << biascoeffsinc_[jc*pred_size+n] << " ";
    os << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

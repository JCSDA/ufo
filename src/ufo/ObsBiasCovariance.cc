/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <Eigen/Core>
#include <cmath>
#include <fstream>
#include <memory>
#include <random>
#include <set>

#include "ufo/ObsBiasCovariance.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/predictors/PredictorBase.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasCovariance::ObsBiasCovariance(ioda::ObsSpace & odb, const eckit::Configuration & conf)
  : conf_(conf), odb_(odb), prednames_(0), jobs_(0), variances_(0), preconditioner_(0),
    ht_rinv_h_(0), obs_num_(0), analysis_variances_(0), minimal_required_obs_number_(0) {
  oops::Log::trace() << "ObsBiasCovariance::Constructor starting" << std::endl;

  // Get the jobs(channels)
  if (conf_.has("obs bias.jobs")) {
    const std::set<int> jobs = oops::parseIntSet(conf_.getString("obs bias.jobs"));
    jobs_.assign(jobs.begin(), jobs.end());
  }

  // Predictor factory
  if (conf_.has("obs bias.predictors")) {
    std::vector<eckit::LocalConfiguration> confs;
    conf_.get("obs bias.predictors", confs);
    for (std::size_t j = 0; j < confs.size(); ++j) {
      std::shared_ptr<PredictorBase> pred(PredictorFactory::create(confs[j], jobs_));
      prednames_.push_back(pred->name());
    }
  }

  if (prednames_.size()*jobs_.size() > 0) {
    // Get the minimal required filtered obs number
    minimal_required_obs_number_ =
      conf_.getUnsigned("obs bias.covariance.minimal required obs number");

    // Override the variance range if provided
    if (conf_.has("obs bias.covariance.variance range")) {
      const std::vector<double>
        range = conf_.getDoubleVector("obs bias.covariance.variance range");
      ASSERT(range.size() == 2);
      smallest_variance_ = range[0];
      largest_variance_ = range[1];
    }

    // Override the preconditioning step size  if provided
    if (conf_.has("obs bias.covariance.step size")) {
      step_size_ = conf_.getDouble("obs bias.covariance.step size");
    }

    // Override the largest analysis variance if provided
    if (conf_.has("obs bias.covariance.largest analysis variance")) {
      largest_analysis_variance_ =
        conf_.getDouble("obs bias.covariance.largest analysis variance");
    }

    // Initialize the variances to upper limit
    variances_.resize(prednames_.size() * jobs_.size());
    std::fill(variances_.begin(), variances_.end(), largest_variance_);

    // Initialize the hessian contribution to zero
    ht_rinv_h_.resize(prednames_.size() * jobs_.size());
    std::fill(ht_rinv_h_.begin(), ht_rinv_h_.end(), 0.0);

    // Initialize the preconditioner to default step size
    preconditioner_.resize(prednames_.size() * jobs_.size());
    std::fill(preconditioner_.begin(), preconditioner_.end(), step_size_);

    // Initialize obs_num_ to ZERO
    obs_num_.resize(jobs_.size());
    std::fill(obs_num_.begin(), obs_num_.end(), 0);

    // Initialize analysis error variances to the upper limit
    analysis_variances_.resize(prednames_.size() * jobs_.size());
    std::fill(analysis_variances_.begin(), analysis_variances_.end(), largest_variance_);

    // Initializes from given prior
    if (conf_.has("obs bias.covariance.prior")) {
      // Get default inflation ratio
      const double inflation_ratio =
        conf.getDouble("obs bias.covariance.prior.inflation.ratio");

      // Check the large inflation ratio when obs number < minimal_required_obs_number
      const double large_inflation_ratio =
        conf.getDouble("obs bias.covariance.prior.inflation.ratio for small dataset");

      // read in Variances prior (analysis_variances_) and number of obs. (obs_num_)
      // from previous cycle
      this->read(conf_);

      // set variances for bias predictor coeff. based on diagonal info
      // of previous analysis error variance
      std::size_t ii;
      for (std::size_t j = 0; j < jobs_.size(); ++j) {
        const double inflation = (obs_num_[j] <= minimal_required_obs_number_) ?
                                 large_inflation_ratio : inflation_ratio;
        for (std::size_t p = 0; p < prednames_.size(); ++p) {
          ii = j*prednames_.size() + p;
          if (inflation > inflation_ratio)
            analysis_variances_[ii] = inflation * analysis_variances_[ii] + smallest_variance_;
          variances_[ii] = inflation * analysis_variances_[ii] + smallest_variance_;
          if (variances_[ii] > largest_variance_) variances_[ii] = largest_variance_;
          if (analysis_variances_[ii] > largest_analysis_variance_)
            analysis_variances_[ii] = largest_analysis_variance_;
        }
      }
    }
  }

  oops::Log::trace() << "ObsBiasCovariance::Constructor is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::read(const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsBiasCovariance::read from file " << std::endl;

  // Default predictor names from GSI
  const std::vector<std::string> gsi_predictors = {"constant",
                                                   "zenith_angle",
                                                   "cloud_liquid_water",
                                                   "lapse_rate_order_2",
                                                   "lapse_rate",
                                                   "cosine_of_latitude_times_orbit_node",
                                                   "sine_of_latitude",
                                                   "emissivity",
                                                   "scan_angle_order_4",
                                                   "scan_angle_order_3",
                                                   "scan_angle_order_2",
                                                   "scan_angle"
                                                   };

  const std::string sensor = conf.getString("obs bias.sensor");

  if (conf.has("obs bias.covariance.prior.datain")) {
    const std::string filename = conf.getString("obs bias.covariance.prior.datain");
    std::ifstream infile(filename);

    if (infile.is_open())
    {
      int ich;            //  sequential number
      std::string nusis;  //  sensor/instrument/satellite
      int nuchan;         //  channel number
      float number;       //  QCed obs number from previous cycle

      float par;
      while (infile >> ich)
      {
        infile >> nusis;
        infile >> nuchan;
        infile >> number;
        if (nusis == sensor) {
          auto ijob = std::find(jobs_.begin(), jobs_.end(), nuchan);
          if (ijob != jobs_.end()) {
            int j = std::distance(jobs_.begin(), ijob);
            obs_num_[j] = static_cast<int>(number);

            for (auto & item : gsi_predictors) {
              infile >> par;
              auto ipred = std::find(prednames_.begin(), prednames_.end(), item);
              if (ipred != prednames_.end()) {
                int p = std::distance(prednames_.begin(), ipred);
                analysis_variances_.at(j*prednames_.size() + p) = static_cast<double>(par);
              }
            }
          } else {
            for (auto & item : gsi_predictors)
              infile >> par;
          }
        } else {
          for (auto & item : gsi_predictors)
            infile >> par;
        }
      }
      infile.close();
      oops::Log::trace() << "ObsBiasCovariance::read from prior file: "
                         << filename << " Done " << std::endl;
    } else {
      oops::Log::error() << "Unable to open file : " << filename << std::endl;
      ABORT("Unable to open bias correction coeffs variance prior file ");
    }
  }

  oops::Log::trace() << "ObsBiasCovariance::read is done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::write(const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsBiasCovariance::write to file " << std::endl;
  oops::Log::trace() << "ObsBiasCovariance::write is not implemented " << std::endl;
  oops::Log::trace() << "ObsBiasCovariance::write is done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::linearize(const ObsBias & bias, const eckit::Configuration & innerConf) {
  oops::Log::trace() << "ObsBiasCovariance::linearize starts" << std::endl;
  if (bias) {
    // Retrieve the QC flags and do statistics from second outer loop
    const int jouter = innerConf.getInt("iteration");
    if (jouter >= 1) {
      // Reset the number of filtered Obs.
      std::fill(obs_num_.begin(), obs_num_.end(), 0);

      // Retrieve the QC flags of previous outer loop and get the number of effective obs.
      const std::string group_name = "EffectiveQC" + std::to_string(jouter-1);
      const std::vector<std::string> vars = odb_.obsvariables().variables();
      std::vector<int> qc_flags(odb_.nlocs(), 999);
      for (std::size_t j = 0; j < vars.size(); ++j) {
        if (odb_.has(group_name , vars[j])) {
          odb_.get_db(group_name, vars[j],  qc_flags);
          obs_num_[j] = std::count(qc_flags.begin(), qc_flags.end(), 0);
        } else {
          throw eckit::UserError("Unable to find QC flags : " + vars[j] + "@" + group_name);
        }
      }

      // Sum across the processors
      odb_.distribution().sum(obs_num_);

      const float missing = util::missingValue(missing);

      // compute the hessian contribution from Jo bias terms channel by channel
      // retrieve the effective error (after QC) for this channel
      ioda::ObsVector r_inv(odb_, "EffectiveError", true);

      // compute \mathrm{R}^{-1}
      std::size_t nvars = r_inv.nvars();
      for (size_t vv = 0; vv < nvars; ++vv) {
        for (size_t ii = 0; ii < r_inv.nlocs(); ++ii) {
          if (r_inv[ii*nvars + vv] != missing) {
            r_inv[ii*nvars + vv] = 1.0f / pow(r_inv[ii*nvars + vv], 2);
          } else {
            r_inv[ii*nvars + vv] = 0.0f;
          }
        }
      }

      // Sum the total number of effective obs. across tasks
      odb_.distribution().sum(obs_num_);

      // compute \mathrm{H}_\beta^\intercal \mathrm{R}^{-1} \mathrm{H}_\beta
      // -----------------------------------------
      std::fill(ht_rinv_h_.begin(), ht_rinv_h_.end(), 0.0);
      for (std::size_t p = 0; p < prednames_.size(); ++p) {
        // retrieve the predictors
        const ioda::ObsVector predx(odb_, prednames_[p] + "Predictor", true);

        // for each variable
        ASSERT(r_inv.nlocs() == predx.nlocs());
        std::size_t nvars = predx.nvars();
        // only keep the diagnoal
        for (size_t vv = 0; vv < nvars; ++vv) {
          for (size_t ii = 0; ii < predx.nlocs(); ++ii)
            ht_rinv_h_[vv*prednames_.size() + p] +=
              pow(predx[ii*nvars + vv], 2) * r_inv[ii*nvars + vv];
        }
      }

      // Sum the hessian contributions across the tasks
      odb_.distribution().sum(ht_rinv_h_);
    }

    // reset variances for bias predictor coeff. based on current data count
    for (std::size_t j = 0; j < jobs_.size(); ++j) {
      if (obs_num_[j] <= minimal_required_obs_number_) {
        for (std::size_t p = 0; p < prednames_.size(); ++p)
          variances_[j*prednames_.size() + p] = smallest_variance_;
      }
    }

    // set a coeff. factor for variances of control variables
    std::size_t index;
    for (std::size_t j = 0; j < jobs_.size(); ++j) {
      for (std::size_t p = 0; p < prednames_.size(); ++p) {
        index = j*prednames_.size() + p;
        preconditioner_[index] = step_size_;
        // L = \mathrm{A}^{-1}
        if (obs_num_[j] > 0)
          preconditioner_[index] = 1.0 / (1.0 / variances_[index] + ht_rinv_h_[index]);
        if (obs_num_[j] > minimal_required_obs_number_) {
          if (ht_rinv_h_[index] > 0.0) {
            analysis_variances_[index] = 1.0 / (1.0 / variances_[index] + ht_rinv_h_[index]);
          } else {
            analysis_variances_[index] = largest_analysis_variance_;
          }
        }
      }
    }
  }
  oops::Log::trace() << "ObsBiasCovariance::linearize is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::multiply(const ObsBiasIncrement & dx1,
                                 ObsBiasIncrement & dx2) const {
  oops::Log::trace() << "ObsBiasCovariance::multiply starts" << std::endl;

  for (std::size_t jj = 0; jj < variances_.size(); ++jj) {
    dx2[jj] = dx1[jj] * variances_[jj];
  }

  oops::Log::trace() << "ObsBiasCovariance::multiply is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::inverseMultiply(const ObsBiasIncrement & dx1,
                                        ObsBiasIncrement & dx2) const {
  oops::Log::trace() << "ObsBiasCovariance::inverseMultiply starts" << std::endl;

  for (std::size_t jj = 0; jj < variances_.size(); ++jj) {
    dx2[jj] = dx1[jj] / variances_[jj];
  }

  oops::Log::trace() << "ObsBiasCovariance::inverseMultiply is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::randomize(ObsBiasIncrement & dx) const {
  oops::Log::trace() << "ObsBiasCovariance::randomize starts" << std::endl;
  if (dx) {
    static util::NormalDistribution<double> dist(variances_.size());
    for (std::size_t jj = 0; jj < variances_.size(); ++jj) {
      dx[jj] = dist[jj] * std::sqrt(variances_[jj]);
    }
  }
  oops::Log::trace() << "ObsBiasCovariance::randomize is done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <fstream>
#include <memory>
#include <random>
#include <set>

#include "ufo/ObsBiasCovariance.h"

#include "ioda/ObsSpace.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/predictors/PredictorBase.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasCovariance::ObsBiasCovariance(const ioda::ObsSpace & odb, const eckit::Configuration & conf)
  : conf_(conf), odb_(odb), prednames_(0), jobs_(0), variances_(0),
    obs_num_prior_(0), variances_prior_(0), minimal_required_obs_number_(0) {
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

    // Initialize variances to upper limit
    variances_.resize(prednames_.size() * jobs_.size());
    std::fill(variances_.begin(), variances_.end(), largest_variance_);

    // Initialize obs_num_prior to ZERO
    obs_num_prior_.resize(jobs_.size());
    std::fill(obs_num_prior_.begin(), obs_num_prior_.end(), 0);

    // Initialize analysis error variances to the upper limit
    variances_prior_.resize(prednames_.size() * jobs_.size());
    std::fill(variances_prior_.begin(), variances_prior_.end(), largest_variance_);

    // Initializes from given prior
    if (conf_.has("obs bias.covariance.prior")) {
      // Get default inflation ratio
      const double inflation_ratio =
        conf.getDouble("obs bias.covariance.prior.inflation.ratio");

      // Check the large inflation ratio when obs number < minimal_required_obs_number
      const double large_inflation_ratio =
        conf.getDouble("obs bias.covariance.prior.inflation.ratio for small dataset");

      // read in Variances prior (variances_prior_) and number of obs. (obs_num_prior_)
      // from previous cycle
      this->read(conf_);

      // set variances for bias predictor coeff. based on diagonal info
      // of previous analysis error variance
      std::size_t ii;
      for (std::size_t j = 0; j < jobs_.size(); ++j) {
        const double inflation = (obs_num_prior_[j] <= minimal_required_obs_number_) ?
                                 large_inflation_ratio : inflation_ratio;
        for (std::size_t p = 0; p < prednames_.size(); ++p) {
          ii = j*prednames_.size() + p;
          variances_[ii] = inflation * variances_prior_[ii] + smallest_variance_;
          if (variances_[ii] > largest_variance_) variances_[ii] = largest_variance_;
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
      std::map<std::string, double> elem;
      while (infile >> ich)
      {
        infile >> nusis;
        infile >> nuchan;
        infile >> number;
        if (nusis == sensor) {
          auto ijob = std::find(jobs_.begin(), jobs_.end(), nuchan);
          if (ijob != jobs_.end()) {
            int j = std::distance(jobs_.begin(), ijob);
            obs_num_prior_[j] = static_cast<int>(number);

            for (auto & item : gsi_predictors) {
              infile >> par;
              auto ipred = std::find(prednames_.begin(), prednames_.end(), item);
              if (ipred != prednames_.end()) {
                int p = std::distance(prednames_.begin(), ipred);
                variances_prior_.at(j*prednames_.size() + p) = static_cast<double>(par);
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
    if ((jouter - 1) >= 0) {
      // Reset the number of filtered Obs.
      std::fill(obs_num_prior_.begin(), obs_num_prior_.end(), 0);

      // Retrieve the QC flags of previous outer loop and get the number of effective obs.
      const std::string group_name = "EffectiveQC" + std::to_string(jouter-1);
      const std::vector<std::string> vars = odb_.obsvariables().variables();
      std::vector<int> qc_flags(odb_.nlocs(), 999);
      for (std::size_t j = 0; j < vars.size(); ++j) {
        if (odb_.has(group_name , vars[j])) {
          odb_.get_db(group_name, vars[j],  qc_flags);
          obs_num_prior_[j] = std::count(qc_flags.begin(), qc_flags.end(), 0);
        } else {
          throw eckit::UserError("Unable to find QC flags : " + vars[j] + "@" + group_name);
        }
      }
      // Sum across the processros
      if (odb_.isDistributed())
        odb_.comm().allReduceInPlace(obs_num_prior_.begin(), obs_num_prior_.end(),
                                     eckit::mpi::sum());
    }

    for (std::size_t j = 0; j < jobs_.size(); ++j) {
      if (obs_num_prior_[j] <= minimal_required_obs_number_) {
        for (std::size_t p = 0; p < prednames_.size(); ++p)
          variances_[j*prednames_.size() + p] = smallest_variance_;
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
    dx2[jj] = dx1[jj] * variances_.at(jj);
  }

  oops::Log::trace() << "ObsBiasCovariance::multiply is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::inverseMultiply(const ObsBiasIncrement & dx1,
                                        ObsBiasIncrement & dx2) const {
  oops::Log::trace() << "ObsBiasCovariance::inverseMultiply starts" << std::endl;

  for (std::size_t jj = 0; jj < variances_.size(); ++jj) {
    dx2[jj] = dx1[jj] / variances_.at(jj);
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

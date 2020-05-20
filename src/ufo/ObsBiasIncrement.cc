/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <set>

#include "ufo/ObsBiasIncrement.h"

#include "ioda/ObsSpace.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ioda::ObsSpace & odb, const eckit::Configuration & conf)
  : biasbase_(), predbases_(0), jobs_(0), odb_(odb), conf_(conf) {
  oops::Log::trace() << "ObsBiasIncrement::create starting." << std::endl;

  /// Get the jobs(channels)
  if (conf_.has("ObsBias.jobs")) {
    const std::set<int> jobs = oops::parseIntSet(conf_.getString("ObsBias.jobs"));
    jobs_.assign(jobs.begin(), jobs.end());
  }

  /// Predictor factory
  if (conf_.has("ObsBias.predictors")) {
    std::vector<eckit::LocalConfiguration> confs;
    conf_.get("ObsBias.predictors", confs);
    typedef std::unique_ptr<PredictorBase> predictor;
    for (std::size_t j = 0; j < confs.size(); ++j) {
      predbases_.push_back(predictor(PredictorFactory::create(confs[j], jobs_)));
      prednames_.push_back(predbases_[j]->name());
    }
  }

  /// Bias model factory
  biasbase_.reset(LinearObsBiasFactory::create(odb_, conf_, prednames_, jobs_));

  oops::Log::trace() << "ObsBiasIncrement::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other, const bool copy)
  : odb_(other.odb_), conf_(other.conf_), biasbase_(),
    predbases_(other.predbases_), prednames_(other.prednames_), jobs_(other.jobs_) {
  oops::Log::trace() << "ObsBiasIncrement::copy ctor starting" << std::endl;

  /// Create a new bias model object
  biasbase_.reset(LinearObsBiasFactory::create(odb_, conf_, prednames_, jobs_));

  /// Copy the bias model coeff data
  if (copy && biasbase_) *biasbase_ = other;

  oops::Log::trace() << "ObsBiasIncrement::copy ctor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other,
                                   const eckit::Configuration & conf)
  : odb_(other.odb_), conf_(conf), biasbase_(), predbases_(), prednames_(), jobs_() {
  oops::Log::trace() << "ObsBiasIncrement::copy ctor starting." << std::endl;
  /// Get the jobs(channels)
  if (conf_.has("ObsBias.jobs")) {
    const std::set<int> jobs = oops::parseIntSet(conf_.getString("ObsBias.jobs"));
    jobs_.assign(jobs.begin(), jobs.end());
  }

  /// Predictor factory
  if (conf_.has("ObsBias.predictors")) {
    std::vector<eckit::LocalConfiguration> confs;
    conf_.get("ObsBias.predictors", confs);
    typedef std::unique_ptr<PredictorBase> predictor;
    for (std::size_t j = 0; j < confs.size(); ++j) {
      predbases_.push_back(predictor(PredictorFactory::create(confs[j], jobs_)));
      prednames_.push_back(predbases_[j]->name());
    }
  }

  /// Bias model factory
  biasbase_.reset(LinearObsBiasFactory::create(odb_, conf_, prednames_, jobs_));

  /// Copy the data
  if (biasbase_) *biasbase_ = other;

  oops::Log::trace() << "ObsBiasIncrement::copy ctor done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::diff(const ObsBias & b1, const ObsBias & b2) {
  if (biasbase_) biasbase_->diff(b1, b2);
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::zero() {
  if (biasbase_) biasbase_->zero();
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator=(const ObsBiasIncrement & rhs) {
  if (rhs) {
    predbases_.clear();
    jobs_.clear();

    predbases_ = rhs.predbases_;
    jobs_      = rhs.jobs_;
    *biasbase_ = rhs;
  }
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator+=(const ObsBiasIncrement & rhs) {
  if (biasbase_) *biasbase_ += rhs;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator-=(const ObsBiasIncrement & rhs) {
  if (biasbase_) *biasbase_ -= rhs;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator*=(const double fact) {
  if (biasbase_) *biasbase_ *= fact;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::axpy(const double fact, const ObsBiasIncrement & rhs) {
  if (biasbase_) biasbase_->axpy(fact, rhs);
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::dot_product_with(const ObsBiasIncrement & rhs) const {
  double zz = 0.0;
  if (biasbase_) zz = biasbase_->dot_product_with(rhs);
  return zz;
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::norm() const {
  double zz = 0.0;
  if (biasbase_) zz = biasbase_->norm();
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::computeObsBiasTL(const GeoVaLs & geovals,
                                        const Eigen::MatrixXd & preds,
                                        ioda::ObsVector & ybiasinc) const {
  if (biasbase_) {
    biasbase_->computeObsBiasTL(geovals, preds, ybiasinc);
  }
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::computeObsBiasAD(GeoVaLs & geovals,
                                        const Eigen::MatrixXd & preds,
                                        const ioda::ObsVector & ybiasinc) {
  if (biasbase_) {
    biasbase_->computeObsBiasAD(geovals, preds, ybiasinc);
  }
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::print(std::ostream & os) const {
  if (biasbase_) os << *biasbase_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

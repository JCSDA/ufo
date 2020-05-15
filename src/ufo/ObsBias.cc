/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/ObsBias.h"

#include <Eigen/Core>
#include <set>

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

#include "ufo/ObsBiasIncrement.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ioda::ObsSpace & odb, const eckit::Configuration & conf)
  : biasbase_(), predbases_(0), jobs_(0), odb_(odb), conf_(conf) {
  oops::Log::trace() << "ObsBias::create starting." << std::endl;

  /// Get the jobs(channels)
  if (conf_.has("ObsBias.jobs")) {
    const std::set<int> jobs = oops::parseIntSet(conf_.getString("ObsBias.jobs"));
    jobs_.assign(jobs.begin(), jobs.end());
  }

  // Predictor factory
  if (conf_.has("ObsBias.predictors")) {
    std::vector<eckit::LocalConfiguration> confs;
    conf_.get("ObsBias.predictors", confs);
    for (std::size_t j = 0; j < confs.size(); ++j) {
      std::shared_ptr<PredictorBase> pred(PredictorFactory::create(confs[j]));
      predbases_.push_back(pred);
      prednames_.push_back(pred->name());
      geovars_ += oops::Variables(pred->requiredGeovars());
      for (auto  item : pred->requiredHdiagnostics()) {
        if (item.size() > 3 && item.substr(item.size() - 3) == "_CH") {
          hdiags_ += oops::Variables({item.substr(0, item.length() - 3)}, jobs_);
        } else {
          hdiags_ += oops::Variables({item});
        }
      }
    }
  }

  // bias model factory
  biasbase_.reset(ObsBiasFactory::create(conf_, prednames_, jobs_));

  /// read or initialize bias coefficients
  this->read(conf);

  oops::Log::trace() << "ObsBias::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : odb_(other.odb_), conf_(other.conf_), biasbase_(), predbases_(other.predbases_),
    prednames_(other.prednames_), jobs_(other.jobs_),
    geovars_(other.geovars_), hdiags_(other.hdiags_) {
  oops::Log::trace() << "ObsBias::copy ctor starting." << std::endl;

  // Creat a new bias model object
  biasbase_.reset(ObsBiasFactory::create(conf_, prednames_, jobs_));

  // copy the bias model coeff data
  if (copy && biasbase_) *biasbase_ = *other.biasbase_;

  oops::Log::trace() << "ObsBias::copy ctor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator+=(const ObsBiasIncrement & dx) {
  if (biasbase_) *biasbase_+=dx;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator=(const ObsBias & rhs) {
  if (rhs) {
    ASSERT(biasbase_);
    conf_ = rhs.conf_;
    *biasbase_ = *rhs.biasbase_;
    predbases_ = rhs.predbases_;
    prednames_ = rhs.prednames_;
    jobs_ = rhs.jobs_;
    geovars_ = rhs.geovars_;
    hdiags_ = rhs.hdiags_;
  }
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBias::read(const eckit::Configuration & conf) {
  if (biasbase_) {
    std::string sensor = conf.getString("ObsBias.sensor");
    biasbase_->read(sensor);
  }
}

// -----------------------------------------------------------------------------

void ObsBias::write(const eckit::Configuration & conf) const {
  if (biasbase_) biasbase_->write(conf);
}

// -----------------------------------------------------------------------------

void ObsBias::computeObsBias(ioda::ObsVector & ybias,
                             const Eigen::MatrixXd & predData) const {
  if (biasbase_) biasbase_->computeObsBias(ybias, predData);
}

// -----------------------------------------------------------------------------
Eigen::MatrixXd ObsBias::computePredictors(const GeoVaLs & geovals,
                                           const ObsDiagnostics & ydiags) const {
  const std::size_t nlocs  = odb_.nlocs();
  const std::size_t npreds = predbases_.size();
  const std::size_t njobs  = jobs_.size();

  Eigen::MatrixXd predData(npreds*njobs, nlocs);

  if (biasbase_) {
    // temporary workspace
    Eigen::MatrixXd tmp(njobs, nlocs);

    for (std::size_t r = 0; r < npreds; ++r) {
      // initialize with zero
      tmp.setConstant(0.0);

      // calculate the predictor
      predbases_[r]->compute(odb_, geovals, ydiags, jobs_, tmp);

      // save
      for (std::size_t i = 0; i < njobs; ++i) {
        predData.row(r+i*npreds) = tmp.row(i);
      }
    }
  }

  oops::Log::trace() << "ObsBias::computePredictors done." << std::endl;
  return predData;
}

// -----------------------------------------------------------------------------

void ObsBias::saveObsBiasTerms(ioda::ObsSpace & odb,
                               const std::string & group,
                               const Eigen::MatrixXd & predData) const {
  if (biasbase_) {
    biasbase_->saveObsBiasTerms(odb, group, predData);
  }
}

// -----------------------------------------------------------------------------

double ObsBias::norm() const {
  double zz = 0.0;
  if (biasbase_) zz = biasbase_->norm();
  return zz;
}

// -----------------------------------------------------------------------------

std::size_t ObsBias::size() const {
  std::size_t zz = 0;
  if (biasbase_) zz = biasbase_->size();
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBias::print(std::ostream & os) const {
  if (biasbase_) os << *biasbase_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

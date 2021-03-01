/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBiasIncrement.h"

#include <iomanip>
#include <set>

#include "eckit/mpi/Comm.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ioda::ObsSpace & odb, const eckit::Configuration & conf)
  : odb_(odb), jobs_(0) {
  oops::Log::trace() << "ObsBiasIncrement::create starting." << std::endl;

  // Get the jobs(channels)
  if (conf.has("jobs")) {
    const std::set<int> jobs = oops::parseIntSet(conf.getString("jobs"));
    jobs_.assign(jobs.begin(), jobs.end());
  }

  // Predictor factory
  if (conf.has("predictors")) {
    std::vector<eckit::LocalConfiguration> confs;
    conf.get("predictors", confs);
    typedef std::unique_ptr<PredictorBase> predictor;
    for (std::size_t j = 0; j < confs.size(); ++j) {
      std::unique_ptr<PredictorBase> predictor(PredictorFactory::create(confs[j], jobs_));
      prednames_.push_back(predictor->name());
    }
  }

  // initialize bias coefficient perturbations
  biascoeffsinc_.resize(prednames_.size()*jobs_.size());
  std::fill(biascoeffsinc_.begin(), biascoeffsinc_.end(), 0.0);

  oops::Log::trace() << "ObsBiasIncrement::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other, const bool copy)
  : odb_(other.odb_),
    prednames_(other.prednames_), jobs_(other.jobs_) {
  oops::Log::trace() << "ObsBiasIncrement::copy ctor starting" << std::endl;

  // initialize bias coefficient perturbations
  biascoeffsinc_.resize(prednames_.size()*jobs_.size());

  // Copy the bias coefficients data, or fill in with zeros
  if (copy) {
    biascoeffsinc_ = other.biascoeffsinc_;
  } else {
    std::fill(biascoeffsinc_.begin(), biascoeffsinc_.end(), 0.0);
  }

  oops::Log::trace() << "ObsBiasIncrement::copy ctor done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::diff(const ObsBias & b1, const ObsBias & b2) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj]= b1[jj] - b2[jj];
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::zero() {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj]= 0.0;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator=(const ObsBiasIncrement & rhs) {
  if (rhs) {
    prednames_     = rhs.prednames_;
    jobs_          = rhs.jobs_;
    biascoeffsinc_ = rhs.biascoeffsinc_;
  }
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator+=(const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] += rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator-=(const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] -= rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator*=(const double fact) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] *= fact;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::axpy(const double fact, const ObsBiasIncrement & rhs) {
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj)
    biascoeffsinc_[jj] += fact * rhs[jj];
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::dot_product_with(const ObsBiasIncrement & rhs) const {
  double zz = 0.0;
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj) {
    zz += biascoeffsinc_[jj] * rhs[jj];
  }
  return zz;
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::norm() const {
  double zz = 0.0;
  for (std::size_t jj = 0; jj < biascoeffsinc_.size(); ++jj) {
    zz += biascoeffsinc_[jj] * biascoeffsinc_[jj];
  }
  if (biascoeffsinc_.size() > 0) zz = std::sqrt(zz/biascoeffsinc_.size());
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::allSumInPlace() {
  odb_.comm().allReduceInPlace(biascoeffsinc_.begin(), biascoeffsinc_.end(), eckit::mpi::sum());
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::print(std::ostream & os) const {
  if (this->serialSize() > 0) {
    // map bias coeffs to eigen matrix (writable)
    Eigen::Map<const Eigen::MatrixXd>
      coeffs(biascoeffsinc_.data(), prednames_.size(), jobs_.size());

    os << "ObsBiasIncrement::print " << std::endl;
    os << "---------------------------------------------------------------" << std::endl;
    auto njobs = jobs_.size();
    for (std::size_t p = 0; p < prednames_.size(); ++p) {
      os << std::fixed << std::setw(20) << prednames_[p]
         << ":  Min= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).minCoeff()
         << ",  Max= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).maxCoeff()
         << ",  Norm= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).norm()
         << std::endl;
    }
    os << "---------------------------------------------------------------" << std::endl;
    os << "ObsBiasIncrement::print done" << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

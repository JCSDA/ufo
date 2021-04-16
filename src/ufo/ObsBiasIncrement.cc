/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBiasIncrement.h"

#include <iomanip>
#include <memory>

#include "eckit/mpi/Comm.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ioda::ObsSpace & odb,
                                   const Parameters_ & params)
  : odb_(odb), vars_(odb.obsvariables()) {
  oops::Log::trace() << "ObsBiasIncrement::create starting." << std::endl;

  // Predictor factory
  for (const eckit::LocalConfiguration &conf : params.variationalBC.value().predictors.value()) {
    std::unique_ptr<PredictorBase> predictor(PredictorFactory::create(conf, vars_));
    prednames_.push_back(predictor->name());
  }

  // initialize bias coefficient perturbations
  biascoeffsinc_ = Eigen::VectorXd::Zero(prednames_.size() * vars_.size());

  oops::Log::trace() << "ObsBiasIncrement::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other, const bool copy)
  : odb_(other.odb_),
    prednames_(other.prednames_), vars_(other.vars_) {
  oops::Log::trace() << "ObsBiasIncrement::copy ctor starting" << std::endl;

  // Copy the bias coefficients data, or fill in with zeros
  if (copy) {
    biascoeffsinc_ = other.biascoeffsinc_;
  } else {
    biascoeffsinc_ = Eigen::VectorXd::Zero(prednames_.size() * vars_.size());
  }

  oops::Log::trace() << "ObsBiasIncrement::copy ctor done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::diff(const ObsBias & b1, const ObsBias & b2) {
  biascoeffsinc_ = b1.data() - b2.data();
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::zero() {
  biascoeffsinc_ = Eigen::VectorXd::Zero(prednames_.size() * vars_.size());
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator=(const ObsBiasIncrement & rhs) {
  if (rhs) {
    prednames_     = rhs.prednames_;
    vars_          = rhs.vars_;
    biascoeffsinc_ = rhs.biascoeffsinc_;
  }
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator+=(const ObsBiasIncrement & rhs) {
  biascoeffsinc_ += rhs.biascoeffsinc_;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator-=(const ObsBiasIncrement & rhs) {
  biascoeffsinc_ -= rhs.biascoeffsinc_;
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement & ObsBiasIncrement::operator*=(const double fact) {
  biascoeffsinc_ *= fact;
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::axpy(const double fact, const ObsBiasIncrement & rhs) {
  biascoeffsinc_ += fact * rhs.biascoeffsinc_;
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::dot_product_with(const ObsBiasIncrement & rhs) const {
  return biascoeffsinc_.dot(rhs.biascoeffsinc_);
}

// -----------------------------------------------------------------------------

double ObsBiasIncrement::norm() const {
  double zz = 0.0;
  if (biascoeffsinc_.size() > 0) {
    zz = biascoeffsinc_.norm()/std::sqrt(biascoeffsinc_.size());
  }
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::allSumInPlace() {
  std::vector<double> buffer(biascoeffsinc_.data(),
                             biascoeffsinc_.data() + biascoeffsinc_.size());
  odb_.distribution().allReduceInPlace(buffer, eckit::mpi::sum());
  biascoeffsinc_ = Eigen::Map<Eigen::VectorXd>(buffer.data(), buffer.size());
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::serialize(std::vector<double> &) const {
  throw eckit::NotImplemented("ufo::ObsBiasIncrement::serialize not implemented", Here());
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::deserialize(const std::vector<double> &, std::size_t &) {
  throw eckit::NotImplemented("ufo::ObsBiasIncrement::deserialize not implemented", Here());
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::print(std::ostream & os) const {
  if (this->serialSize() > 0) {
    // map bias coeffs to eigen matrix
    Eigen::Map<const Eigen::MatrixXd>
      coeffs(biascoeffsinc_.data(), prednames_.size(), vars_.size());
    os << "ObsBiasIncrement::print " << std::endl;
    os << "---------------------------------------------------------------" << std::endl;
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

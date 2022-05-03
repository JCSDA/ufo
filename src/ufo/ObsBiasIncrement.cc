/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBiasIncrement.h"

#include <iomanip>
#include <memory>

#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "ufo/ObsBias.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ioda::ObsSpace & odb,
                                   const Parameters_ & params)
  : vars_(odb.assimvariables()) {
  oops::Log::trace() << "ufo::ObsBiasIncrement::create starting." << std::endl;

  // Predictor factory
  for (const PredictorParametersWrapper &wrapper :
       params.variationalBC.value().predictors.value()) {
    std::unique_ptr<PredictorBase> predictor(PredictorFactory::create(wrapper.predictorParameters,
                                                                      vars_));
    prednames_.push_back(predictor->name());
  }

  // initialize bias coefficient perturbations
  biascoeffsinc_ = Eigen::VectorXd::Zero(prednames_.size() * vars_.size());

  oops::Log::trace() << "ufo::ObsBiasIncrement::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBiasIncrement::ObsBiasIncrement(const ObsBiasIncrement & other, const bool copy)
  : prednames_(other.prednames_), vars_(other.vars_) {
  oops::Log::trace() << "ufo::ObsBiasIncrement::copy ctor starting" << std::endl;

  // Copy the bias coefficients data, or fill in with zeros
  if (copy) {
    biascoeffsinc_ = other.biascoeffsinc_;
  } else {
    biascoeffsinc_ = Eigen::VectorXd::Zero(prednames_.size() * vars_.size());
  }

  oops::Log::trace() << "ufo::ObsBiasIncrement::copy ctor done." << std::endl;
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

void ObsBiasIncrement::serialize(std::vector<double> & vect) const {
  std::vector<double> vec_obs(biascoeffsinc_.data(),
                              biascoeffsinc_.data() + biascoeffsinc_.size());
  vect.insert(vect.end(), vec_obs.begin(), vec_obs.end());
  oops::Log::trace() << "ObsBiasIncrement::serialize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::deserialize(const std::vector<double> & vect, std::size_t & index) {
  for (unsigned int jj = 0; jj < biascoeffsinc_.size(); ++jj) {
    biascoeffsinc_[jj] = vect[index];
    ++index;
  }
  oops::Log::trace() << "ObsBiasIncrement::deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasIncrement::print(std::ostream & os) const {
  if (this->serialSize() > 0) {
    // map bias coeffs to eigen matrix
    Eigen::Map<const Eigen::MatrixXd>
      coeffs(biascoeffsinc_.data(), prednames_.size(), vars_.size());
    os << "ufo::ObsBiasIncrement::print " << std::endl;
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
    os << "ufo::ObsBiasIncrement::print done" << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

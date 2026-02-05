/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/errors/ObsErrorDiagonal.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"


namespace ufo {

// -----------------------------------------------------------------------------

static ObsErrorMaker<ObsErrorDiagonal> makerDiagUFO_("diagonal");

// -----------------------------------------------------------------------------

ObsErrorDiagonal::ObsErrorDiagonal(const Parameters_ & params,
                                   ioda::ObsSpace & obsgeom,
                                   const eckit::mpi::Comm &timeComm)
  : ObsErrorBase(timeComm),
    stddev_(obsgeom, "ObsError"),
    inverseVariance_(obsgeom),
    params_(params)
{
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  oops::Log::trace() << "ObsErrorDiagonal:ObsErrorDiagonal constructed nobs = "
                     << stddev_.nobs() << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonal::update(const ioda::ObsVector & obserr) {
  stddev_ = obserr;
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  oops::Log::info() << "ObsErrorDiagonal covariance updated " << stddev_.nobs() << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonal::multiply(ioda::ObsVector & dy) const {
  dy /= inverseVariance_;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonal::inverseMultiply(ioda::ObsVector & dy) const {
  dy *= inverseVariance_;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonal::randomize(ioda::ObsVector & dy) const {
  if (params_.zeroMeanPerturbations.value()) {
    randomizeWithZeroEnsembleMean(dy);
  } else {
    randomizeWithoutZeroEnsembleMean(dy);
  }
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonal::save(const std::string & name) const {
  stddev_.save(name);
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorDiagonal::getObsErrors() const {
  return std::make_unique<ioda::ObsVector>(stddev_);
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorDiagonal::getInverseVariance() const {
  return std::make_unique<ioda::ObsVector>(inverseVariance_);
}

// -----------------------------------------------------------------------------
void ObsErrorDiagonal::print(std::ostream & os) const {
  os << "UFO Diagonal observation error covariance, inverse variances: "
     << inverseVariance_ << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonal::randomizeWithoutZeroEnsembleMean(ioda::ObsVector & dy) const {
  dy.random();
  dy *= this->stddev_;
  dy *= params_.pert.value();
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonal::randomizeWithZeroEnsembleMean(ioda::ObsVector & dy) const {
  ioda::ObsVector perturbation(dy);
  ioda::ObsVector sum(dy);
  sum.zero();

  // Generate initial independent perturbations for all ensemble members.
  // Calculate their sum and store this member's perturbations in 'dy'.
  for (int member = 1; member <= params_.numberOfMembers.value().value(); ++member) {
    perturbation.random();
    sum += perturbation;
    if (member == params_.member.value().value())
      dy = perturbation;
  }

  // Subtract the ensemble mean of perturbations from this member's perturbations.
  dy.axpy(-1.0 / params_.numberOfMembers.value().value(), sum);

  // Scale perturbations to the requested amplitude.
  dy *= stddev_;
  dy *= std::sqrt(params_.numberOfMembers.value().value() /
                 (params_.numberOfMembers.value().value() - 1.0)) * params_.pert.value();
}

// -----------------------------------------------------------------------------
}  // namespace ufo

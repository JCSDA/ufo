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

ObsErrorDiagonal::ObsErrorDiagonal(const Parameters_ & options, ioda::ObsSpace & obsgeom,
                                   const eckit::mpi::Comm &timeComm)
  : ObsErrorBase(timeComm),
    stddev_(obsgeom, "ObsError"), inverseVariance_(obsgeom), options_(options)
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
  dy.random();
  dy *= stddev_;
  dy *= options_.pert;
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


}  // namespace ufo

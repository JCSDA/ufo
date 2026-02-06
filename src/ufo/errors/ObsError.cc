/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/errors/ObsError.h"

namespace ufo {

ObsError::ObsError(const eckit::Configuration & obsErrConf, ioda::ObsSpace & odb)
                    : R_() {
  oops::Log::trace() << "ufo::ObsError::ObsError starting" << std::endl;
  Parameters_ params;
  if (!obsErrConf.empty()) {
    params.validateAndDeserialize(obsErrConf);
  } else {
    params = Parameters_();
  }
  R_.reset(ObsErrorFactory::create(params.errorParameters, odb));
  oops::Log::trace() << "ufo::ObsError::ObsError finished" << std::endl;
  }

void ObsError::multiply(ioda::ObsVector & dy) const {
  R_->multiply(dy);
}

void ObsError::inverseMultiply(ioda::ObsVector & dy) const {
  R_->inverseMultiply(dy);
}

void ObsError::randomize(ioda::ObsVector & dy) const {
  R_->randomize(dy);
}

void ObsError::save(const std::string & name) const {
  R_->save(name);
}

void ObsError::update(const ioda::ObsVector &dy) {
  R_->update(dy);
}

std::unique_ptr<ioda::ObsVector> ObsError::getObsErrors() const {
  return R_->getObsErrors();
}

std::unique_ptr<ioda::ObsVector> ObsError::getInverseVariance() const {
  return R_->getInverseVariance();
}

double ObsError::getRMSE() const {
  return R_->getRMSE();
}

void ObsError::localize(ioda::ObsVector & locvector) const {
  R_->localize(locvector);
}

int ObsError::localDim() const {
  return R_->localDim();
}

Eigen::MatrixXf ObsError::localInverseMultiply(const Eigen::MatrixXf &zz) const {
  return R_->localInverseMultiply(zz);
}

Eigen::VectorXd ObsError::local_invVarR() const {
  return R_->local_invVarR();
}

void ObsError::print(std::ostream & os) const {
  os << *R_;
}

}  // namespace ufo

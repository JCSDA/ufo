/*
 * (C) Copyright 2026 Bureau of Meteorology.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "ufo/errors/ObsErrorDiagonalInvGamma.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"


namespace ufo {

// -----------------------------------------------------------------------------

static ObsErrorMaker<ObsErrorDiagonalInvGamma> makerDiagUFO_("diagonal inverse gamma");

// -----------------------------------------------------------------------------

ObsErrorDiagonalInvGamma::ObsErrorDiagonalInvGamma(const Parameters_ & params,
                                   ioda::ObsSpace & obsgeom,
                                   const eckit::mpi::Comm &timeComm)
  : ObsErrorBase(timeComm),
    yobs_(obsgeom, "ObsValue"),
    stddev_(obsgeom, "ObsError"),
    inverseVariance_(obsgeom),
    params_(params)
{
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  oops::Log::trace() << "ObsErrorDiagonalInvGamma:ObsErrorDiagonalInvGamma constructed nobs = "
                     << stddev_.nobs() << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonalInvGamma::update(const ioda::ObsVector & obserr) {
  stddev_ = obserr;
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  oops::Log::trace() << "ObsErrorDiagonalInvGamma covariance updated "
                    << stddev_.nobs() << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonalInvGamma::multiply(ioda::ObsVector & dy) const {
  dy /= inverseVariance_;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonalInvGamma::inverseMultiply(ioda::ObsVector & dy) const {
  dy *= inverseVariance_;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonalInvGamma::localize(ioda::ObsVector & locvector) const {
  oops::Log::trace() << "ufo::ObsErrorDiagonalInvGamma::localize start" << std::endl;

  const double missing = util::missingValue<double>();

  assert(locvector.size() == stddev_.size());
  std::vector<double> localstdev;
  std::vector<double> localinvvar;
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    if (locvector[jj] != missing && locvector[jj] <= 0) {
      throw eckit::BadValue("Localization weights must be positive. Use "
                            "oops::util::missingValue<double>() to indicate "
                            "an observation with a weight of zero.");
    }
    if (locvector[jj] != missing && stddev_[jj] != missing) {
      localstdev.push_back(stddev_[jj] * std::pow(locvector[jj], -0.5));
      localinvvar.push_back(locvector[jj] * std::pow(stddev_[jj], -2.0));
    }
  }
  local_stddev_ = Eigen::Map<Eigen::VectorXd>(localstdev.data(), localstdev.size());
  local_inverseVariance_ =
    Eigen::Map<Eigen::VectorXd>(localinvvar.data(), localinvvar.size());
}

// -----------------------------------------------------------------------------

Eigen::MatrixXf ObsErrorDiagonalInvGamma::localInverseMultiply(const Eigen::MatrixXf & zz) const {
  Eigen::MatrixXf zzRinv(zz.rows(), zz.cols());
  for (int ii = 0; ii < zz.rows(); ++ii) {
    zzRinv(ii, Eigen::placeholders::all) = zz(ii, Eigen::placeholders::all)
                                 .cwiseProduct(local_inverseVariance_.cast<float>().transpose());
  }
  return zzRinv;
}

// -----------------------------------------------------------------------------

int ObsErrorDiagonalInvGamma::localDim() const {
  return local_inverseVariance_.size();
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonalInvGamma::randomize(ioda::ObsVector & dy) const {
  if (params_.oneMeanPerturbations.value()) {
    randomizeWithOneEnsembleMean(dy);
  } else {
    randomizeWithoutOneEnsembleMean(dy);
  }
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonalInvGamma::save(const std::string & name) const {
  stddev_.save(name);
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorDiagonalInvGamma::getObsErrors() const {
  return std::make_unique<ioda::ObsVector>(stddev_);
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorDiagonalInvGamma::getInverseVariance() const {
  return std::make_unique<ioda::ObsVector>(inverseVariance_);
}

// -----------------------------------------------------------------------------
void ObsErrorDiagonalInvGamma::print(std::ostream & os) const {
  os << "UFO Diagonal observation error covariance, inverse variances: "
     << inverseVariance_ << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonalInvGamma::randomizeWithoutOneEnsembleMean(ioda::ObsVector & dy) const {
  // Sample from IG with mean of 1 when number of samples is large
  dy.random("InverseGamma", params_.relvar.value());

  // Invert, since GIG filter requires yobs**2/ypert term
  dy.invert();

  // Scale by yobs since for IG, ypert = yobs*ypert_with_mean_of_1
  dy *= yobs_;

  // Subtract obs value to obtain deviation from obs
  dy -= yobs_;
}

// -----------------------------------------------------------------------------

void ObsErrorDiagonalInvGamma::randomizeWithOneEnsembleMean(ioda::ObsVector & dy) const {
  ioda::ObsVector perturbation(dy);
  ioda::ObsVector sum(dy);
  sum.zero();

  // Due to sampling error, ensemble mean of random inverse gamma values may not be
  // sufficiently close to 1; hence first scale raw perturbations to enforce
  // ensemble mean = 1 condition.
  // Generate initial independent perturbations for all ensemble members.
  // Calculate their sum and store this member's perturbations in 'dy'.
  if (params_.member.value() == boost::none) {
    throw eckit::UserError("Member number not specified");
  } else if (params_.numberOfMembers.value() == boost::none) {
      throw eckit::UserError("Number of members not specified");
  } else {
      for (int member = 1; member <= params_.numberOfMembers.value().value(); ++member) {
        perturbation.random("InverseGamma", params_.relvar.value());
        sum += perturbation;
        if (member == params_.member.value().value())
          dy = perturbation;
      }
  }
  dy /= sum;
  dy *= params_.numberOfMembers.value().value();

  // Invert, since GIG filter requires yobs**2/ypert term
  dy.invert();

  // Scale by obs value
  dy *= yobs_;

  // Subtract obs value to obtain deviation from obs
  dy -= yobs_;
}

// -----------------------------------------------------------------------------
}  // namespace ufo

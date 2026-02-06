/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERROR_H_
#define UFO_ERRORS_OBSERROR_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "oops/util/Printable.h"

#include "ufo/errors/ObsErrorBase.h"

namespace ufo {

// wrapper class for ufo ObsError implementations which appears in the
// ufo ObsTraits as the external facing ObsError because the abstract
// ObsErrorBase can't be used in the ObsTraits.
class ObsError : public util::Printable,
                 private boost::noncopyable {
 public:
  typedef ObsErrorParametersWrapper Parameters_;

  ObsError(const eckit::Configuration &, ioda::ObsSpace &);
  ObsError() = default;

  void multiply(ioda::ObsVector & dy) const;
  void inverseMultiply(ioda::ObsVector & dy) const;
  void randomize(ioda::ObsVector & dy) const;

  void save(const std::string &name) const;
  void update(const ioda::ObsVector &stddev);

  virtual std::unique_ptr<ioda::ObsVector> getObsErrors() const;
  virtual std::unique_ptr<ioda::ObsVector> getInverseVariance() const;
  virtual double getRMSE() const;

  void localize(ioda::ObsVector & locvector) const;
  int localDim() const;
  Eigen::MatrixXf localInverseMultiply(const Eigen::MatrixXf &zz) const;
  Eigen::VectorXd local_invVarR() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsErrorBase> R_;
};

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERROR_H_

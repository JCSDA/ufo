/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORBASE_H_
#define UFO_ERRORS_OBSERRORBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/PolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "ufo/errors/ObsErrorParametersBase.h"

namespace ioda {
  class ObsVector;
}

namespace ufo {

class ObsErrorBase : public util::Printable,
                     private boost::noncopyable {
 public:
  explicit ObsErrorBase(const eckit::mpi::Comm &timeComm) : timeComm_(timeComm) {}
  virtual ~ObsErrorBase() = default;

  virtual void multiply(ioda::ObsVector & dy) const = 0;
  virtual void inverseMultiply(ioda::ObsVector & dy) const = 0;

  /// Generate random perturbation in \p dy.
  virtual void randomize(ioda::ObsVector & dy) const = 0;

  /// Save obs errors to the group \p name.
  virtual void save(const std::string &name) const = 0;

  /// Set the diagonal of the covariance matrix to \p stddev squared.
  virtual void update(const ioda::ObsVector &stddev) = 0;

  /// Return a copy of obs error std. dev. If this ObsVector is modified (e.g. by obs filters),
  /// it should be passed back to update() to ensure the covariance matrix stays consistent.
  virtual std::unique_ptr<ioda::ObsVector> getObsErrors() const = 0;

  /// Return the vector of inverse obs error variances.
  virtual std::unique_ptr<ioda::ObsVector> getInverseVariance() const = 0;

  /// Get mean error for Jo table.
  virtual double getRMSE() const = 0;

  /// Create local R matrix.
  virtual void localize(ioda::ObsVector & locvector) const = 0;

  /// Return dimension of local R matrix.
  virtual int localDim() const = 0;

  /// Multiply a local obs vector \p zz by \f$R^{-1}\f$.
  virtual Eigen::MatrixXf localInverseMultiply(const Eigen::MatrixXf & zz) const = 0;

  /// Get local inverse variance
  virtual Eigen::VectorXd local_invVarR() const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
  const eckit::mpi::Comm & timeComm_;
};

// -----------------------------------------------------------------------------

class ObsErrorFactory;

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// ObsErrorParametersBase.
class ObsErrorParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of ObsErrorParametersBase
  /// controlling the behavior of an Obs Error operator. The type of the subclass is determined
  /// by the value of the "covariance model" key in the Configuration object from which this object
  /// is deserialized.
  oops::PolymorphicParameter<ObsErrorParametersBase, ObsErrorFactory>
    errorParameters{"covariance model", "diagonal", this};
};

// -----------------------------------------------------------------------------

class ObsErrorFactory{
 public:
  static ObsErrorBase * create(const ObsErrorParametersBase &params,
                               ioda::ObsSpace &odb);

  static std::unique_ptr<ObsErrorParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~ObsErrorFactory() = default;

 protected:
  explicit ObsErrorFactory(const std::string &name);

 private:
  // Private make called by the public create. Determines ctor signature
  // for classes derived from ObsErrorBase
  virtual ObsErrorBase * make(const ObsErrorParametersBase &,
                              ioda::ObsSpace &) = 0;

  virtual std::unique_ptr<ObsErrorParametersBase> makeParameters() const = 0;

  static std::map < std::string, ObsErrorFactory* > & getMakers() {
    static std::map < std::string, ObsErrorFactory* > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class ObsErrorMaker : public ObsErrorFactory {
  typedef typename T::Parameters_ Parameters_;

  ObsErrorBase * make(const ObsErrorParametersBase &params,
                      ioda::ObsSpace &odb) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(stronglyTypedParams, odb, odb.commTime());
  }
  std::unique_ptr<ObsErrorParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit ObsErrorMaker(const std::string & name) : ObsErrorFactory(name) {}
};

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORBASE_H_

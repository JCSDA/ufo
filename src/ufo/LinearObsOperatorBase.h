/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_LINEAROBSOPERATORBASE_H_
#define UFO_LINEAROBSOPERATORBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"
#include "ufo/ObsOperatorParametersBase.h"

#include "ufo/utils/VariableNameMap.h"

namespace oops {
class ObsVariables;
class Variables;
}

namespace ioda {
class ObsVector;
}

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;
class ObsBias;
class ObsBiasIncrement;
// -----------------------------------------------------------------------------
/// Base class for linear observation operators
///
/// The implementer should first define a subclass of ObsOperatorParametersBase
/// holding the settings of the operator in question. The LinearObsOperatorBase subclass should
/// then typedef `Parameters_` to the name of the ObsOperatorParametersBase subclass and provide a
/// constructor with the following signature:
///
///    SubclassName(const ioda::ObsSpace &, const Parameters_ &);
class LinearObsOperatorBase : public util::Printable,
                              private boost::noncopyable {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  explicit LinearObsOperatorBase(const ioda::ObsSpace & odb,
                                 const VariableNameMap & nameMap = VariableNameMap(boost::none))
           : odb_(odb), nameMap_(nameMap) {}
  virtual ~LinearObsOperatorBase() {}

/// Obs Operator
  virtual void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) = 0;
  virtual void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const = 0;
  virtual void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const = 0;

/// Operator input required from Model
  virtual const oops::Variables & requiredVars() const = 0;

/// \brief List of variables simulated by this operator.
///
/// The default implementation returns the list of all simulated variables in the ObsSpace.
  virtual oops::ObsVariables simulatedVars() const;

/// \brief The space containing the observations to be simulated by this operator.
  const ioda::ObsSpace &obsspace() const { return odb_; }

 private:
  virtual void print(std::ostream &) const = 0;

 private:
  const ioda::ObsSpace & odb_;
 protected:
  mutable VariableNameMap nameMap_;
};

// -----------------------------------------------------------------------------

class LinearObsOperatorFactory;

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// ObsOperatorParametersBase.
class LinearObsOperatorParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearObsOperatorParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of ObsOperatorParametersBase
  /// controlling the behavior of a linear observation operator. The type of the subclass is
  /// determined by the value of the "name" key in the Configuration object from which this object
  /// is deserialized.
  oops::RequiredPolymorphicParameter<ObsOperatorParametersBase, LinearObsOperatorFactory>
    operatorParameters{"name", this};
};

// -----------------------------------------------------------------------------

/// Linear obs operator factory
class LinearObsOperatorFactory {
 public:
  /// \brief Create and return a new linear observation operator.
  ///
  /// The type of the operator is determined by the `name` attribute of \p params. \p params must
  /// be an instance of the subclass of ObsOperatorParametersBase associated with that operator,
  /// otherwise an exception will be thrown.
  static LinearObsOperatorBase * create(const ioda::ObsSpace &, const ObsOperatorParametersBase &);

  /// \brief Create and return an instance of the subclass of ObsOperatorParametersBase
  /// storing parameters of linear observation operators of the specified type.
  static std::unique_ptr<ObsOperatorParametersBase> createParameters(
      const std::string &name);

  /// \brief Return the names of all operators that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~LinearObsOperatorFactory() = default;

 protected:
  /// \brief Register a maker able to create linear observation operators of type \p name.
  explicit LinearObsOperatorFactory(const std::string &name);

 private:
  virtual LinearObsOperatorBase * make(const ioda::ObsSpace &,
                                       const ObsOperatorParametersBase &) = 0;

  virtual std::unique_ptr<ObsOperatorParametersBase> makeParameters() const = 0;

  static std::map < std::string, LinearObsOperatorFactory * > & getMakers() {
    static std::map < std::string, LinearObsOperatorFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class LinearObsOperatorMaker : public LinearObsOperatorFactory {
  typedef typename T::Parameters_   Parameters_;

  LinearObsOperatorBase * make(const ioda::ObsSpace & odb,
                               const ObsOperatorParametersBase & params) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(odb, stronglyTypedParams);
  }

  std::unique_ptr<ObsOperatorParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit LinearObsOperatorMaker(const std::string & name) : LinearObsOperatorFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_LINEAROBSOPERATORBASE_H_

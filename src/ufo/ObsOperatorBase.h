/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSOPERATORBASE_H_
#define UFO_OBSOPERATORBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "ioda/ObsSpace.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"
#include "ufo/ObsOperatorParametersBase.h"

#include "ufo/utils/VariableNameMap.h"

namespace oops {
  template <typename OBS> class Locations;
  class Variables;
}

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;
  struct ObsTraits;

// -----------------------------------------------------------------------------
/// Base class for observation operators
///
/// The implementer should first define a subclass of ObsOperatorParametersBase
/// holding the settings of the operator in question. The ObsOperatorBase subclass should then
/// typedef `Parameters_` to the name of the ObsOperatorParametersBase subclass and provide a
/// constructor with the following signature:
///
///    SubclassName(const ioda::ObsSpace &, const Parameters_ &);
class ObsOperatorBase : public util::Printable,
                        private boost::noncopyable {
 public:
  typedef oops::Locations<ObsTraits> Locations_;

  explicit ObsOperatorBase(const ioda::ObsSpace & odb,
                           const VariableNameMap & nameMap = VariableNameMap(boost::none))
     : odb_(odb), nameMap_(nameMap) {}
  virtual ~ObsOperatorBase() {}

/// Obs Operator
  virtual void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const = 0;

/// Operator input required from Model
  virtual const oops::Variables & requiredVars() const = 0;

/// \brief Model variable interpolation paths.
///
/// The default implementation places a vertical interpolation path at the latitude, longitude and
/// time of each observation location and asks for all required model variables to be interpolated
/// along these paths.
  virtual Locations_ locations() const;

/// \brief List of variables simulated by this operator.
///
/// The default implementation returns the list of all simulated variables in the ObsSpace.
  virtual oops::Variables simulatedVars() const;

 private:
  virtual void print(std::ostream &) const = 0;
  const ioda::ObsSpace & odb_;
 protected:
  mutable VariableNameMap nameMap_;
};

// -----------------------------------------------------------------------------

class ObsOperatorFactory;

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// ObsOperatorParametersBase.
class ObsOperatorParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsOperatorParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of ObsOperatorParametersBase
  /// controlling the behavior of an observation operator. The type of the subclass is determined
  /// by the value of the "name" key in the Configuration object from which this object
  /// is deserialized.
  oops::RequiredPolymorphicParameter<ObsOperatorParametersBase, ObsOperatorFactory>
    operatorParameters{"name", this};
};

// -----------------------------------------------------------------------------

/// Obs Operator Factory
class ObsOperatorFactory {
 public:
  /// \brief Create and return a new observation operator.
  ///
  /// The type of the operator is determined by the `name` attribute of \p params. \p params must
  /// be an instance of the subclass of ObsOperatorParametersBase associated with that operator,
  /// otherwise an exception will be thrown.
  static ObsOperatorBase * create(const ioda::ObsSpace &, const ObsOperatorParametersBase &params);

  /// \brief Create and return an instance of the subclass of ObsOperatorParametersBase
  /// storing parameters of observation operators of the specified type.
  static std::unique_ptr<ObsOperatorParametersBase> createParameters(
      const std::string &name);

  /// \brief Return the names of all operators that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~ObsOperatorFactory() = default;

 protected:
  /// \brief Register a maker able to create observation operators of type \p name.
  explicit ObsOperatorFactory(const std::string &name);

 private:
  virtual ObsOperatorBase * make(const ioda::ObsSpace &, const ObsOperatorParametersBase &) = 0;

  virtual std::unique_ptr<ObsOperatorParametersBase> makeParameters() const = 0;

  static std::map < std::string, ObsOperatorFactory * > & getMakers() {
    static std::map < std::string, ObsOperatorFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class ObsOperatorMaker : public ObsOperatorFactory {
  typedef typename T::Parameters_   Parameters_;

  ObsOperatorBase * make(const ioda::ObsSpace & odb,
                         const ObsOperatorParametersBase & params) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(odb, stronglyTypedParams);
  }

  std::unique_ptr<ObsOperatorParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit ObsOperatorMaker(const std::string & name) : ObsOperatorFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSOPERATORBASE_H_

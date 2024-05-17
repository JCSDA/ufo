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

#include "ioda/ObsDataVector.h"
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
  class ObsVariables;
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
  typedef ioda::ObsDataVector<int> QCFlags_t;

  explicit ObsOperatorBase(const ioda::ObsSpace & odb,
                           const VariableNameMap & nameMap = VariableNameMap(boost::none))
     : odb_(odb), nameMap_(nameMap) {}
  virtual ~ObsOperatorBase() {}

// Obs Operator
  virtual void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                           const QCFlags_t &) const = 0;

/// \brief Required model variables.
  virtual const oops::Variables & requiredVars() const = 0;

/// \brief Return an object holding one or more collections of paths sampling the observation
/// locations and indicating along which of these sets of paths individual model variables should
/// be interpolated.
///
/// The default implementation places a vertical interpolation path at the latitude and longitude
/// of each observation location and asks for all required model variables to be interpolated along
/// these paths at the observation times.
  virtual Locations_ locations() const;

/// \brief List of variables simulated by this operator.
///
/// The default implementation returns the list of all simulated variables in the ObsSpace.
  virtual oops::ObsVariables simulatedVars() const;

/// \brief Convert values of model variables stored in the sampled format to the reduced format.
///
/// This typically consists in computing, for each location, a weighted average of all the sampled
/// profiles obtained by model field interpolation along the paths sampling that location.
///
/// The sampled and reduced representations of model variables sampled along exactly one path per
/// location are identical. The GeoVaLs class automatically detects such variables, stores
/// explicitly only their sampled representation and treats the reduced representation as an alias
/// for the sampled one. In consequence, this method needs to be reimplemented only in subclasses
/// whose locations() method requests some model variables to be sampled along multiple
/// interpolation paths per location and at least one of the following conditions is met:
///
/// 1. Any of these variables are needed by observation filters or bias predictors (which typically
///    operate on GeoVaLs stored in the reduced format).
///
/// 2. The operator wants to use the reduced representation of any of these variables in
///    simulateObs().
///
/// \param[in]
///   List of variables whose reduced representation should be computed.
/// \param[inout] geovals
///   On input, an object storing the values of model variables in the sampled format.
///   This function is expected to fill in the reduced representation of at least the variables
///   `vars` (it may optionally do so for other variables as well).
///
/// The default implementation instructs `geovals` to store the reduced representation of each
/// variable already available in the sampled format and checks if it is aliased (and therefore
/// identical) with the sampled representation. If that is not the case, the function throws an
/// exception.
///
/// The reason why by default this function operates on all variables rather than just those listed
/// in `vars` is that obs filters and bias predictors, which operate on GeoVaLs stored in the
/// reduced format, sometimes fail to request explicitly all the model variables they use.
/// Operators that request some variables to be sampled along multiple paths per location may want
/// to be less lenient and reduce only the variables in `vars`, since variable reduction may then
/// be costly or even ill-defined.
  virtual void computeReducedVars(const oops::Variables & vars, GeoVaLs & geovals) const;

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

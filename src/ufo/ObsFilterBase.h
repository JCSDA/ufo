/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/parameters/RequiredPolymorphicParameter.h"
#include "oops/util/Printable.h"

#include "ufo/ObsFilterParametersBase.h"

namespace ioda {
class ObsSpace;
class ObsVector;
template <typename DATA> class ObsDataVector;
}
namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

/// \brief Stage at which to run a filter.
///
/// AUTO indicates that the stage will be determined automatially based on the data requested
/// by the filter.
/// PRE, PRIOR and POST indicate that the filter will be run at the relevant stage
/// but a check will be performed to ensure that all of the data structures required by the filter
/// are present.
enum class FilterStage {AUTO,
                        PRE,
                        PRIOR,
                        POST};

/// \brief Base class for generic implementations of filters processing observations.
///
/// Use this class as a base class for generic implementations
/// and interface::ObsFilterBase as a base class for OBS-specific implementations.
///
/// Note: implementations of this interface can opt to extract their settings either from
/// a Configuration object or from a subclass of ObsFilterParametersBase.
///
/// In the former case, they should provide a constructor with the following signature:
///
///    ObsFilter(ObsSpace_ &, const eckit::Configuration &,
///              ObsDataPtr_<int>, ObsDataPtr_<float>);
///
/// In the latter case, the implementer should first define a subclass of ObsFilterParametersBase
/// holding the settings of the filter in question. The implementation of the ObsFilter interface
/// should then typedef `Parameters_` to the name of that subclass and provide a constructor with
/// the following signature:
///
///    ObsFilter(ObsSpace_ &, const Parameters_ &,
///              ObsDataPtr_<int>, ObsDataPtr_<float>);
// TODO(someone): check whether we can just keep Parameters here and remove the
//                Configuration option.
class ObsFilterBase : public util::Printable,
                      private boost::noncopyable {
 public:
  ObsFilterBase() {}
  virtual ~ObsFilterBase() {}

  /// \brief Perform any observation processing steps that do not require access to GeoVaLs or
  /// outputs produced by the observation operator.
  virtual void preProcess() = 0;

  /// \brief Perform any observation processing steps that require access to GeoVaLs, but not to
  /// outputs produced by the observation operator.
  virtual void priorFilter(const GeoVaLs &gv) = 0;

  /// \brief Perform any observation processing steps that require access to both GeoVaLs and
  /// outputs produced by the observation operator.
  ///
  /// \param gv
  ///   GeoVaLs.
  /// \param ov
  ///   Model equivalents produced by the observation operator.
  /// \param bv
  ///   Bias of departure produced by the observation operator.
  /// \param dv
  ///   Observation diagnostics produced by the observation operator.
  virtual void postFilter(const GeoVaLs & gv,
                          const ioda::ObsVector &ov,
                          const ioda::ObsVector &bv,
                          const ObsDiagnostics &dv) = 0;

  /// \brief Check the required filter data are present prior to running this filter.
  virtual void checkFilterData(const FilterStage filterStage) = 0;

  /// \brief Return the list of GeoVaLs required by this filter.
  virtual oops::Variables requiredVars() const = 0;

  /// \brief Return the list of observation diagnostics required by this filter.
  virtual oops::ObsVariables requiredHdiagnostics() const = 0;
};

// =============================================================================

class FilterFactory;

// -----------------------------------------------------------------------------

/// \brief Contains a polymorphic parameter holding an instance of a subclass of
/// ObsFilterParametersBase.
class ObsFilterParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsFilterParametersWrapper, Parameters)
 public:
  /// After deserialization, holds an instance of a subclass of ObsFilterParametersBase
  /// controlling the behavior of an observation filter. The type of the subclass is determined
  /// by the value of the "filter" key in the Configuration object from which this object
  /// is deserialized.
  oops::RequiredPolymorphicParameter<ObsFilterParametersBase, FilterFactory> filterParameters{
      "filter", this};

  /// Indices of iterations at which this filter should be applied.
  oops::OptionalParameter<std::string> applyAtIterations{"apply at iterations", this};
};

// =============================================================================

/// ObsFilter factory
class FilterFactory {
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ioda::ObsDataVector<DATA> >;

 public:
  /// \brief Create and return a new observation filter.
  ///
  /// The type of the filter is determined by the `filter` attribute of \p params. \p params
  /// must be an instance of the subclass of ObsFilterParametersBase associated with that filter,
  /// otherwise an exception will be thrown.
  static std::unique_ptr<ObsFilterBase> create(ioda::ObsSpace &,
                                               const ObsFilterParametersBase & params,
                                               ObsDataPtr_<int> flags = ObsDataPtr_<int>(),
                                               ObsDataPtr_<float> obserr = ObsDataPtr_<float>());

  /// \brief Create and return an instance of the subclass of ObsFilterParametersBase
  /// storing parameters of observation filters of the specified type.
  static std::unique_ptr<ObsFilterParametersBase> createParameters(
      const std::string &name);

  /// \brief Return the names of all filters that can be created by one of the
  /// registered makers.
  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~FilterFactory() = default;

 protected:
  /// \brief Register a maker able to create observation filters of type \p name.
  explicit FilterFactory(const std::string &name);

 private:
  virtual std::unique_ptr<ObsFilterBase> make(ioda::ObsSpace &,
                                              const ObsFilterParametersBase &,
                                              ObsDataPtr_<int>, ObsDataPtr_<float>) = 0;

  virtual std::unique_ptr<ObsFilterParametersBase> makeParameters() const = 0;

  static std::map < std::string, FilterFactory * > & getMakers() {
    static std::map < std::string, FilterFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class FilterMaker : public FilterFactory {
  /// Defined as T::Parameters_ if T defines a Parameters_ type; otherwise as
  /// icObsFilterParameters.
  typedef oops::TParameters_IfAvailableElseFallbackType_t<T, GenericObsFilterParameters>
          Parameters_;

  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ioda::ObsDataVector<DATA> >;

  std::unique_ptr<ObsFilterBase> make(ioda::ObsSpace & os,
                                      const ObsFilterParametersBase & params,
                                      ObsDataPtr_<int> flags,
                                      ObsDataPtr_<float> obserr) override {
        const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
        return std::make_unique<T>(os,
               oops::parametersOrConfiguration<oops::HasParameters_<T>::value>(stronglyTypedParams),
               flags,
               obserr);
  }

  std::unique_ptr<ObsFilterParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit FilterMaker(const std::string & name) : FilterFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

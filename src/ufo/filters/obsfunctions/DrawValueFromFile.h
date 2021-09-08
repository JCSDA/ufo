/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_DRAWVALUEFROMFILE_H_
#define UFO_FILTERS_OBSFUNCTIONS_DRAWVALUEFROMFILE_H_

#include <algorithm>     // transform
#include <list>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>       // pair
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/dataextractor/DataExtractor.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {


struct InterpMethodParameterTraitsHelper {
  typedef InterpMethod EnumType;
  static constexpr char enumTypeName[] = "InterpMethod";
  static constexpr util::NamedEnumerator<InterpMethod> namedValues[] = {
    { InterpMethod::EXACT, "exact" },
    { InterpMethod::NEAREST, "nearest" },
    { InterpMethod::LEAST_UPPER_BOUND, "least upper bound" },
    { InterpMethod::GREATEST_LOWER_BOUND, "greatest lower bound" },
    { InterpMethod::LINEAR, "linear" },
    { InterpMethod::BILINEAR, "bilinear" }
  };
};

struct ExtrapolationModeParameterTraitsHelper {
  typedef ExtrapolationMode EnumType;
  static constexpr char enumTypeName[] = "ExtrapolationMethod";
  static constexpr util::NamedEnumerator<ExtrapolationMode> namedValues[] = {
    { ExtrapolationMode::ERROR, "error" },
    { ExtrapolationMode::NEAREST, "nearest" },
    { ExtrapolationMode::MISSING, "missing" }
  };
};

}  // namespace ufo


namespace oops {

template <>
struct ParameterTraits<ufo::InterpMethod> :
    public EnumParameterTraits<ufo::InterpMethodParameterTraitsHelper>
{};

template <>
struct ParameterTraits<ufo::ExtrapolationMode> :
    public EnumParameterTraits<ufo::ExtrapolationModeParameterTraitsHelper>
{};

}  // namespace oops


namespace ufo {

/// \brief How to identify the relevant elements of the interpolated array along a dimension
/// indexed by a particular variable.
class InterpolationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(InterpolationParameters, Parameters)

 public:
  /// Name of the indexing variable (e.g. `latitude@MetaData`).
  oops::RequiredParameter<std::string> name{"name", this};

  /// Method used to map the value of a variable to a range of slices of the interpolated array
  /// along the dimension indexed by that variable.
  ///
  /// \see InterpMethod for the list of supported methods.
  oops::RequiredParameter<InterpMethod> method{"method", this};

  /// Extrapolation mode for the given variable and interpolation method i.e. behaviour for
  /// out-of-bounds extract.
  ///
  /// \sea ExtrapolationMode for a list of supported modes.
  oops::Parameter<ExtrapolationMode> extrapMode{
    "extrapolation mode", ExtrapolationMode::ERROR, this};
};


/// \brief Options controlling the DrawValueFromFile ObsFunction (excluding the `group` option).
class DrawValueFromFileParametersWithoutGroup : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DrawValueFromFileParametersWithoutGroup, Parameters)

 public:
  /// Path to the file containing the data to interpolate.
  oops::RequiredParameter<std::string> fpath{"file", this};

  /// List of interpolation variables and associated methods.
  /// Note that channel numbers is handled implicitly by the "channels" (see below).
  oops::RequiredParameter<std::vector<InterpolationParameters>> interpolation{"interpolation",
                                                                              this};
  /// List of channel numbers (then deriving an observation error per channel)
  /// If this option is provided, the channel number is implicitly prepended to the list of
  /// interpolation variables and matched exactly.
  oops::OptionalParameter<std::set<int>> chlist{"channels", this};
};


/// \brief Options controlling the DrawValueFromFile ObsFunction
class DrawValueFromFileParameters : public DrawValueFromFileParametersWithoutGroup {
  OOPS_CONCRETE_PARAMETERS(DrawValueFromFileParameters,
                           DrawValueFromFileParametersWithoutGroup)

 public:
  /// The file should contain exactly one variable from this group. This is the variable that
  /// will be interpolated.
  oops::RequiredParameter<std::string> group{"group", this};
};


/// \brief Produce values by interpolating an array loaded from a file, indexed by
/// coordinates whose names correspond to ObsSpace variables.
///
/// \tparam T
///   Type of values produced by this ObsFunction. Must be `float`, `int` or `std::string`.
///
/// \details See DataExtractor for details on the format of this file.
///
/// ### example configurations: ###
///
/// \code{.yaml}
///     - filter: Variable Assignment
///       assignments:
///       - name: interpolated_value@DerivedObsValue
///         function:
///           name: DrawValueFromFile@ObsFunction
///           channels: 1-3
///           options:
///             file: <filepath>
///             channels: 1-3
///             group: DerivedObsValue
///             interpolation:
///             - name: satellite_id@MetaData
///               method: exact
///             - name: processing_center@MetaData
///               method: exact
///             - name: air_pressure@MetaData
///               method: linear
/// \endcode
///
/// Note that channel number extraction is implicit, using the channels selected and performed as
/// an exact match before any user defined interpolation takes place.
template <typename T>
class DrawValueFromFile : public ObsFunctionBase<T> {
 public:
  explicit DrawValueFromFile(const eckit::LocalConfiguration &);

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<T> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  Variables allvars_;
  std::unordered_map<std::string, InterpMethod> interpMethod_;
  std::unordered_map<std::string, ExtrapolationMode> extrapMode_;
  std::string fpath_;
  DrawValueFromFileParameters options_;
  std::vector<int> channels_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_DRAWVALUEFROMFILE_H_

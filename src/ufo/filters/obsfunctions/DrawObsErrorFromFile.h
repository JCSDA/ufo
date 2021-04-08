/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_
#define UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_

#include <algorithm>     // transform
#include <list>
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
#include "ufo/utils/ObsErrorNetcdfHandler.h"
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
    { InterpMethod::LINEAR, "linear" }
  };
};

}  // namespace ufo


namespace oops {

template <>
struct ParameterTraits<ufo::InterpMethod> :
    public EnumParameterTraits<ufo::InterpMethodParameterTraitsHelper>
{};

}  // namespace oops


namespace ufo {

/// \brief Options controlling the interpolation method
class InterpolationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(InterpolationParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> name{"name", this};
  oops::RequiredParameter<InterpMethod> method{"method", this};
};


/// \brief Options controlling DrawObsErrorFromFile ObsFunction
class DrawObsErrorFromFileParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DrawObsErrorFromFileParameters, Parameters)

 public:
  /// Filepath to the observation variance.
  oops::RequiredParameter<std::string> fpath{"file", this};
  /// List of interpolation variables and associated methods (exact, linear and nearest supported)
  oops::RequiredParameter<std::vector<InterpolationParameters>> interpolation{"interpolation",
                                                                              this};
};


/// \brief Derive observation error from a file representing the variance or the covariance
/// matrix.
/// \details See ObsErrorNetcdfErrorHandler for details on the format of this file.
///
/// ### example configurations: ###
///
/// \code{.yaml}
///     - Filter: Perform Action
///       filter variables:
///       - name: air_temperature
///       action:
///         name: assign error
///         error function:
///           name: DrawObsErrorFromFile@ObsFunction
///           options:
///             file: <filepath>
///             interpolation:
///             - name: channel_number@MetaData
///               method: exact
///             - name: satellite_id@MetaData
///               method: exact
///             - name: processing_center@MetaData
///               method: exact
///             - name: air_pressure@MetaData
///               method: linear
/// \endcode
///
class DrawObsErrorFromFile : public ObsFunctionBase {
 public:
  typedef std::list<std::pair<std::string, boost::variant<std::vector<int>,
                                                          std::vector<float>,
                                                          std::vector<std::string>>
                             >> ObData;
  static const std::string classname() {return "DrawObsErrorFromFile";}

  explicit DrawObsErrorFromFile(const eckit::LocalConfiguration);
  ~DrawObsErrorFromFile();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  Variables allvars_;
  std::unordered_map<std::string, InterpMethod> interpMethod_;
  std::string fpath_;
  DrawObsErrorFromFileParameters options_;

  static std::string get_full_name(const ufo::Variable &variable) {
    std::string name = variable.variable();
    if (variable.group().size() > 0) name += ("@"+variable.group());
    return name;
  }

  /// \brief This is a convenience function for updating our container for useful observation data
  template <typename T>
  static void updateObData(const ObsFilterData &in, const ufo::Variable &var, ObData &obData) {
    std::vector<T> dat;
    in.get(var, dat);
    std::string name = get_full_name(var);
    obData.emplace_back(name, std::move(dat));
  }

  /// \brief Add datetime observation information data to our container.
  /// \details We simply convert the datetimes to strings as our implementation is not discriminate
  /// between the two types.
  static void updateObDataDateTime(const ObsFilterData &in, const ufo::Variable &var,
                                   ObData &obData) {
    std::vector<util::DateTime> dat;
    std::vector<std::string> datConv;
    in.get(var, dat);
    datConv.resize(dat.size());

    // Convert the vec. of datetime. to strings
    std::transform(dat.begin(), dat.end(), datConv.begin(),
                   [](util::DateTime dt){return dt.toString();});
    std::string name = get_full_name(var);
    obData.emplace_back(name, std::move(datConv));
  }
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_DRAWOBSERRORFROMFILE_H_

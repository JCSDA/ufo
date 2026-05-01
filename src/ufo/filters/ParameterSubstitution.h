/*
 * (C) Crown Copyright 2024 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_PARAMETERSUBSTITUTION_H_
#define UFO_FILTERS_PARAMETERSUBSTITUTION_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/ObsFilterBase.h"
#include "ufo/ObsFilterParametersBase.h"

namespace ufo {

/// Parameters controlling the ParameterSubstitution filter.
class ParameterSubstitutionParameters : public ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ParameterSubstitutionParameters, ObsFilterParametersBase)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> sectionToRepeat{"section to repeat",
      "Filter configuration that will be repeated with variations",
      this};
  oops::RequiredParameter<std::vector<eckit::LocalConfiguration>> repetitions{"repetitions",
      "Configuration options to change for each repetition of the filter",
      this};
};

/// \brief Enables a filter to be run multiple times with chosen parameters varied each time.
///
/// The base filter configuration is specified in the `section to repeat` parameter.
/// All parameters that the user wishes to substitute must be present in the filter's configuration.
/// These parameters can be assigned arbitrary values in this section, but it is recommended
/// to use `{}` to signify a parameter that will be substituted.
/// For example, to indicate that the parameter `min_horizontal_spacing` will be repeated, the
/// user can use the following line:
///
///    min_horizontal_spacing: {}
///
/// The `repetitions` parameter contains a list of parameters to repeat.
/// Each of these parameters is assigned a list of the values that will be used when repeating
/// the filter. On the first iteration of the filter, the first value of each parameter is
/// substituted into the relevant part of the base filter configuration, and the filter is run.
/// The same then occurs for any subsequent repetitions.
/// Changes made in an earlier repetition do not carry over to subsequent repetitions.
///
/// For example, to indicate that a filter should be run twice, first with its `shuffle`
/// parameter set to `true`, and then to `false`, the following syntax should be used:
///
///    repetitions:
///    - shuffle:
///      - true
///      - false
///
/// It is possible to repeat complex parameters such as a `where` block.
/// To ensure a compact and readable yaml, it is recommended to use the 'yaml flow' syntax, e.g.
///
///    repetitions:
///    - where:
///      - [{variable: {name: MetaData/superObservation}, is_in: 0}]
///      - [{variable: {name: MetaData/superObservation}, is_in: 1}]
///
/// Note the need to use `[{` and `}]` around the values of interest.
///
/// Exceptions are thrown in the following circumstances:
/// * Attempting to replace an invalid filter parameter.
/// * Specifying an inconsistent number of repetitions for different parameters.
///
/// An example yaml is as follows:
///
///    - filter: Parameter Substitution
///      section to repeat:
///        filter: Poisson Disk Thinning
///        min_horizontal_spacing: {}
///        exclusion_volume_shape: ellipsoid
///        shuffle: {}
///      repetitions:
///      - min_horizontal_spacing:
///        - { "0": 2000, "1": 1000 }
///        - { "0": 3000, "1": 2000 }
///        - { "0": 4000, "1": 3000 }
///      - shuffle
///        - true
///        - false
///        - true
///
/// In this case, the `Poisson Disk Thinning` filter is run three times; each time
/// different values of `min_horizontal_spacing` and `shuffle` are used.
/// The `exclusion volume shape` parameter remains the same each repetition.
class ParameterSubstitution : public ObsFilterBase,
  private util::ObjectCounter<ParameterSubstitution> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ParameterSubstitutionParameters Parameters_;

  static const std::string classname() {return "ufo::ParameterSubstitution";}

  ParameterSubstitution(ioda::ObsSpace &,
                        const Parameters_ &,
                        ioda::ObsDataVector<int> &,
                        ioda::ObsDataVector<float> &);

  ~ParameterSubstitution() {}

  void preProcess() override;
  void priorFilter(const GeoVaLs &) override;
  void postFilter(const GeoVaLs &,
                  const ioda::ObsVector &,
                  const ioda::ObsVector &,
                  const ObsDiagnostics &) override;
  void checkFilterData(const FilterStage) override;
  oops::Variables requiredVars() const override;
  oops::ObsVariables requiredHdiagnostics() const override;
  void print(std::ostream & os) const override;

 private:
  const Parameters_ parameters_;
  std::vector<std::unique_ptr<ObsFilterBase>> filters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PARAMETERSUBSTITUTION_H_

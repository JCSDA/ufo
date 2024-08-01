/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_CREATEDIAGNOSTICFLAGS_H_
#define UFO_FILTERS_CREATEDIAGNOSTICFLAGS_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/ObsProcessorBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief Name and initialization options of a flag to be set up by the CreateDiagnosticFlags
/// filter.
class DiagnosticFlagParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DiagnosticFlagParameters, Parameters)

 public:
  /// Flag name.
  oops::RequiredParameter<std::string> name{"name", this};

  /// Initial value for the flag.
  oops::Parameter<bool> initialValue{"initial value", false, this};

  /// Determines what happens if the flag already exists.
  ///
  /// By default, the flag is not reinitialized, i.e. its current value is preserved.
  /// Set `force reinitialization` to `true` to reset the flag to `initial value`.
  oops::Parameter<bool> forceReinitialization{"force reinitialization", false, this};
};

/// \brief Options controlling the operation of the CreateDiagnosticFlags filter.
///
/// The filter creates a Boolean variable `DiagnosticFlags/<flag>/<var>` for each flag `<flag>`
/// specified in the `flags` list and each filter variable `<var>`. The initial values of these
/// variables can be customized using the `initial value` and `force reinitialization` options.
class CreateDiagnosticFlagsParameters : public oops::ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(CreateDiagnosticFlagsParameters, ObsFilterParametersBase)

 public:
  /// Flag to use bitmap diagnostic flags (if it's set, "flags" option is ignored)
  oops::Parameter<bool> bitMap{"bitmap diagnostic flags", false, this};

  /// Determines what happens if the bitmap flag already exists.
  ///
  /// By default, the flag is not reinitialized, i.e. its current value is preserved.
  /// Set `force reinitialization` to `true` to reset the flag to `initial value`.
  oops::Parameter<bool> forceReinitialization{"force bitmap reinitialization", false, this};

  /// The list of flags to create.
  oops::Parameter<std::vector<DiagnosticFlagParameters>> flags{"flags", {}, this};

  /// Simulated variables (and channels) for which the flags will be created.
  /// If not specified, defaults to all simulated variables in the ObsSpace.
  ///
  /// For example, if the filter is asked to create the flag `Unfolded` and the `filter variables`
  /// option is set to `[eastward_wind, northward_wind]` while the ObsSpace contains simulated
  /// variables `air_temperature`, `eastward_wind` and `northward_wind`, the filter will create
  /// Boolean variables `DiagnosticFlags/Unfolded/eastward_wind` and
  /// `DiagnosticFlags/Unfolded/northward_wind`, but not
  /// `DiagnosticFlags/Unfolded/air_temperature`.
  oops::OptionalParameter<std::vector<Variable>> filterVariables{
    "filter variables", this};

  /// If set to true, creates diagnostic flag with variable name
  /// `DiagnosticFlags/<flag>/observationReport` for each flag `<flag>` specified in the
  /// `flags` list.
  oops::Parameter<bool> observationReportFlags{"create observation report flags", false, this};

  /// If set to true, the filter will be executed only after the obs operator has been invoked.
  oops::Parameter<bool> deferToPost{"defer to post", false, this};
};

/// \brief Create and initialize new diagnostic flags.
///
/// \see CreateDiagnosticFlagsParameters for the list of accepted YAML options.
class CreateDiagnosticFlags : public ObsProcessorBase,
                              private util::ObjectCounter<CreateDiagnosticFlags> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef CreateDiagnosticFlagsParameters Parameters_;

  static const std::string classname() {return "ufo::CreateDiagnosticFlags";}

  CreateDiagnosticFlags(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                        std::shared_ptr<ioda::ObsDataVector<int>> qcflags,
                        std::shared_ptr<ioda::ObsDataVector<float>> obserr);
  ~CreateDiagnosticFlags() override;

 private:
  /// Creates a variable DiagnosticsFlags/
  void doFilter() override;

  void print(std::ostream &) const override;

  oops::ObsVariables getFilterVariables() const;

  template<class T>
  void createFlag(const std::string & flagName, const std::string & varName,
                  bool forceReinitialization, T initialValue) const;

 private:
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_CREATEDIAGNOSTICFLAGS_H_

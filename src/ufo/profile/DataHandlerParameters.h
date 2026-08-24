/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_DATAHANDLERPARAMETERS_H_
#define UFO_PROFILE_DATAHANDLERPARAMETERS_H_

#include <map>
#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/profile/ModelParameters.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

  /// \brief Options controlling the operation of the EntireSampleDataHandler
  /// and ProfileDataHandler classes.
  class DataHandlerParameters : public oops::Parameters {
     OOPS_CONCRETE_PARAMETERS(DataHandlerParameters, Parameters)

   public:  // functions
    /// Determine whether a variable group is optional or not.
    bool getOptional(const std::string &groupname) const;

    /// Determine number of entries per profile for a variable group.
    size_t getEntriesPerProfile(const std::string &groupname) const;

   public:  // variables
    /// Groups of variables whose presence in the input sample is optional
    /// (if not present, all variables are initially set to zero).
    oops::Parameter<std::vector<std::string>> groups_optional
      {"groups_optional", {"Corrections", "Counters"}, this};

    /// Groups of variables which have one value per profile.
    /// Needed for some unit tests.
    oops::Parameter<std::vector<std::string>> groups_singlevalue
      {"groups_singlevalue", {}, this};

    /// Groups of variables which are on model levels.
    /// Needed for some unit tests.
    oops::Parameter<std::vector<std::string>> groups_modellevels
      {"groups_modellevels",
          {"ModelLevelsDerivedValue",
              "ModelLevelsDiagnosticFlags/FinalQCRejection",
              "ModelLevelsDiagnosticFlags/PermanentStationRejection",
              "ModelLevelsDiagnosticFlags/PartialLayerUsedInAveraging"},
          this};

    /// Groups of variables which are on model rho levels.
    /// Needed for some unit tests.
    oops::Parameter<std::vector<std::string>> groups_modelrholevels
      {"groups_modelrholevels",
          {"ModelRhoLevelsDerivedValue", "ModelRhoLevelsFlags"}, this};

    /// Number of errors, accumulated over checks, that cause the observation to have failed.
    oops::Parameter<int> nErrorsFail {"nErrorsFail", 1, this};

    /// Maximum number of profile levels to be processed (a legacy of the OPS code).
    /// No maximum is assigned if this parameter is not specified.
    oops::OptionalParameter<int> maxlev {"maxlev", this};

    /// If not sorting observations, ensure number of profiles is consistent
    oops::Parameter<bool> ValidateTotalNumProf {"ValidateTotalNumProf", true, this};

    /// Default vertical coordinate to use in the slant path location algorithm.
    /// This can be overridden for each variable by using the \p alternativeVerticalCoordinate
    /// option.
    oops::Parameter<std::string> defaultVerticalCoordinate
      {"defaultVerticalCoordinate", "air_pressure", this};

    /// Alternative vertical coordinate(s) to use in the slant path location algorithm.
    /// The first string in each pair is the name of the variable whose slanted profile is to
    /// be determined, and the second string is the vertical coordinate that should be used to
    /// find the slant path locations for the variable.
    /// This will typically be useful for models whose variables appear on staggered vertical levels
    /// (e.g. with vertical coordinates 'air_pressure' and 'air_pressure_levels').
    oops::Parameter<std::map<std::string, std::string>> alternativeVerticalCoordinate
      {"alternativeVerticalCoordinate",
          {{"eastward_wind", "air_pressure_levels"}, {"northward_wind", "air_pressure_levels"},
                {"ExnerPA", "air_pressure_levels"}, {"LogPA", "air_pressure_levels"}}, this};

    /// Parameters related to the model.
    ModelParameters ModParameters{this};
  };
}  // namespace ufo

#endif  // UFO_PROFILE_DATAHANDLERPARAMETERS_H_


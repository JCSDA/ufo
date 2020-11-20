/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_DATAHANDLERPARAMETERS_H_
#define UFO_PROFILE_DATAHANDLERPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

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
    bool getOptional(const std::string &groupname) const
    {
      bool optional = false;
      if (std::find(groups_optional.value().begin(), groups_optional.value().end(), groupname)
          != groups_optional.value().end())
        optional = true;
      return optional;
    }

    /// Determine number of entries per profile for a variable group.
    size_t getEntriesPerProfile(const std::string &groupname) const
    {
      size_t entriesPerProfile = 0;
      // Variables with one entry per profile.
      if (std::find(groups_singlevalue.value().begin(), groups_singlevalue.value().end(), groupname)
          != groups_singlevalue.value().end()) {
        entriesPerProfile = 1;
      }
      return entriesPerProfile;
    }

   public:  // variables
    /// Groups of variables whose presence in the input sample is optional
    /// (if not present, all variables are initially set to zero).
    oops::Parameter<std::vector<std::string>> groups_optional
      {"groups_optional", {"Corrections", "Counters"}, this};

    /// Groups of variables which have one value per profile.
    oops::Parameter<std::vector<std::string>> groups_singlevalue
      {"groups_singlevalue", {"Counters"}, this};

    /// Number of errors, accumulated over checks, that cause the observation to have failed.
    oops::Parameter<int> nErrorsFail {"nErrorsFail", 1, this};

    /// Maximum number of profile levels to be processed (a legacy of the OPS code).
    /// No maximum is assigned if this parameter is not specified.
    oops::OptionalParameter<int> maxlev {"maxlev", this};

    /// If not sorting observations, ensure number of profiles is consistent
    oops::Parameter<bool> ValidateTotalNumProf {"ValidateTotalNumProf", true, this};
  };
}  // namespace ufo

#endif  // UFO_PROFILE_DATAHANDLERPARAMETERS_H_


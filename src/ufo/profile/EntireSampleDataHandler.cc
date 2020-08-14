/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "ufo/profile/EntireSampleDataHandler.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {
  EntireSampleDataHandler::EntireSampleDataHandler(ioda::ObsSpace &obsdb,
                                                   const ProfileConsistencyCheckParameters &options)
    : obsdb_(obsdb),
      options_(options)
  {}

  void EntireSampleDataHandler::writeQuantitiesToObsdb()
  {
    // Write out all variables in the QCFlags and Corrections groups.
    for (const auto& it_data : entireSampleData_) {
      std::string fullname = it_data.first;
      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(fullname, varname, groupname);

      if (groupname == "QCFlags") {
        putDataVector(fullname, get<int>(fullname));
      } else if (groupname == "Corrections") {
        putDataVector(fullname, get<float>(fullname));
      }
    }

    // Write out the NumAnyErrors counter, which is used in multiple QC checks.
    const std::vector <int> &NumAnyErrors =
      get<int>(ufo::VariableNames::name_counter_NumAnyErrors);
    putDataVector(ufo::VariableNames::name_counter_NumAnyErrors, NumAnyErrors);
  }
}  // namespace ufo

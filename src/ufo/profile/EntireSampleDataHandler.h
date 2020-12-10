/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_ENTIRESAMPLEDATAHANDLER_H_
#define UFO_PROFILE_ENTIRESAMPLEDATAHANDLER_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "boost/variant.hpp"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/missingValues.h"

#include "ufo/profile/DataHandlerParameters.h"

#include "ufo/utils/metoffice/MetOfficeQCFlags.h"
#include "ufo/utils/StringUtils.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  /// \brief Retrieve and store data for entire sample.
  /// This class uses lazy loading; vectors of variables are retrieved once requested
  /// and cached after that.
  /// Variables in certain groups are optional, meaning that if they are not present on
  /// the obsdb they will be filled with a default value if requested.
  class EntireSampleDataHandler {
   public:
    EntireSampleDataHandler(ioda::ObsSpace &obsdb,
                            const DataHandlerParameters &options);

    /// Retrieve a vector containing the requested variable for the entire data sample.
    ///    -# If the variable has previously been placed in a vector, return the vector.
    ///    -# If the variable is present in the input data set, fill the vector with those values.
    ///    -# If the variable is not present in the input data set, and 'optional' is true,
    ///       fill the vector with zeros.
    ///    -# If the variable is not present in the input data set, and 'optional' is false,
    ///       do not fill the vector.
    /// Also store the name of the variable, enabling it to be retrieved later.
    template <typename T>
      std::vector<T>& get(const std::string &fullname)
      {
        // Determine variable and group names, optional, and number of entries per profile.
        std::string varname;
        std::string groupname;
        ufo::splitVarGroup(fullname, varname, groupname);
        bool optional = options_.getOptional(groupname);
        size_t entriesPerProfile = options_.getEntriesPerProfile(groupname);

        std::vector <T> vec_all;  // Vector storing data for entire sample.
        if (entireSampleData_.find(fullname) != entireSampleData_.end()) {
          // If the vector is already present, return it.
          // If the type T is incorrect then boost::get will return an exception.
          // Provide additional information if that occurs.
          try {
            return boost::get<std::vector<T>> (entireSampleData_[fullname]);
          } catch (boost::bad_get) {
            throw eckit::BadParameter("Template parameter passed to boost::get "
                                      "probably has the wrong type", Here());
          }
        } else if (obsdb_.has(groupname, varname) || optional) {
          // Initially fill the vector with the default value for the type T.
          if (entriesPerProfile == 0) {
            vec_all.assign(obsdb_.nlocs(), defaultValue(vec_all));
          } else {
            vec_all.assign(entriesPerProfile * obsdb_.nrecs(), defaultValue(vec_all));
          }
          // Retrieve variable from the obsdb if present, overwriting the default value.
          if (obsdb_.has(groupname, varname)) obsdb_.get_db(groupname, varname, vec_all);
        }

        // If the vector contains entirely missing values, clear it.
        const T missingValue = util::missingValue(missingValue);  // Missing value for type T.
        bool allMissing = true;  // Signifies all elements in the vector are missing.
        for (size_t idx = 0; allMissing && idx < vec_all.size(); ++idx)
          allMissing = vec_all[idx] == missingValue;
        if (allMissing) {
          oops::Log::debug() << "All elements of " << fullname << " are missing" << std::endl;
          vec_all.clear();
        }

        // Add vector to map (even if it is empty).
        entireSampleData_.emplace(fullname, std::move(vec_all));
        return boost::get<std::vector<T>> (entireSampleData_[fullname]);
      }

    /// Write various quantities to the obsdb so they can be used in future QC checks.
    /// The particular variables written out are hardcoded but this could be changed to a
    /// configurable list if requred.
    void writeQuantitiesToObsdb();

   private:
    /// Put entire data vector on obsdb.
    template <typename T>
      void putDataVector(const std::string &fullname,
                         const std::vector <T> &datavec)
      {
        // Do not store the vector if it is empty.
        if (datavec.empty()) return;

        std::string varname;
        std::string groupname;
        ufo::splitVarGroup(fullname, varname, groupname);
        obsdb_.put_db(groupname, varname, datavec);
      }

    /// Observation database.
    ioda::ObsSpace &obsdb_;

    /// Configurable parameters.
    const DataHandlerParameters &options_;

    /// Default value used to fill vector of integers.
    int defaultValue(const std::vector <int> &vec) {return 0;}

    /// Default value used to fill vector of floats.
    float defaultValue(const std::vector <float> &vec) {return 0.0f;}

    /// Default value used to fill vector of strings.
    std::string defaultValue(const std::vector <std::string> &vec) {return "";}

    /// Container of each variable in the entire data set.
    std::unordered_map <std::string, boost::variant
      <std::vector <int>, std::vector <float>, std::vector <std::string>>> entireSampleData_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_ENTIRESAMPLEDATAHANDLER_H_

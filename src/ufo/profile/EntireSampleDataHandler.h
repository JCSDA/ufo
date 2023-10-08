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

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

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
    EntireSampleDataHandler(const ObsFilterData &data,
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
        const bool optional = options_.getOptional(groupname);
        const size_t entriesPerProfile = options_.getEntriesPerProfile(groupname);

        std::vector <T> vec_all;  // Vector storing data for entire sample.
        auto it_entireSampleData = entireSampleData_.find(fullname);
        if (it_entireSampleData != entireSampleData_.end()) {
          // If the vector is already present, return it.
          // If the type T is incorrect then boost::get will return an exception.
          // Provide additional information if that occurs.
          try {
            return boost::get<std::vector<T>> (it_entireSampleData->second);
          } catch (boost::bad_get) {
            throw eckit::BadParameter("Template parameter passed to boost::get for " +
                                      fullname + " probably has the wrong type", Here());
          }
        } else if (data_.has(Variable(fullname)) || optional) {
          // Initially fill the vector with the default value for the type T.
          if (entriesPerProfile == 0) {
            vec_all.assign(obsdb_.nlocs(), defaultValue(vec_all, groupname));
          } else {
            vec_all.assign(entriesPerProfile * obsdb_.nrecs(), defaultValue(vec_all, groupname));
          }
          // Retrieve variable from the obsdb if present, overwriting the default value.
          if (data_.has(Variable(fullname))) {
            if (entriesPerProfile == 0) {
              data_.get(Variable(fullname), vec_all);
            } else {
              obsdb_.get_db(groupname, varname, vec_all);
            }
          }
        }

        // Add vector to map.
        entireSampleData_.emplace(fullname, std::move(vec_all));
        return boost::get<std::vector<T>> (entireSampleData_[fullname]);
      }

    /// Write various quantities to the obsdb so they can be used in future QC checks.
    /// The particular variables written out are hardcoded but this could be changed to a
    /// configurable list if requred.
    void writeQuantitiesToObsdb();

    /// Initialise vector in the entire sample for a variable that is not currently
    /// stored. Fill the vector with the default value for the data type.
    template <typename T>
      void initialiseVector(const std::string fullname)
      {
        auto it_entireSampleData = entireSampleData_.find(fullname);
        if (it_entireSampleData == entireSampleData_.end() ||
            (it_entireSampleData != entireSampleData_.end() &&
             get<T>(fullname).size() == 0)) {
          std::string varname;
          std::string groupname;
          ufo::splitVarGroup(fullname, varname, groupname);
          const size_t entriesPerProfile = options_.getEntriesPerProfile(groupname);
          std::vector <T> vec_all;  // Vector storing data for entire sample.
          if (entriesPerProfile == 0) {
            vec_all.assign(obsdb_.nlocs(), defaultValue(vec_all, groupname));
          } else {
            vec_all.assign(entriesPerProfile * obsdb_.nrecs(),
                           defaultValue(vec_all, groupname));
          }
          entireSampleData_[fullname] = vec_all;
        }
      }

   private:  // functions
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

   private:  // variables
    /// ObsFilterData
    const ObsFilterData &data_;

    /// Observation database.
    ioda::ObsSpace &obsdb_;

    /// Configurable parameters.
    const DataHandlerParameters &options_;

    /// Default value used to fill vector of integers.
    int defaultValue(const std::vector <int> &vec, const std::string &groupname);

    /// Default value used to fill vector of floats.
    float defaultValue(const std::vector <float> &vec, const std::string &groupname);

    /// Default value used to fill vector of strings.
    std::string defaultValue(const std::vector <std::string> &vec, const std::string &groupname);

    /// Default value used to fill vector of booleans.
    bool defaultValue(const std::vector <bool> &vec, const std::string &groupname);

    /// Container of each variable in the entire data set.
    std::unordered_map <std::string, boost::variant
                        <std::vector <int>,
                         std::vector <float>,
                         std::vector <std::string>,
                         std::vector <bool>>> entireSampleData_;

    /// Missing value (int)
    const int missingValueInt = util::missingValue<int>();

    /// Missing value (float)
    const float missingValueFloat = util::missingValue<float>();

    /// Missing value (string)
    const std::string missingValueString = util::missingValue<std::string>();
  };
}  // namespace ufo

#endif  // UFO_PROFILE_ENTIRESAMPLEDATAHANDLER_H_

/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKBASE_H_
#define UFO_PROFILE_PROFILECHECKBASE_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/types/FloatCompare.h"

#include "oops/base/Variables.h"
#include "oops/util/CompareNVectors.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/PropertiesOfNVectors.h"

#include "ufo/filters/ConventionalProfileProcessingParameters.h"

#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/VariableNames.h"

#include "ufo/utils/metoffice/MetOfficeQCFlags.h"

namespace ufo {

  /// \brief Profile QC checker base class
  class ProfileCheckBase {
   public:
    explicit ProfileCheckBase(const ConventionalProfileProcessingParameters &options);

    virtual ~ProfileCheckBase() {}

    /// Run check
    virtual void runCheck(ProfileDataHandler &profileDataHandler) = 0;

    /// Fill variables in validator using a ProfileDataHandler.
    /// This function will only be called for subclasses of ProfileCheckBase
    /// whose runOnEntireSample() implementation returns false.
    /// If runOnEntireSample() returns true, the subclass needs
    /// to handle the storage of validation data on its own.
    /// For an example see the ProfileAveragePressure class.
    virtual void fillValidationData(ProfileDataHandler &profileDataHandler) {}

    /// Get result of check (default fail)
    virtual bool getResult() {return false;}

    /// Run this check on the entire sample?
    virtual bool runOnEntireSample() {return false;}

    /// This check requires HofX to have been calculated.
    virtual bool requiresHofX() {return false;}

    /// List of names of required GeoVaLs.
    virtual oops::Variables getGeoVaLNames() {return oops::Variables();}

    /// List of names of GeoVaLs used in check validation.
    virtual oops::Variables getValidationGeoVaLNames() {return oops::Variables();}

    /// List of names of required obs diagnostics.
    virtual oops::Variables getObsDiagNames() {return oops::Variables();}

   protected:  // functions
    /// Apply correction to vector of values
    template <typename T>
      void correctVector(const std::vector <T> &v1,
                         const std::vector <T> &v2,
                         std::vector <T> &vout)
      {
        ASSERT(v1.size() == v2.size());
        vout.assign(v1.size(), 0);
        std::transform(v1.begin(), v1.end(), v2.begin(), vout.begin(), std::plus<T>());
      }

    /// Set a QC flag on one profile level.
    /// This is the base case for one vector.
    template <typename T>
      void SetQCFlag(const int& flag,
                     const size_t& jlev,
                     std::vector <T> &vec)
      {
        if (vec.size() > jlev) vec[jlev] |= flag;
      }

    /// Set a QC flag on one profile level.
    /// This is the recursive case that accepts an arbitrary number of vectors
    /// using a variadic template.
    template <typename T, typename... Args>
      void SetQCFlag(const int& flag,
                     const size_t& jlev,
                     std::vector <T> &vec1,
                     Args&... vecs)
    {
      if (vec1.size() > jlev) vec1[jlev] |= flag;
      SetQCFlag(flag, jlev, vecs...);
    }

    /// Add an "OPS_" prefix to the names of variables that are used in the comparison
    /// with results from the Met Office OPS system.
    std::string addOPSPrefix(const std::string & fullname);

   protected:  // variables
    /// Configurable parameters
    const ConventionalProfileProcessingParameters &options_;

    /// Missing value (int)
    const int missingValueInt = util::missingValue<int>();

    /// Missing value (float)
    const float missingValueFloat = util::missingValue<float>();
  };

  /// Profile check factory
  class ProfileCheckFactory
  {
   public:
    static std::unique_ptr<ProfileCheckBase> create(const std::string&,
                                                    const ConventionalProfileProcessingParameters&);
    virtual ~ProfileCheckFactory() = default;
   protected:
    explicit ProfileCheckFactory(const std::string &);
   private:
    virtual std::unique_ptr<ProfileCheckBase> make
      (const ConventionalProfileProcessingParameters&) = 0;

    static std::map <std::string, ProfileCheckFactory*> & getMakers()
      {
        static std::map <std::string, ProfileCheckFactory*> makers_;
        return makers_;
      }
  };

  template<class T>
    class ProfileCheckMaker : public ProfileCheckFactory
    {
      virtual std::unique_ptr<ProfileCheckBase>
        make(const ConventionalProfileProcessingParameters &options)
      {
        return std::unique_ptr<ProfileCheckBase>(new T(options));
      }
     public:
      explicit ProfileCheckMaker(const std::string & name)
        : ProfileCheckFactory(name) {}
    };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKBASE_H_

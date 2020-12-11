/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckBase.h"

#include <map>
#include <string>

#include "oops/util/Logger.h"

namespace ufo {
  ProfileCheckBase::ProfileCheckBase(const ProfileConsistencyCheckParameters &options,
                                     ProfileDataHandler &profileDataHandler,
                                     ProfileCheckValidator &profileCheckValidator)
    : options_(options),
      profileDataHandler_(profileDataHandler),
      profileCheckValidator_(profileCheckValidator)
  {}

  ProfileCheckFactory::ProfileCheckFactory(const std::string & name)
  {
    if (getMakers().find(name) != getMakers().end())
      throw eckit::BadParameter(name + " already registered in ufo::ProfileCheckFactory.", Here());
    getMakers()[name] = this;
  }

  std::unique_ptr<ProfileCheckBase>
  ProfileCheckFactory::create(const std::string& name,
                              const ProfileConsistencyCheckParameters &options,
                              ProfileDataHandler &profileDataHandler,
                              ProfileCheckValidator &profileCheckValidator)
  {
    oops::Log::trace() << "ProfileCheckBase::create starting" << std::endl;
    typename std::map<std::string, ProfileCheckFactory*>::iterator jloc =
      getMakers().find(name);
    if (jloc == getMakers().end()) {
      std::string makerNameList;
      for (const auto& makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
      throw eckit::BadParameter(name + " does not exist in ufo::ProfileCheckFactory. "
                                "Possible values:" + makerNameList, Here());
    }
    std::unique_ptr<ProfileCheckBase> ptr = jloc->second->make(options,
                                                               profileDataHandler,
                                                               profileCheckValidator);
    oops::Log::trace() << "ProfileCheckBase::create done" << std::endl;
    return ptr;
  }
}  // namespace ufo


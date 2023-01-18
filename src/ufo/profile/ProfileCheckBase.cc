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
  ProfileCheckBase::ProfileCheckBase(const ConventionalProfileProcessingParameters &options)
    : options_(options)
  {}

  std::string ProfileCheckBase::addOPSPrefix(const std::string & fullname)
  {
    std::string varname;
    std::string groupname;
    ufo::splitVarGroup(fullname, varname, groupname);
    return groupname + std::string("/OPS_") + varname;
  }

  ProfileCheckFactory::ProfileCheckFactory(const std::string & name)
  {
    if (getMakers().find(name) != getMakers().end())
      throw eckit::BadParameter(name + " already registered in ufo::ProfileCheckFactory.", Here());
    getMakers()[name] = this;
  }

  std::unique_ptr<ProfileCheckBase>
  ProfileCheckFactory::create(const std::string& name,
                              const ConventionalProfileProcessingParameters &options)
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
    std::unique_ptr<ProfileCheckBase> ptr = jloc->second->make(options);
    oops::Log::trace() << "ProfileCheckBase::create done" << std::endl;
    return ptr;
  }
}  // namespace ufo


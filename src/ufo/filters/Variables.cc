/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/Variables.h"

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

#include "ufo/utils/SplitVarGroup.h"

// -----------------------------------------------------------------------------
namespace ufo {

// -----------------------------------------------------------------------------

Variables::Variables()
  : fullnames_(0) {
  oops::Log::trace() << "ufo::Variables created empty" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & conf, const std::string & allgroup)
  : fullnames_(0) {
  oops::Log::trace() << "ufo::Variables(config) start " << conf << std::endl;
  if (conf.has("variables")) {
    std::vector<std::string> fullnames, varnames, grpnames;
    std::string var, grp;
    conf.get("variables", fullnames);
    // get variable name and group, check that both are specified
    for (size_t jvar = 0; jvar < fullnames.size(); ++jvar) {
      splitVarGroup(fullnames[jvar], var, grp);
      if (var == "" || (allgroup == "" && grp == "")) {
        oops::Log::error() << "Both name and group should be specified for variable" << std::endl;
        ABORT("Both name and group should be specified for variable");
      }
      if (allgroup != "" && grp != "") {
        oops::Log::error() << "Group for variable should only be specified once, but two "
                           << "groups found: " << grp << " and " << allgroup << std::endl;
        ABORT("Group for variable should only be specified once");
      }
      if (grp == "") grp = allgroup;
      varnames.push_back(var);
      grpnames.push_back(grp);
    }
    // read channels if available
    if (conf.has("channels")) {
      std::string chlist = conf.getString("channels");
      std::set<int> channelset = oops::parseIntSet(chlist);
      std::vector<int> channels;
      std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels));
      // assuming the same channel subsetting applies to all variables
      for (size_t jvar = 0; jvar < varnames.size(); ++jvar) {
        for (size_t jch = 0; jch < channels.size(); ++jch) {
          fullnames_.push_back(varnames[jvar]+"_"+std::to_string(channels[jch])+"@"+
                               grpnames[jvar]);
        }
      }
    } else {
      for (size_t jvar = 0; jvar < varnames.size(); ++jvar) {
        fullnames_.push_back(varnames[jvar]+"@"+grpnames[jvar]);
      }
    }
  } else if (conf.has("variable")) {
    std::string fullname(conf.getString("variable"));
    std::string var, grp;
    splitVarGroup(fullname, var, grp);
    if (var == "" || (allgroup == "" && grp == "")) {
      oops::Log::error() << "Both name and group should be specified for variable" << std::endl;
      ABORT("Both name and group should be specified for variable");
    }
    if (allgroup != "" && grp != "") {
      oops::Log::error() << "Group for variable should only be specified once, but two "
                         << "groups found: " << grp << " and " << allgroup << std::endl;
      ABORT("Group for variable should only be specified once");
    }
    if (grp == "") grp = allgroup;
    // read channels if available
    if (conf.has("channel")) {
      int channel = conf.getInt("channel");
      fullnames_.push_back(var+"_"+std::to_string(channel)+"@"+grp);
    } else {
      fullnames_.push_back(fullname);
    }
  } else {
    oops::Log::error() << "No variables section in " << conf << std::endl;
    ABORT("No variables section in config, unable to initialize ufo::Variables");
  }
  oops::Log::trace() << "ufo::Variables(conf) done" << std::endl;
}

// -----------------------------------------------------------------------------

Variables & Variables::operator+=(const Variables & rhs) {
  fullnames_.insert(fullnames_.end(), rhs.fullnames_.begin(), rhs.fullnames_.end());
  return *this;
}

// -----------------------------------------------------------------------------

Variables & Variables::operator+=(const std::string & rhs) {
  fullnames_.push_back(rhs);
  return *this;
}

// -----------------------------------------------------------------------------

const std::string Variables::variable(const size_t jvar) const {
  ASSERT(jvar < fullnames_.size());
  std::string var, grp;
  splitVarGroup(fullnames_[jvar], var, grp);
  return var;
}

// -----------------------------------------------------------------------------

const std::string Variables::group(const size_t jvar) const {
  ASSERT(jvar < fullnames_.size());
  std::string var, grp;
  splitVarGroup(fullnames_[jvar], var, grp);
  return grp;
}

// -----------------------------------------------------------------------------

oops::Variables Variables::allFromGroup(const std::string & group) const {
  oops::Variables vars;
  for (size_t jj = 0; jj < fullnames_.size(); ++jj) {
    std::string var, grp;
    splitVarGroup(fullnames_[jj], var, grp);
    if (grp == group) vars.push_back(var);
  }
  return vars;
}

// -----------------------------------------------------------------------------

bool Variables::has(const std::string & var) const {
  bool found = false;
  for (size_t jj = 0; jj < fullnames_.size(); ++jj) {
    found = found || fullnames_[jj] == var;
  }
  return found;
}

// -----------------------------------------------------------------------------

size_t Variables::find(const std::string & var) const {
  size_t ii = fullnames_.size();
  for (size_t jj = 0; jj < fullnames_.size(); ++jj) {
    if (fullnames_[jj] == var) ii = jj;
  }
  ASSERT(ii < fullnames_.size());
  return ii;
}

// -----------------------------------------------------------------------------

void Variables::removeDuplicates() {
  std::sort(fullnames_.begin(), fullnames_.end());
  fullnames_.erase(std::unique(fullnames_.begin(), fullnames_.end() ), fullnames_.end());
}

// -----------------------------------------------------------------------------

Variables::~Variables() {}

// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  os << fullnames_.size() << " variables: ";
  for (size_t jj = 0; jj < fullnames_.size(); ++jj) {
    if (jj > 0) os << ", ";
    os << fullnames_[jj];
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

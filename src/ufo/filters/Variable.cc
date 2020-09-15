/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/Variable.h"

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

#include "ufo/utils/StringUtils.h"


// -----------------------------------------------------------------------------
namespace ufo {

// -----------------------------------------------------------------------------

Variable::Variable(const eckit::Configuration & conf)
  : varname_(), grpname_(), channels_(),
    options_(conf.getSubConfiguration("options")) {
  oops::Log::trace() << "ufo::Variable(config) start " << conf << std::endl;
  std::string fullname;
  conf.get("name", fullname);
  splitVarGroup(fullname, varname_, grpname_);
  // read channels if available
  if (conf.has("channels")) {
    std::string chlist = conf.getString("channels");
    std::set<int> channelset = oops::parseIntSet(chlist);
    std::copy(channelset.begin(), channelset.end(), std::back_inserter(channels_));
  }
  oops::Log::trace() << "ufo::Variable(conf) done" << std::endl;
}

// -----------------------------------------------------------------------------

Variable::Variable(const std::string & fullname, const eckit::LocalConfiguration options)
  : varname_(), grpname_(), channels_(), options_(options) {
  oops::Log::trace() << "ufo::Variable(name) start" << std::endl;
  splitVarGroup(fullname, varname_, grpname_);
  oops::Log::trace() << "ufo::Variable(name) done" << std::endl;
}

// -----------------------------------------------------------------------------

Variable::Variable(const std::string & fullname, const std::vector<int> & channels)
  : varname_(), grpname_(), channels_(channels), options_() {
  oops::Log::trace() << "ufo::Variable(name, channels) start " << std::endl;
  splitVarGroup(fullname, varname_, grpname_);
  oops::Log::trace() << "ufo::Variable(name, channels) done" << std::endl;
}

// -----------------------------------------------------------------------------

Variable::Variable(const Variable & var, const std::string & group)
  : varname_(var.varname_), grpname_(group), channels_(var.channels_),
    options_(var.options_) {
}

// -----------------------------------------------------------------------------

Variable::~Variable() {
}

// -----------------------------------------------------------------------------

size_t Variable::size() const {
  if (channels_.size() == 0) {
    return 1;
  } else {
    return channels_.size();
  }
}

// -----------------------------------------------------------------------------

Variable Variable::operator[](const size_t jch) const {
  ASSERT(jch < this->size());
  std::string var = varname_;
  if (channels_.size() > 0)  var += "_" + std::to_string(channels_[jch]);
  if (grpname_ != "")        var += "@" + grpname_;
  return Variable(var, options_);
}


// -----------------------------------------------------------------------------

const std::string & Variable::variable() const {
  return varname_;
}

// -----------------------------------------------------------------------------

std::string Variable::variable(const size_t jch) const {
  ASSERT(jch < this->size());
  if (jch == 0 && channels_.size() == 0) {
    return varname_;
  } else {
    return varname_ + "_" + std::to_string(channels_[jch]);
  }
}

// -----------------------------------------------------------------------------

const std::string & Variable::group() const {
  return grpname_;
}

// -----------------------------------------------------------------------------

const std::vector<int> & Variable::channels() const {
  return channels_;
}

// -----------------------------------------------------------------------------

oops::Variables Variable::toOopsVariables() const {
  oops::Variables vars;
  for (size_t jj = 0; jj < this->size(); ++jj) {
    vars.push_back(this->variable(jj));
  }
  return vars;
}


// -----------------------------------------------------------------------------

void Variable::print(std::ostream & os) const {
  if (channels_.size() == 0) {
    os << varname_ + "@" + grpname_;
  } else {
    for (size_t jj = 0; jj < channels_.size(); ++jj) {
      if (jj > 0) os << ", ";
      os << varname_ + "_" + std::to_string(channels_[jj]) + "@" + grpname_;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

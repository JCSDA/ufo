/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/Variables.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"

// -----------------------------------------------------------------------------
namespace ufo {

// -----------------------------------------------------------------------------

Variables::Variables()
  : vars_() {
  oops::Log::trace() << "ufo::Variables created empty" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const std::vector<eckit::LocalConfiguration> & confs)
  : vars_() {
  oops::Log::trace() << "ufo::Variables(config) start " << std::endl;
  for (size_t jvar = 0; jvar < confs.size(); ++jvar) {
    vars_.push_back(Variable(confs[jvar]));
  }
  oops::Log::trace() << "ufo::Variables(conf) done" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const oops::ObsVariables & oopsvars)
  : vars_() {
  oops::Log::trace() << "ufo::Variables(oops::Vars) start" << std::endl;
  if (oopsvars.channels().size() > 0) {
    // note: for variables with channels will use only the first variable
    // should work for all current cases
    // find string before channel:
    size_t pos = oopsvars[0].find_last_of('_');
    vars_.push_back(Variable(oopsvars[0].substr(0, pos), oopsvars.channels()));
  } else {
    for (size_t jvar = 0; jvar < oopsvars.size(); ++jvar) {
      vars_.push_back(Variable(oopsvars[jvar]));
    }
  }
}

// -----------------------------------------------------------------------------

Variables::Variables(const oops::Variables & oopsvars)
  : vars_() {
  oops::Log::trace() << "ufo::Variables(oops::Vars) start" << std::endl;
  for (size_t jvar = 0; jvar < oopsvars.size(); ++jvar) {
    vars_.push_back(Variable(oopsvars[jvar].name()));
  }
  oops::Log::trace() << "ufo::Variables(oops::Vars) end" << std::endl;
}

// -----------------------------------------------------------------------------

Variables::Variables(const ufo::Variables & vars, const std::string & group)
  : vars_() {
  oops::Log::trace() << "ufo::Variables(ufovars, group) start " << std::endl;
  for (size_t jvar = 0; jvar < vars.size(); ++jvar) {
    vars_.push_back(Variable(vars[jvar], group));
  }
}

// -----------------------------------------------------------------------------

Variables::Variables(const std::vector<ufo::Variable> & vars)
  : vars_(vars) {
  oops::Log::trace() << "ufo::Variables(std::vector<ufo::Variable>) start " << std::endl;
}

// -----------------------------------------------------------------------------

Variables::~Variables() {
}

// -----------------------------------------------------------------------------

Variables & Variables::operator+=(const Variables & rhs) {
  vars_.insert(vars_.end(), rhs.vars_.begin(), rhs.vars_.end());
  return *this;
}

// -----------------------------------------------------------------------------

Variables & Variables::operator+=(const Variable & rhs) {
  vars_.push_back(rhs);
  return *this;
}

// -----------------------------------------------------------------------------

size_t Variables::size() const {
  return vars_.size();
}

// -----------------------------------------------------------------------------

const Variable & Variables::operator[](const size_t jj) const {
  return vars_[jj];
}

// -----------------------------------------------------------------------------

size_t Variables::nvars() const {
  size_t nvars = 0;
  for (size_t ivar = 0; ivar < vars_.size(); ++ivar) {
    nvars += vars_[ivar].size();
  }
  return nvars;
}

// -----------------------------------------------------------------------------

Variable Variables::variable(const size_t jj) const {
  size_t curr_indx = 0;
  for (size_t ivar = 0; ivar < vars_.size(); ++ivar) {
    if (jj < curr_indx + vars_[ivar].size())
      return vars_[ivar][jj-curr_indx];
    else
      curr_indx += vars_[ivar].size();
  }
  ABORT("Variable index exceeds collective variable arrays size");
  abort(); /* Prevent g++ warning of missing return */
}

// -----------------------------------------------------------------------------

Variables Variables::allFromGroup(const std::string & group) const {
  Variables vars;
  for (size_t ivar = 0; ivar < vars_.size(); ++ivar) {
    if (vars_[ivar].group() == group) {
      vars += vars_[ivar];
    } else if (vars_[ivar].group() == ObsFunctionTraits<float>::groupName) {
      ObsFunction<float> obsfunc(vars_[ivar]);
      ufo::Variables funcvars = obsfunc.requiredVariables();
      vars += funcvars.allFromGroup(group);
    } else if (vars_[ivar].group() == ObsFunctionTraits<int>::groupName) {
      ObsFunction<int> obsfunc(vars_[ivar]);
      ufo::Variables funcvars = obsfunc.requiredVariables();
      vars += funcvars.allFromGroup(group);
    } else if (vars_[ivar].group() == ObsFunctionTraits<std::string>::groupName) {
      ObsFunction<std::string> obsfunc(vars_[ivar]);
      ufo::Variables funcvars = obsfunc.requiredVariables();
      vars += funcvars.allFromGroup(group);
    } else if (vars_[ivar].group() == ObsFunctionTraits<util::DateTime>::groupName) {
      ObsFunction<util::DateTime> obsfunc(vars_[ivar]);
      ufo::Variables funcvars = obsfunc.requiredVariables();
      vars += funcvars.allFromGroup(group);
    }
  }
  return vars;
}

// -----------------------------------------------------------------------------

oops::ObsVariables Variables::toOopsObsVariables() const {
  oops::ObsVariables vars;
  for (size_t ivar = 0; ivar < vars_.size(); ++ivar) {
    for (size_t jj = 0; jj < vars_[ivar].size(); ++jj) {
      vars.push_back(vars_[ivar].variable(jj));
    }
  }
  return vars;
}

// -----------------------------------------------------------------------------

oops::Variables Variables::toOopsVariables() const {
  oops::Variables vars;
  for (size_t ivar = 0; ivar < vars_.size(); ++ivar) {
    for (size_t jj = 0; jj < vars_[ivar].size(); ++jj) {
      vars.push_back(vars_[ivar].variable(jj));
    }
  }
  return vars;
}

// -----------------------------------------------------------------------------

bool Variables::hasGroup(const std::string & group) const {
  bool found = false;
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (vars_[jj].group() == group)
      found = true;

    if (vars_[jj].group() == ObsFunctionTraits<float>::groupName) {
      ObsFunction<float> obsfunc(vars_[jj]);
      ufo::Variables funcvars = obsfunc.requiredVariables();
      found = found || funcvars.hasGroup(group);
    } else if (vars_[jj].group() == ObsFunctionTraits<int>::groupName) {
      ObsFunction<int> obsfunc(vars_[jj]);
      ufo::Variables funcvars = obsfunc.requiredVariables();
      found = found || funcvars.hasGroup(group);
    } else if (vars_[jj].group() == ObsFunctionTraits<std::string>::groupName) {
      ObsFunction<std::string> obsfunc(vars_[jj]);
      ufo::Variables funcvars = obsfunc.requiredVariables();
      found = found || funcvars.hasGroup(group);
    } else if (vars_[jj].group() == ObsFunctionTraits<util::DateTime>::groupName) {
      ObsFunction<util::DateTime> obsfunc(vars_[jj]);
      ufo::Variables funcvars = obsfunc.requiredVariables();
      found = found || funcvars.hasGroup(group);
    }
  }
  return found;
}

// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  os << vars_.size() << " variables: ";
  for (size_t jj = 0; jj < vars_.size(); ++jj) {
    if (jj > 0) os << ", ";
    os << vars_[jj];
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo

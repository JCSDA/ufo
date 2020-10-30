/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_VARIABLES_H_
#define UFO_FILTERS_VARIABLES_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"
#include "ufo/filters/Variable.h"

namespace ufo {

// -----------------------------------------------------------------------------

class Variables: public util::Printable {
 public:
  static const std::string classname() {return "ufo::Variables";}

  Variables();
  explicit Variables(const std::vector<eckit::LocalConfiguration> &);
  explicit Variables(const oops::Variables &);
  Variables(const ufo::Variables &, const std::string &);
  ~Variables();

  Variables & operator+=(const Variables &);
  Variables & operator+=(const Variable &);

  size_t size() const;
  const Variable & operator[](const size_t) const;

// the below two functions are for compatibility with oops::Variables and should
// eventually be removed
  size_t nvars() const;
  Variable variable(const size_t) const;

  Variables allFromGroup(const std::string &) const;
  oops::Variables toOopsVariables() const;

  bool hasGroup(const std::string &) const;
  operator bool() const {return !vars_.empty();}

 private:
  void print(std::ostream &) const;

  std::vector<Variable> vars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLES_H_

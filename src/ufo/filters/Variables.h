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

#include "oops/util/Printable.h"
#include "ufo/filters/Variable.h"

namespace eckit {
  class LocalConfiguration;
}

namespace oops {
  class ObsVariables;
}

namespace ufo {

// -----------------------------------------------------------------------------

class Variables: public util::Printable {
 public:
  static const std::string classname() {return "ufo::Variables";}

  Variables();
  explicit Variables(const std::vector<eckit::LocalConfiguration> &);
  explicit Variables(const oops::Variables &);
  explicit Variables(const oops::ObsVariables &);
  Variables(const ufo::Variables &, const std::string &);
  explicit Variables(const std::vector<Variable> &);
  ~Variables();

  Variables & operator+=(const Variables &);
  Variables & operator+=(const Variable &);

  /// \brief Return the number of constituent Variable objects (some of which may contain multiple
  /// channels).
  size_t size() const;
  /// \brief Return a given constituent Variable (which may contain multiple channels).
  const Variable & operator[](const size_t) const;

// the below two functions are for compatibility with oops::Variables and should
// eventually be removed
  /// \brief Return the number of constituent "primitive" (single-channel) variables.
  size_t nvars() const;
  /// \brief Return a given constituent "primitive" (single-channel) variable.
  Variable variable(const size_t) const;

  Variables allFromGroup(const std::string &) const;
  oops::ObsVariables toOopsObsVariables() const;
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

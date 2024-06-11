/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_VARIABLE_H_
#define UFO_FILTERS_VARIABLE_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/Printable.h"

namespace oops {
  class ObsVariables;
  class Variables;
  class Variable;
}

namespace ufo {

// -----------------------------------------------------------------------------

class Variable: public util::Printable {
 public:
  static const std::string classname() {return "ufo::Variable";}

  explicit Variable(const eckit::Configuration &);
  explicit Variable(const std::string &,
                    const eckit::LocalConfiguration conf = eckit::LocalConfiguration());
  Variable(const std::string &, const std::vector<int> &);
  Variable(const Variable &, const std::string &);
  ~Variable();

  size_t size() const;
  Variable operator[](const size_t) const;
  const std::string & variable() const;
  oops::Variable toOopsVariable() const;
  std::string variable(const size_t) const;
  const std::string & group() const;
  const std::vector<int> & channels() const;

  /// Return the full variable name including the group name (but no channel number), e.g.
  /// `MetaData/pressure`.
  std::string fullName() const;

  oops::Variables toOopsVariables() const;
  oops::ObsVariables toOopsObsVariables() const;
  const eckit::LocalConfiguration & options() const {return options_;}

 private:
  void print(std::ostream &) const;
  std::string varname_;
  std::string grpname_;
  std::vector<int> channels_;
  eckit::LocalConfiguration options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLE_H_

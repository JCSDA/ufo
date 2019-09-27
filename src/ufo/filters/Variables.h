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

namespace ufo {

// -----------------------------------------------------------------------------

class Variables: public util::Printable {
 public:
  static const std::string classname() {return "ufo::Variables";}

  Variables();
  explicit Variables(const eckit::Configuration &, const std::string & allgroup = "");
  Variables(const std::string &, const std::vector<int> &);
  ~Variables();

  Variables & operator+=(const Variables &);
  Variables & operator+=(const std::string &);

  size_t size() const {return fullnames_.size();}
  const std::string & operator[](const size_t kk) const {return fullnames_.at(kk);}

  const std::string variable(const size_t) const;
  const std::string group(const size_t) const;
  oops::Variables allFromGroup(const std::string &) const;

  bool has(const std::string &) const;
  bool hasGroup(const std::string &) const;
  size_t find(const std::string &) const;
  operator bool() const {return !fullnames_.empty();}

  void removeDuplicates();

 private:
  void print(std::ostream &) const;

  std::vector<std::string> fullnames_;  // full variable names (with channel and @)
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLES_H_

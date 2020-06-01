/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_STRINGUTILS_H_
#define UFO_UTILS_STRINGUTILS_H_

#include <set>
#include <string>

namespace ufo {
  void splitVarGroup(const std::string &, std::string &, std::string &);
  void splitInstSat(const std::string &, std::string &, std::string &);
  bool isFloat(const std::string &);
  bool readFloat(const std::string &, float &);
}  // namespace ufo

#endif  // UFO_UTILS_STRINGUTILS_H_

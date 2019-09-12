/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PROCESSWHERE_H_
#define UFO_FILTERS_PROCESSWHERE_H_

#include <set>
#include <string>
#include <vector>

namespace eckit {class Configuration;}

namespace ufo {
  class ObsFilterData;
  class Variables;

ufo::Variables getAllWhereVariables(const eckit::Configuration &);
std::vector<bool> processWhere(const eckit::Configuration &, const ObsFilterData &);

}  // namespace ufo

#endif  // UFO_FILTERS_PROCESSWHERE_H_

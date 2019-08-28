/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PROCESSWHERE_H_
#define UFO_FILTERS_PROCESSWHERE_H_

#include <string>
#include <vector>

namespace eckit {class Configuration;}
namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}
namespace oops {class Variables;}

namespace ufo {
  class GeoVaLs;

oops::Variables preProcessWhere(const eckit::Configuration &, const std::string &);
std::vector<bool> processWhere(ioda::ObsSpace &, const GeoVaLs &,
                               const eckit::Configuration &);
std::vector<bool> processWhere(ioda::ObsSpace &, const ioda::ObsDataVector<float> &,
                               const eckit::Configuration &);
}  // namespace ufo

#endif  // UFO_FILTERS_PROCESSWHERE_H_

/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_COPYVARS2ODV_H_
#define UFO_FILTERS_COPYVARS2ODV_H_

#include <string>
#include <vector>

namespace eckit {class Configuration;}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
  class ObsVector;
}

namespace oops {class Variables;}

namespace ufo {
  class GeoVaLs;

const oops::Variables collectFilterVars(const eckit::Configuration &);
ioda::ObsDataVector<float> copyVars2ODV(const GeoVaLs &, ioda::ObsDataVector<float> &);
ioda::ObsDataVector<float> copyVars2ODV(const ioda::ObsVector &, ioda::ObsDataVector<float> &,
                                        const std::string &);
ioda::ObsDataVector<float> copyVars2ODV(ioda::ObsSpace &, ioda::ObsDataVector<float> &);
}  // namespace ufo

#endif  // UFO_FILTERS_COPYVARS2ODV_H_

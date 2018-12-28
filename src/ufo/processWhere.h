/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROCESSWHERE_H_
#define UFO_PROCESSWHERE_H_

#include <vector>

namespace eckit {class Configuration;}
namespace ioda {class ObsSpace;}

namespace ufo {

std::vector<bool> processWhere(ioda::ObsSpace &, const eckit::Configuration &);

}  // namespace ufo

#endif  // UFO_PROCESSWHERE_H_

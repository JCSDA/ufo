/*
 * (C) Copyright 2019 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */


#ifndef UFO_OPERATORS_TIMEOPER_OBSTIMEOPERUTIL_H_
#define UFO_OPERATORS_TIMEOPER_OBSTIMEOPERUTIL_H_

#include <algorithm>
#include <ostream>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace ufo {

class ObsTimeOperParameters;

std::vector<std::vector<float>> timeWeightCreate(const ioda::ObsSpace & odb_,
                                                 const ObsTimeOperParameters & parameters);

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_TIMEOPER_OBSTIMEOPERUTIL_H_

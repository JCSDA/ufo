/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_CONSTANTS_H_
#define UFO_UTILS_CONSTANTS_H_

#include <cmath>

//------------------------------------------------------------------------------------------------------

namespace ufo {

//------------------------------------------------------------------------------------------------------

/// Some useful constants
struct Constants {
    static constexpr double deg2rad() { return M_PI / 180.; }
    static constexpr double rad2deg() { return 180. * M_1_PI; }
};

//------------------------------------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_UTILS_CONSTANTS_H_

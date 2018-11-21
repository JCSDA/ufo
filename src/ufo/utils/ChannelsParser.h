/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_CHANNELSPARSER_H_
#define UFO_UTILS_CHANNELSPARSER_H_

#include <string>
#include <vector>

namespace ufo {

  std::vector<int> parseChannels(const std::string&);

}  // namespace ufo

#endif  // UFO_UTILS_CHANNELSPARSER_H_

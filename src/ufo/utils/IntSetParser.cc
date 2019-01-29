/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/IntSetParser.h"

#include <algorithm>

#include <sstream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

namespace ufo {

// -----------------------------------------------------------------------------

/// Function to split string on delimiter

std::vector<std::string> splitString(const std::string& str, char delim)
{
  std::vector<std::string> result;
  std::stringstream s(str);
  std::string substr;
  while (std::getline(s, substr, delim)) {
    result.push_back(substr);
  }
  return result;
}

// -----------------------------------------------------------------------------

/// Function to parse channels (supports commas for separating channels
//  and channel ranges and dashes for channel ranges).
//  For example: 1-5, 9, 13-45
//  Returns a std::set, no need to sort or remove duplicates and find/insert are in log(n)

std::set<int> parseIntSet(const std::string & str) {
  std::set<int> channels;

// split string by commas to get individual channels or ranges
  std::vector<std::string> ranges = splitString(str, ',');

  for (int irange = 0; irange < ranges.size(); irange++) {
    // split the element by dashes (in case it is a range)
    std::vector<std::string> range = splitString(ranges[irange], '-');
    ASSERT((range.size() == 1) || (range.size() == 2));
    // add a single channel
    if (range.size() == 1) {
      // add a single channel
      channels.insert(std::stoi(range[0]));
    } else if (range.size() == 2) {
      // add a range
      int start = std::stoi(range[0]);
      int stop  = std::stoi(range[1]);
      for (int ch = start; ch <= stop; ch++) {
        channels.insert(ch);
      }
    }
  }

  return channels;
}

// -----------------------------------------------------------------------------

void splitVarGroup(const std::string & vargrp, std::string & var, std::string & grp) {
  const size_t at = vargrp.find("@");
  var = vargrp.substr(0, at);
  if (at != std::string::npos) {
    grp = vargrp.substr(at + 1, std::string::npos);
    const size_t no_at = grp.find("@");
    ASSERT(no_at == std::string::npos);
  }
}

// -----------------------------------------------------------------------------
}  // namespace ufo

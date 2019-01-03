/*
 * (C) Copyright 2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <set>
#include <vector>

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ufo/utils/IntSetParser.h"

int main(int argc,  char ** argv) {
// Get configuration file from command line
  ASSERT(argc >= 2);
  eckit::PathName configfile = argv[argc - 1];

// Read configuration
  eckit::YAMLConfiguration config(configfile);

  oops::Log::info() << "Configuration input file is: " << configfile << std::endl;
  oops::Log::info() << "Full configuration is:"  << config << std::endl;

// Read channels list
  std::string chlist = config.getString("channels");
  std::set<int> channels = ufo::parseIntSet(chlist);
// Read expected output of parseChannels
  std::vector<int> expected = config.getIntVector("parsed_channels");

  ASSERT(std::equal(expected.begin(), expected.end(), channels.begin()));
  ASSERT(std::equal(channels.begin(), channels.end(), expected.begin()));
  return 0;
};


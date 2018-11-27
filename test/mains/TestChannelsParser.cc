/*
 * (C) Copyright 2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */


#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ufo/utils/ChannelsParser.h"

int main(int argc,  char ** argv) {

// Get configuration file from command line
  ASSERT(argc >= 2);
  eckit::PathName configfile = argv[argc - 1];

// Read configuration
  eckit::YAMLConfiguration config(configfile);

  oops::Log::info() << "Configuration input file is: " << configfile << std::endl;
  oops::Log::info() << "Full configuration is:"  << config << std::endl;

  std::string chlist = config.getString("channels"); //"1, 7-11, 22, 5-8";
  std::vector<int> channels = ufo::parseChannels(chlist);
  std::vector<int> expected = config.getIntVector("parsed_channels"); //{1, 5, 6, 7, 8, 9, 10, 11, 22};
  ASSERT(std::equal(expected.begin(), expected.end(), channels.begin()));
  return 0;
};


/*
 * (C) Copyright 2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "eckit/exception/Exceptions.h"

#include "ufo/utils/ChannelsParser.h"

int main(int argc,  char ** argv) {
  std::string s = "1, 7-11, 22, 5-8";
  std::vector<int> channels = ufo::parseChannels(s);
  std::vector<int> expected = {1, 5, 6, 7, 8, 9, 10, 11, 22};
  ASSERT(std::equal(expected.begin(), expected.end(), channels.begin()));
  return 0;
};


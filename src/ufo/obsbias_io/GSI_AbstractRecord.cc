/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/obsbias_io/GSI_AbstractRecord.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace ufo {

// -----------------------------------------------------------------------------

AbstractRecord::AbstractRecord()
  : vector_(gsi_predictors.size(), 0.0) {}

// -----------------------------------------------------------------------------

AbstractRecord::~AbstractRecord() {}

// -----------------------------------------------------------------------------

void AbstractRecord::setID(const std::size_t seq,
                           const std::string & sensor,
                           const std::size_t channel) {
  seq_ = seq;
  sensor_ = sensor;
  channel_ = channel;
}

// -----------------------------------------------------------------------------

void AbstractRecord::fillVector(const std::vector< std::string > & predictors,
                           const std::vector< double > & data) {
  for (std::size_t i = 0; i < predictors.size(); ++i) {
    const auto iter = std::find(gsi_predictors.begin(), gsi_predictors.end(), predictors[i]);
    if (iter != gsi_predictors.end()) {
      vector_.at(std::distance(gsi_predictors.begin(), iter)) = data[i];
    }
  }
}

// -----------------------------------------------------------------------------

std::vector< double >
AbstractRecord::readByChannel(std::fstream & inFile,
                              const std::string & sensor,
                              const int channel,
                              const std::vector< std::string > & predictors) {
  inFile.clear();
  inFile.seekg(0, std::ios::beg);

  while (readNext(inFile)) {
    if (sensor == sensor_ && channel == channel_) {
      break;
    }
  }

  std::vector< double > data(predictors.size(), 0.0);

  for (std::size_t i = 0; i < gsi_predictors.size(); ++i) {
    const auto iter = std::find(predictors.begin(), predictors.end(), gsi_predictors[i]);
    if (iter != predictors.end()) {
      data.at(std::distance(predictors.begin(), iter)) = vector_[i];
    }
  }

  return data;
}

// -----------------------------------------------------------------------------

const std::string & AbstractRecord::getSensorID() const {
  return sensor_;
}

// -----------------------------------------------------------------------------

const std::size_t AbstractRecord::getChannel() const {
  return channel_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/obsbias_io/GSI_Record.h"

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

Record::Record()
  : biasCoeffs_(gsi_predictors.size(), 0.0),
    tlap_(util::missingValue(0.0)),
    tsum_(util::missingValue(0.0)),
    ntlapupdate_(util::missingValue(1)) {}

// -----------------------------------------------------------------------------

Record::~Record() {}

// -----------------------------------------------------------------------------

void Record::setID(const std::size_t seq,
           const std::string & sensor,
           const std::size_t channel) {
  seq_ = seq;
  sensor_ = sensor;
  channel_ = channel;
}

// -----------------------------------------------------------------------------

void Record::setPredictors(const std::vector< std::string > & names,
                   const std::vector< double > & data) {
  for (std::size_t i = 0; i < names.size(); ++i) {
    const auto iter = std::find(gsi_predictors.begin(), gsi_predictors.end(), names[i]);
    if (iter != gsi_predictors.end()) {
      biasCoeffs_.at(std::distance(gsi_predictors.begin(), iter)) = data[i];
    }
  }
}

// -----------------------------------------------------------------------------

void Record::setValue(const std::string & name,
              const double value) {
  if (name == std::string("tlap")) {
    tlap_ = value;
  } else if (name == std::string("tsum")) {
    tsum_ = value;
  } else if (name == std::string("ntlapupdate")) {
    ntlapupdate_ = value;
  } else {
    oops::Log::error() << "Unable to set value for "
                       << name << std::endl;
    ABORT("Unable to set value for " + name);
  }
}

// -----------------------------------------------------------------------------

bool Record::readNext(std::fstream & inFile) {
  if (inFile >> seq_ &&
      inFile >> sensor_ &&
      inFile >> channel_ &&
      inFile >> tlap_ &&
      inFile >> tsum_ &&
      inFile >> ntlapupdate_) {
    for (std::size_t i = 0; i < gsi_predictors.size(); ++i) {
      inFile >> biasCoeffs_.at(i);
    }
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------

std::vector< double >
Record::readByChannels(std::fstream & inFile,
               const std::string & sensor,
               const std::vector< int > & channels,
               const std::string & name) {
  inFile.clear();
  inFile.seekg(0, std::ios::beg);

  std::vector< double > data(channels.size(), 0.0);

  while (readNext(inFile)) {
    const auto iter = std::find(channels.begin(), channels.end(), channel_);
    if (sensor == sensor_ && iter != channels.end()) {
      if (name == std::string("tlap")) {
        data.at(std::distance(channels.begin(), iter)) = tlap_;
      } else if (name == std::string("tsum")) {
        data.at(std::distance(channels.begin(), iter)) = tsum_;
      } else if (name == std::string("ntlapupdate")) {
        data.at(std::distance(channels.begin(), iter)) = ntlapupdate_;
      } else {
        oops::Log::error() << "Unable to read data of "
                           << name << std::endl;
        ABORT("Unable to read data of " + name);
      }
    }
  }
  return data;
}

// -----------------------------------------------------------------------------

std::vector< double >
Record::readByPredictors(std::fstream & inFile,
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
      data.at(std::distance(predictors.begin(), iter)) = biasCoeffs_[i];
    }
  }

  return data;
}

// -----------------------------------------------------------------------------

const std::string & Record::getSensorID() const {
  return sensor_;
}

// -----------------------------------------------------------------------------

const std::size_t Record::getChannel() const {
  return channel_;
}

// -----------------------------------------------------------------------------

void Record::writeTo(std::fstream & outFile) {
  outFile << seq_ << " ";
  outFile << sensor_ << " ";
  outFile << channel_ << " ";
  outFile << tlap_ << " ";
  outFile << tsum_ << " ";
  outFile << ntlapupdate_;
  outFile << std::endl;

  for (const auto & item : biasCoeffs_) {
    outFile << item << "  ";
  }
  outFile << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

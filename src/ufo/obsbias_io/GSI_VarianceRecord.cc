/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/obsbias_io/GSI_VarianceRecord.h"

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

VarianceRecord::VarianceRecord()
  : variances_(gsi_predictors.size(), 0.0),
    countOfQCedObs_(util::missingValue(0.0)) {}

// -----------------------------------------------------------------------------

VarianceRecord::~VarianceRecord() {}

// -----------------------------------------------------------------------------

void VarianceRecord::fillVector(const std::vector< std::string > & predictors,
                                const std::vector< double > & data) {
  for (std::size_t i = 0; i < predictors.size(); ++i) {
    const auto iter = std::find(gsi_predictors.begin(), gsi_predictors.end(), predictors[i]);
    if (iter != gsi_predictors.end()) {
      variances_.at(std::distance(gsi_predictors.begin(), iter)) = data[i];
    }
  }
}

// -----------------------------------------------------------------------------

void VarianceRecord::setValueByVarName(const std::string & name,
                                       const double value) {
  if (name == std::string("countOfQCedObs")) {
    countOfQCedObs_ = value;
  } else {
    oops::Log::error() << "Unable to set value for "
                       << name << std::endl;
    ABORT("Unable to set value for " + name);
  }
}

// -----------------------------------------------------------------------------

bool VarianceRecord::readNext(std::fstream & inFile) {
  if (inFile >> seq_ &&
      inFile >> sensor_ &&
      inFile >> channel_ &&
      inFile >> countOfQCedObs_) {
    for (std::size_t i = 0; i < gsi_predictors.size(); ++i) {
      inFile >> variances_.at(i);
    }
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------

std::vector< double >
VarianceRecord::readByVarName(std::fstream & inFile,
                              const std::string & sensor,
                              const std::vector< int > & channels,
                              const std::string & name) {
  inFile.clear();
  inFile.seekg(0, std::ios::beg);

  std::vector< double > data(channels.size(), 0.0);

  while (readNext(inFile)) {
    const auto iter = std::find(channels.begin(), channels.end(), channel_);
    if (sensor == sensor_ && iter != channels.end()) {
      if (name == std::string("countOfQCedObs")) {
        data.at(std::distance(channels.begin(), iter)) = countOfQCedObs_;
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
VarianceRecord::readByChannel(std::fstream & inFile,
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
      data.at(std::distance(predictors.begin(), iter)) = variances_[i];
    }
  }

  return data;
}

// -----------------------------------------------------------------------------

void VarianceRecord::writeTo(std::fstream & outFile) {
  outFile << seq_ << " ";
  outFile << sensor_ << " ";
  outFile << channel_ << " ";
  outFile << countOfQCedObs_;
  outFile << std::endl;

  for (const auto & item : variances_) {
    outFile << item << "  ";
  }
  outFile << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

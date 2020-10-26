/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIAS_IO_OBSBIASIO_H_
#define UFO_OBSBIAS_IO_OBSBIASIO_H_

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "ufo/obsbias_io/GSI_Record.h"
#include "ufo/obsbias_io/GSI_VarianceRecord.h"

namespace ufo {
// -----------------------------------------------------------------------------

template< typename T >
class ObsBiasIO {
 public:
  ObsBiasIO(const std::string & fileName,
            const std::ios_base::openmode & mode)
    : fileIO_{fileName, mode}, records_() {
    if (!fileIO_.is_open()) {
      oops::Log::error() << "Unable to open file to read/write "
                         << fileName << std::endl;
      ABORT("Unable to open file to read/write " + fileName);
    }
  }

// -----------------------------------------------------------------------------

  ~ObsBiasIO() {
    if (fileIO_.is_open()) {
      fileIO_.close();
      oops::Log::debug() << "ObsBiasIO:: read/write done" << std::endl;
    }
  }

// -----------------------------------------------------------------------------

  std::vector< double >
  readByVarName(const std::string & sensor,
                const std::vector< int > & channels,
                const std::string & name) {
    T record;
    return (record.readByVarName(fileIO_, sensor, channels, name));
  }

// -----------------------------------------------------------------------------

  std::vector< double >
  readByChannel(const std::string & sensor,
                const int channel,
                const std::vector< std::string > & predictors) {
    T record;
    return (record.readByChannel(fileIO_, sensor, channel, predictors));
  }

  void addByVarName(const std::string & sensor,
                    const std::vector< int > & channels,
                    const std::string & name,
                    const std::vector< double > & data) {
    assert(data.size() == channels.size());

    bool notFound = true;
    for (std::size_t i = 0; i < channels.size(); ++i) {
      for (auto & record : records_) {
        if (sensor == record.getSensorID() && channels[i] == record.getChannel()) {
          record.setValueByVarName(name, data[i]);
          notFound = false;
          break;
        }
      }
      if (notFound) {
        records_.emplace_back();
        records_.back().setID(records_.size(), sensor, channels[i]);
        records_.back().setValueByVarName(name, data[i]);
      }
    }
  }

// -----------------------------------------------------------------------------

  void addByChannel(const std::string & sensor,
                    const int channel,
                    const std::vector< std::string > & predictors,
                    const std::vector< double > & data) {
    bool notFound = true;
    for (auto & record : records_) {
      if (sensor == record.getSensorID() && channel == record.getChannel()) {
        record.fillVector(predictors, data);
        notFound = false;
        break;
      }
    }

    if (notFound) {
      records_.emplace_back();
      records_.back().setID(records_.size(), sensor, channel);
      records_.back().fillVector(predictors, data);
    }
  }

// -----------------------------------------------------------------------------

  void commit() {
    for (auto & record : records_) {
      record.writeTo(fileIO_);
    }
  }

// -----------------------------------------------------------------------------

 private:
  std::fstream fileIO_;
  std::vector< T > records_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_IO_OBSBIASIO_H_

/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIAS_IO_GSI_VARIANCERECORD_H_
#define UFO_OBSBIAS_IO_GSI_VARIANCERECORD_H_

#include <fstream>
#include <string>
#include <vector>

#include "ufo/obsbias_io/GSI_Record.h"

namespace ufo {

// -----------------------------------------------------------------------------

class VarianceRecord: public Record {
 public:
  VarianceRecord();
  virtual ~VarianceRecord();

  void setPredictors(const std::vector< std::string > &,
                     const std::vector< double > &);

  void setValue(const std::string &, const double);

  std::vector< double > readByChannels(std::fstream &,
                                       const std::string &,
                                       const std::vector< int > &,
                                       const std::string &);

  std::vector< double > readByPredictors(std::fstream &,
                                         const std::string &,
                                         const int,
                                         const std::vector< std::string > &);

  void writeTo(std::fstream &);

 private:
  bool readNext(std::fstream &);
  double countOfQCedObs_;
  std::vector< double > variances_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_IO_GSI_VARIANCERECORD_H_

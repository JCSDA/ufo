/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIAS_IO_GSI_RECORD_H_
#define UFO_OBSBIAS_IO_GSI_RECORD_H_

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "ufo/obsbias_io/GSI_AbstractRecord.h"

namespace ufo {

// -----------------------------------------------------------------------------

class Record : public AbstractRecord {
 public:
  Record();
  virtual ~Record();

  void setValueByVarName(const std::string &,
                         const double) override;

  std::vector< double > readByVarName(std::fstream &,
                                      const std::string &,
                                      const std::vector< int > &,
                                      const std::string &) override;

  void writeTo(std::fstream &) override;

 private:
  bool readNext(std::fstream &) override;
  double tlap_, tsum_;
  std::size_t ntlapupdate_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_IO_GSI_RECORD_H_

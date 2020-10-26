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

namespace ufo {

// -----------------------------------------------------------------------------

// Default predictor names from GSI
// temporary solution, we should have a self-explanatory obsbias file
static const std::vector<std::string> gsi_predictors =
  {"constant",
   "zenith_angle",
   "cloud_liquid_water",
   "lapse_rate_order_2",
   "lapse_rate",
   "cosine_of_latitude_times_orbit_node",
   "sine_of_latitude",
   "emissivity",
   "scan_angle_order_4",
   "scan_angle_order_3",
   "scan_angle_order_2",
   "scan_angle"
  };

// -----------------------------------------------------------------------------

class Record {
 public:
  Record();
  virtual ~Record();

  void setID(const std::size_t,
             const std::string &,
             const std::size_t);

  void setPredictors(const std::vector< std::string > &,
                     const std::vector< double > &);

  void setValue(const std::string &,
                const double);

  std::vector< double >
  readByChannels(std::fstream &,
                 const std::string &,
                 const std::vector< int > &,
                 const std::string &);

  std::vector< double >
  readByPredictors(std::fstream &,
                   const std::string &,
                   const int,
                   const std::vector< std::string > &);

  const std::string & getSensorID() const;

  const std::size_t getChannel() const;

  void writeTo(std::fstream &);

 protected:
  std::size_t seq_;
  std::string sensor_;
  int channel_;

 private:
  bool readNext(std::fstream &);
  double tlap_, tsum_;
  std::size_t ntlapupdate_;
  std::vector< double > biasCoeffs_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_IO_GSI_RECORD_H_

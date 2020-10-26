/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSBIAS_IO_GSI_ABSTRACTRECORD_H_
#define UFO_OBSBIAS_IO_GSI_ABSTRACTRECORD_H_

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

class AbstractRecord {
 public:
  AbstractRecord();
  virtual ~AbstractRecord();

  void setID(const std::size_t,
             const std::string &,
             const std::size_t);

  void fillVector(const std::vector< std::string > &,
                  const std::vector< double > &);

  virtual void setValueByVarName(const std::string &,
                                 const double) = 0;

  virtual std::vector< double >
  readByVarName(std::fstream &,
                const std::string &,
                const std::vector< int > &,
                const std::string &) = 0;

  std::vector< double >
  readByChannel(std::fstream &,
                const std::string &,
                const int,
                const std::vector< std::string > &);

  const std::string & getSensorID() const;

  const std::size_t getChannel() const;

  virtual void writeTo(std::fstream &) = 0;

 protected:
  std::size_t seq_;
  std::string sensor_;
  int channel_;
  std::vector< double > vector_;

 private:
  virtual bool readNext(std::fstream &) = 0;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OBSBIAS_IO_GSI_ABSTRACTRECORD_H_

/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_METOFFICE_METOFFICEBMATRIXSTATIC_H_
#define UFO_UTILS_METOFFICE_METOFFICEBMATRIXSTATIC_H_

#include <Eigen/Dense>

#include <algorithm>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/utils/metoffice/MetOfficeBMatrixStatic.interface.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

// -----------------------------------------------------------------------------
/// MetOfficeBMatrixStatic: Met Office static model covariance
/// This class provides access to the static b matrix used for radiance
/// processing by the Met Office.  The objects main method is to multiply
/// an eigen matrix by the bmatrix
// -----------------------------------------------------------------------------

class MetOfficeBMatrixStatic : public util::Printable,
                               private util::ObjectCounter<MetOfficeBMatrixStatic> {
 public:
  static const std::string classname() {return "ufo::MetOfficeBMatrixStatic";}

  explicit MetOfficeBMatrixStatic(const eckit::Configuration &);

  size_t getindex(const float) const;
  size_t getsize(void) const;
  void multiply(const float, const Eigen::MatrixXf &, Eigen::MatrixXf &) const;
  void scale(const size_t elem, const float stdev);

 private:
  void print(std::ostream &) const override;
  F90obfilter keyMetOfficeBMatrixStatic_;  // key to Fortran for B
  size_t nbands_;                          // number of latitude bands for B
  size_t nelements_;                       // number of elements in each dimension of B
  std::vector<float> southlimits_;         // southern latitude limit per band
  std::vector<float> northlimits_;         // northern latitude limit per band
  std::vector<Eigen::MatrixXf> elements_;  // container for B contents
};

}  // namespace ufo

#endif  // UFO_UTILS_METOFFICE_METOFFICEBMATRIXSTATIC_H_

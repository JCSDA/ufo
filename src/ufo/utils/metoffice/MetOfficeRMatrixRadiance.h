/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_METOFFICE_METOFFICERMATRIXRADIANCE_H_
#define UFO_UTILS_METOFFICE_METOFFICERMATRIXRADIANCE_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/utils/metoffice/MetOfficeRMatrixRadiance.interface.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

// -----------------------------------------------------------------------------
/// MetOfficeRMatrixStatic: Met Office static model covariance
/// This class provides access to the static r matrix used for radiance
/// processing by the Met Office.
// -----------------------------------------------------------------------------

class MetOfficeRMatrixRadiance : public util::Printable,
                               private util::ObjectCounter<MetOfficeRMatrixRadiance> {
 public:
  static const std::string classname() {return "ufo::MetOfficeRMatrixRadiance";}

  explicit MetOfficeRMatrixRadiance(const eckit::Configuration &);

  void add(const std::vector<int> &, const Eigen::MatrixXf &, Eigen::MatrixXf &) const;

 private:
  void print(std::ostream &) const override;
  F90obfilter keyMetOfficeRMatrixRadiance_;
  size_t nchans_;
  size_t wmoid_;
  size_t rtype_;
  std::vector<int> channels_;
  std::vector<float> errors_;
};

}  // namespace ufo

#endif  // UFO_UTILS_METOFFICE_METOFFICERMATRIXRADIANCE_H_

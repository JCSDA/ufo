/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <assert.h>
#include <vector>

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "ufo/utils/metoffice/MetOfficeRMatrixRadiance.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// \brief Constructor
MetOfficeRMatrixRadiance::MetOfficeRMatrixRadiance(const eckit::Configuration & config):
    nchans_(0), wmoid_(0), rtype_(0), channels_(), errors_()
{
  oops::Log::trace() << "MetOfficeRMatrixRadiance constructor starting" << std::endl;

  // Read rmatrix file into Fortran object
  ufo_metoffice_rmatrixradiance_setup_f90(keyMetOfficeRMatrixRadiance_, config,
                                          nchans_, wmoid_, rtype_);

  // Only diagonal setup at the moment
  // OPS 1=full; 2=diagonal; 3=band diagonal
  if (rtype_ != 2) {
    ABORT("R-matrix type not currently in use - only diagonal");
  }

  // Map Fortran data to c++
  std::vector<int> chans_data(nchans_);
  std::vector<float> elements_data(nchans_);
  ufo_metoffice_rmatrixradiance_getelements_f90(keyMetOfficeRMatrixRadiance_, nchans_,
                                              chans_data.data(), elements_data.data());
  channels_ = chans_data;
  errors_ = elements_data;

  // Remove Fortran object as no longer needed
  ufo_metoffice_rmatrixradiance_delete_f90(keyMetOfficeRMatrixRadiance_);

  oops::Log::trace() << "MetOfficeRMatrixRadiance constructor end" << std::endl;
}
// -----------------------------------------------------------------------------
/// \brief Add r matrix variance onto input array
void MetOfficeRMatrixRadiance::add(const std::vector<int> & chans_used,
                                   const Eigen::MatrixXf & in,
                                   Eigen::MatrixXf & out) const {
  int matrows = in.rows();
  int matcols = in.cols();
  oops::Log::debug() << "matrows, matcols = " << matrows << " " << matcols << std::endl;
  assert(matrows == chans_used.size());
  assert(matcols == chans_used.size());
  out = in;
  if (rtype_ == 2) {
    for (size_t ichan = 0; ichan < chans_used.size(); ++ichan) {
      auto it = std::find(channels_.begin(), channels_.end(), chans_used[ichan]);
      if (it == channels_.end()) {
        oops::Log::error() << "Channel not found in R-matrix: "
                           << chans_used[ichan] << std::endl;
        ABORT("Invalid channel specified for R-matrix");
      } else {
        size_t index = it - channels_.begin();
        out(ichan, ichan) += errors_[index] * errors_[index];
      }
    }
  } else {
    ABORT("R-matrix type not currently in use - only diagonal");
  }
}
// -----------------------------------------------------------------------------
/// \brief Print
void MetOfficeRMatrixRadiance::print(std::ostream & os) const {
  os << "MetOfficeRMatrixRadiance: print starting" << std::endl;
  os << "nchans_ = " << nchans_ << std::endl;
  os << "wmoid_ = " << wmoid_ << std::endl;
  os << "rtype_ = " << rtype_ << std::endl;
  if (nchans_ < 30) {
    os << "channels = ";
    for (std::vector<int>::const_iterator i = channels_.begin(); i != channels_.end(); ++i)
        os << *i << ' ';
    os << std::endl;
    os << "errors (stdevs) = ";
    for (std::vector<float>::const_iterator i = errors_.begin(); i != errors_.end(); ++i)
        os << *i << ' ';
    os << std::endl;
  } else {
    os << "channels[0] = " << channels_[0] << std::endl;
    os << "errors[0] (stdev) = " << errors_[0] << std::endl;
  }
  os << "MetOfficeRMatrixRadiance: print end" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace ufo

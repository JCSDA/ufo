/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/metoffice/MetOfficeBMatrixStatic.h"
#include "oops/util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// \brief Constructor
MetOfficeBMatrixStatic::MetOfficeBMatrixStatic(const eckit::Configuration & config):
    nbands_(0), nelements_(0), southlimits_(), northlimits_(), elements_()
{
  oops::Log::trace() << "MetOfficeBMatrixStatic constructor starting" << std::endl;

  // Read bmatrix file into Fortran object
  ufo_metoffice_bmatrixstatic_setup_f90(keyMetOfficeBMatrixStatic_, config, nbands_, nelements_);

  // Map Fortran data to c++
  northlimits_.resize(nbands_);
  southlimits_.resize(nbands_);
  size_t ntotal = nelements_ * nelements_ * nbands_;
  std::vector<float> Btotal(ntotal);
  ufo_metoffice_bmatrixstatic_getelements_f90(keyMetOfficeBMatrixStatic_,
                                              nelements_, nbands_, southlimits_.data(),
                                              northlimits_.data(), Btotal.data());

  // For each of nbands, store B matrix in container
  for (size_t i = 0; i < nbands_; ++i) {
    Eigen::Map<Eigen::MatrixXf> bmap(Btotal.data()+i*nelements_*nelements_,
                                     nelements_, nelements_);
    elements_.push_back(bmap);
  }

  // Remove the Fortran object because it is no longer needed
  ufo_metoffice_bmatrixstatic_delete_f90(keyMetOfficeBMatrixStatic_);

  oops::Log::trace() << "MetOfficeBMatrixStatic constructor end" << std::endl;
}
// -----------------------------------------------------------------------------
/// \brief Return bmatrix size (number of rows or columns of square matrix)
size_t MetOfficeBMatrixStatic::getsize(void) const {
  return nelements_;
}
// -----------------------------------------------------------------------------
/// \brief Find bmatrix band index for a given latitude
size_t MetOfficeBMatrixStatic::getindex(const float latitude) const {
  auto lower = std::lower_bound(northlimits_.begin(), northlimits_.end(), latitude);
  return std::distance(northlimits_.begin(), lower);
}
// -----------------------------------------------------------------------------
/// \brief Multiply input matrix by bmatrix array based on latitude
void MetOfficeBMatrixStatic::multiply(const float lat,
                                      const Eigen::MatrixXf & in,
                                      Eigen::MatrixXf & out) const {
  size_t index = this->getindex(lat);
  out = elements_[index] * in;
}

// -----------------------------------------------------------------------------
/// \brief Scale elements of bmatrix array to user-defined standard deviation
void MetOfficeBMatrixStatic::scale(const size_t elem, const float stdev) {
  for (size_t iband = 0; iband < nbands_; ++iband) {
    float scaling = stdev/std::sqrt(elements_[iband](elem, elem));
    elements_[iband].row(elem) *= scaling;
    elements_[iband].col(elem) *= scaling;
  }
}
// -----------------------------------------------------------------------------
/// \brief Print
void MetOfficeBMatrixStatic::print(std::ostream & os) const {
  os << "MetOfficeBMatrixStatic: start print" << std::endl;
  os << "nbands_ = " << nbands_ << std::endl;
  os << "nelements_ = " << nelements_ << std::endl;
  os << "southlimits_[0] = " << southlimits_[0] << std::endl;
  os << "northlimits_[0] = " << northlimits_[0] << std::endl;
  os << "MetOfficeBMatrixStatic: end print" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace ufo

/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */
#include "ufo/obsbias/ObsBiasLinearCombination.h"

#include <Eigen/Core>
#include <fstream>

#include "ioda/ObsVector.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

static ObsBiasMaker<ObsBiasLinearCombination> makerBiasLinearCombination_("LinearCombination");

// -----------------------------------------------------------------------------

ObsBiasLinearCombination::ObsBiasLinearCombination(const eckit::Configuration & conf,
                                                   const std::vector<std::string> & prednames,
                                                   const std::vector<int> & jobs)
  : ObsBiasBase(conf), biascoeffs_(0), prednames_(prednames), jobs_(jobs) {
  oops::Log::trace() << "ObsBiasLinearCombination create.starting" << std::endl;

  // Initialize the biascoeffs
  biascoeffs_.resize(prednames_.size() * jobs_.size(), 0.0);

  oops::Log::trace() << "ObsBiasLinearCombination created." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombination::read(const std::string & sensor) {
  oops::Log::trace() << "ObsBiasLinearCombination::read from file: "
                     << input_filename_ << " Starting "<< std::endl;

// Default predictor names from GSI
  const std::vector<std::string> gsi_predictors = {"constant",
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
  std::ifstream infile(input_filename_);

  std::size_t ich;     //  sequential number
  std::string nusis;   //  sensor/instrument/satellite
  std::size_t nuchan;  //  channel number
  float tlap, tsum;
  std::size_t ntlapupdate;

  if (infile.is_open())
  {
    biascoeffs_.clear();
    float par;
    while (!infile.eof())
    {
      infile >> ich;
      infile >> nusis;
      infile >> nuchan;
      infile >> tlap;
      infile >> tsum;
      infile >> ntlapupdate;
      if (nusis == sensor &&
           std::find(jobs_.begin(), jobs_.end(), nuchan) != jobs_.end()) {
        for (auto & item : gsi_predictors) {
          infile >> par;
          if (std::find(prednames_.begin(), prednames_.end(), item)
              != prednames_.end()) {
            biascoeffs_.push_back(static_cast<double>(par));
          }
        }
      } else {
        for (auto item : gsi_predictors) {
          infile >> par;
        }
      }
    }
    infile.close();
    oops::Log::trace() << "ObsBiasLinearCombination::read from file: "
                      << input_filename_ << " Done " << std::endl;
  } else {
    oops::Log::error() << "Unable to open file : " << input_filename_ << std::endl;
    ABORT("Unable to open bias correction parameters file ");
  }

  oops::Log::trace() << "ObsBiasLinearCombination::read done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombination::write(const eckit::Configuration & conf) const {
  oops::Log::trace() << "ObsBiasLinearCombination::write to file : "
                     << output_filename_ << std::endl;

  oops::Log::trace() << "ObsBiasLinearCombination::write to file not implmented" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombination::computeObsBias(ioda::ObsVector & ybias,
                                              const Eigen::MatrixXd & predData) const {
  oops::Log::trace() << "ObsBiasLinearCombination::compute starting" << std::endl;
  const std::size_t nlocs  = ybias.nlocs();
  const std::size_t npreds = prednames_.size();
  const std::size_t njobs  = jobs_.size();

  ASSERT(biascoeffs_.size() == npreds*njobs);
  ASSERT(predData.rows() == npreds*njobs);
  ASSERT(predData.cols() == nlocs);
  ASSERT(ybias.nvars() == njobs);

  /* predData memory layout (njobs*npreds X nlocs)
   *     Loc     0      1      2       3
   *           --------------------------
   * ch1 pred1 | 0      1      2       4
   *     pred2 | 5      6      7       8
   *     pred3 | 9     10     11      12 
   * ch2 pred1 |13     14     15      16
   *     pred2 |17     18     19      20
   *     ....  |
   */

  ybias.zero();

  /* ybias memory layout (nlocs X njobs)
   *     ch1    ch2    ch3     ch4
   * Loc --------------------------
   *  0 | 0      1      2       3
   *  1 | 4      5      6       7
   *  2 | 8      9     10      11 
   *  3 |12     13     14      15
   *  4 |16     17     18      19
   * ...|
   */

  /* map bias coeff to eigen matrix npreds X njobs (read only)
   * bias coeff memory layout (njobs*npreds X nlocs)
   *        ch1    ch2    ch3     ch4
   *       --------------------------
   * pred1 | 0      1      2       4
   * pred2 | 5      6      7       8
   * pred3 | 9     10     11      12 
   * ....  |
   */
  Eigen::Map<const Eigen::MatrixXd> coeffs(biascoeffs_.data(), npreds, njobs);

  Eigen::VectorXd tmp;  // nlocs X 1
  for (std::size_t jch = 0; jch < njobs; ++jch) {
    //  ( nlocs X 1 ) =  ( nlocs X npreds ) * (  npreds X 1 )
    tmp = predData.transpose().block(0, jch*npreds , nlocs, npreds) *
          coeffs.block(0, jch, npreds, 1);
    for (std::size_t jrow = 0; jrow < nlocs; ++jrow) {
      ybias[jrow*njobs+jch] = tmp(jrow);
    }
  }

  oops::Log::trace() << "ObsBiasLinearCombination::compute done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombination::saveObsBiasTerms(ioda::ObsSpace & odb,
                                                const std::string & group,
                                                const Eigen::MatrixXd & predData) const {
  oops::Log::trace() << "ObsBias::saveObsBiasTerms startng." << std::endl;
  const std::size_t nlocs = odb.nlocs();
  const std::size_t npreds = prednames_.size();
  const std::size_t njobs = jobs_.size();

  ASSERT(predData.rows() == npreds*njobs && predData.cols() == nlocs);

  //  map bias coeff to eigen matrix npreds X njobs (read only)
  Eigen::Map<const Eigen::MatrixXd> coeffs(biascoeffs_.data(), npreds, njobs);

  // save ObsBiasTerms (bias_coeff x predictor) for QC
  std::string varname;
  std::vector<double> vec(nlocs);
  for (std::size_t jch = 0; jch < njobs; ++jch) {
    for (std::size_t jpred = 0; jpred < npreds; ++jpred) {
      varname = prednames_[jpred] + "_" + std::to_string(jobs_[jch]);
      Eigen::VectorXd::Map(&vec[0], nlocs) = predData.row(jpred+jch*npreds) * coeffs(jpred, jch);
      odb.put_db(group, varname, vec);
    }
  }

  oops::Log::trace() << "ObsBias::saveObsBiasTerms done." << std::endl;
}

// -----------------------------------------------------------------------------

double ObsBiasLinearCombination::norm() const {
  double zz = 0.0;
  for (unsigned int jj = 0; jj < biascoeffs_.size(); ++jj) {
    zz += biascoeffs_[jj] * biascoeffs_[jj];
  }
  if (biascoeffs_.size() > 0) zz = std::sqrt(zz/this->size());
  return zz;
}

// -----------------------------------------------------------------------------

ObsBiasLinearCombination & ObsBiasLinearCombination::operator+=(const ObsBiasIncrement & dx) {
  for (unsigned int jj = 0; jj < biascoeffs_.size(); ++jj)
    biascoeffs_[jj] += dx[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBiasLinearCombination & ObsBiasLinearCombination::operator=(const ObsBias & rhs) {
  for (unsigned int jj = 0; jj < biascoeffs_.size(); ++jj)
    biascoeffs_[jj] = rhs[jj];
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBiasLinearCombination::print(std::ostream & os) const {
  os << "ObsBiasLinearCombination::print NOT implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

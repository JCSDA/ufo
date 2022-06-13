/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/gnssro/QC/ROobserror.h"

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------

ROobserror::ROobserror(ioda::ObsSpace & obsdb,
                       const eckit::Configuration & config,
                       std::shared_ptr<ioda::ObsDataVector<int> > qc,
                       std::shared_ptr<ioda::ObsDataVector<float> > oberr)
  : FilterBase(obsdb, config, qc, oberr)
{
  oops::Log::trace() << "ROobserror contructor starting" << std::endl;
  const oops::Variables filvar = filtervars_[0].toOopsVariables();;
  oops::Log::trace() << "ROobserror contructor =  "<< filvar << std::endl;
  ufo_roobserror_create_f90(key_, obsdb, config, filvar);
  oops::Log::trace() << "ROobserror contructor key = " << key_ << std::endl;

  // Declare the geovals that are needed by the fortran
  allvars_ += Variable("air_temperature@GeoVaLs");
  allvars_ += Variable("geopotential_height@GeoVaLs");

  // Get the number of horizontal geovals (used by ROPP-2D)
  // Default to 1
  n_horiz = config.getInt("n_horiz", 1);
}

// -----------------------------------------------------------------------------

ROobserror::~ROobserror() {
  oops::Log::trace() << "ROobserror destructor key = " << key_ << std::endl;
  ufo_roobserror_delete_f90(key_);
}

// -----------------------------------------------------------------------------

void ROobserror::applyFilter(const std::vector<bool> & apply,
                             const Variables & filtervars,
                             std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ROobserror using priorFilter" << std::endl;
  // Check that we have valid data to apply the filter to
  if (obsdb_.nlocs() > 0 && data_.getGeoVaLs()->nlocs() > 0) {
    // Get the geovals
    Eigen::ArrayXXf air_temperature = get_geovals("air_temperature@GeoVaLs");
    Eigen::ArrayXXf geopot_height = get_geovals("geopotential_height@GeoVaLs");

    // Call the fortran routines to do the processing
    flags_->save("FortranQC");    // should pass values to fortran properly
    obserr_->save("FortranERR");  // should pass values to fortran properly
    ufo_roobserror_prior_f90(key_,
                             air_temperature.rows(), air_temperature.cols(), air_temperature.data(),
                             geopot_height.rows(), geopot_height.cols(), geopot_height.data());
    flags_->read("FortranQC");    // should get values from fortran properly
    obserr_->read("FortranERR");  // should get values from fortran properly
  }
}


Eigen::ArrayXXf ROobserror::get_geovals(const std::string& var_name) const {
    // Get the geovals
    // Note that ROPP has more geovals than observation locations, and this converts to
    // the correct number (there are n_horiz geovals for every observation)
    size_t nlocs = obsdb_.nlocs() * static_cast<size_t>(n_horiz);
    ASSERT(nlocs == data_.getGeoVaLs()->nlocs());
    size_t nlevs = data_.nlevs(Variable(var_name));
    Eigen::ArrayXXf all_geovals(nlocs, nlevs);
    std::vector<float> single_geoval(nlocs);
    for (int ilev=0; ilev < static_cast<int>(nlevs); ilev++) {
        data_.getGeoVaLs()->getAtLevel(single_geoval, Variable(var_name).variable(), ilev);
        all_geovals.col(ilev) = Eigen::VectorXf::Map(single_geoval.data(), single_geoval.size());
    }
    return all_geovals;
}


// -----------------------------------------------------------------------------

void ROobserror::print(std::ostream & os) const {
  os << "ROobserror::print not yet implemented " << key_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

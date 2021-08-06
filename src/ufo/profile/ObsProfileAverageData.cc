/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/FloatCompare.h"
#include "oops/util/missingValues.h"

#include "ufo/profile/ObsProfileAverageData.h"
#include "ufo/profile/SlantPathLocations.h"
#include "ufo/utils/OperatorUtils.h"  // for getOperatorVariables

namespace ufo {

  ObsProfileAverageData::ObsProfileAverageData(const ioda::ObsSpace & odb,
                                               const eckit::Configuration & config)
    : odb_(odb)
  {
    // Ensure observations have been grouped into profiles.
    if (odb_.obs_group_vars().empty())
      throw eckit::UserError("Group variables configuration is empty", Here());

    // Ensure observations have been sorted by air pressure in descending order.
    if (odb_.obs_sort_var() != "air_pressure")
      throw eckit::UserError("Sort variable must be air_pressure", Here());
    if (odb_.obs_sort_order() != "descending")
      throw eckit::UserError("Profiles must be sorted in descending order", Here());

    // Check the ObsSpace has been extended. If this is not the case
    // then it will not be possible to access profiles in the original and
    // extended sections of the ObsSpace.
    if (!odb_.has("MetaData", "extended_obs_space"))
      throw eckit::UserError("The extended obs space has not been produced", Here());

    // Check the observed pressure is present. This is necessary in order to
    // determine the slant path locations.
    if (!odb_.has("MetaData", "air_pressure"))
      throw eckit::UserError("air_pressure@MetaData not present", Here());

    // Set up configuration options.
    options_.validateAndDeserialize(config);

    // Add model air pressure to the list of variables used in this operator.
    // This GeoVaL is used to determine the slant path locations.
    modelVerticalCoord_ = options_.modelVerticalCoordinate;
    requiredVars_ += oops::Variables({modelVerticalCoord_});

    // Add any simulated variables to the list of variables used in this operator.
    getOperatorVariables(config, odb_.obsvariables(), operatorVars_, operatorVarIndices_);
    requiredVars_ += operatorVars_;

    // If required, set up vectors for OPS comparison.
    if (options_.compareWithOPS.value())
      this->setUpAuxiliaryReferenceVariables();
  }

  const oops::Variables & ObsProfileAverageData::requiredVars() const {
    return requiredVars_;
  }

  const oops::Variables & ObsProfileAverageData::simulatedVars() const {
    return operatorVars_;
  }

  const std::vector<int> & ObsProfileAverageData::operatorVarIndices() const {
    return operatorVarIndices_;
  }

  void ObsProfileAverageData::cacheGeoVaLs(const GeoVaLs & gv) const {
    // Only perform the caching once.
    if (!cachedGeoVaLs_) cachedGeoVaLs_.reset(new GeoVaLs(gv));
  }

  std::vector<std::size_t> ObsProfileAverageData::getSlantPathLocations
  (const std::vector<std::size_t> & locsOriginal,
   const std::vector<std::size_t> & locsExtended) const
  {
    const std::vector<std::size_t> slant_path_location =
      ufo::getSlantPathLocations(odb_,
                                 *cachedGeoVaLs_,
                                 locsOriginal,
                                 modelVerticalCoord_,
                                 options_.numIntersectionIterations.value() - 1);

    // If required, compare slant path locations and slant pressure with OPS output.
    if (options_.compareWithOPS.value()) {
      // Vector of slanted pressures, used for comparisons with OPS.
      std::vector<float> slant_pressure;
      // Number of levels for model pressure.
      const std::size_t nlevs_p = cachedGeoVaLs_->nlevs(modelVerticalCoord_);
      // Vector used to store different pressure GeoVaLs.
      std::vector <float> pressure_gv(nlevs_p);
      for (std::size_t mlev = 0; mlev < nlevs_p; ++mlev) {
        cachedGeoVaLs_->getAtLocation(pressure_gv, modelVerticalCoord_, slant_path_location[mlev]);
        slant_pressure.push_back(pressure_gv[mlev]);
      }
      this->compareAuxiliaryReferenceVariables(locsExtended,
                                               slant_path_location,
                                               slant_pressure);
    }

    return slant_path_location;
  }

  void ObsProfileAverageData::setUpAuxiliaryReferenceVariables() {
    if (!(odb_.has("MetOfficeHofX", "slant_path_location") &&
          odb_.has("MetOfficeHofX", "slant_pressure")))
      throw eckit::UserError("At least one reference variable is not present", Here());
    // Get reference values of the slant path locations and pressures.
    slant_path_location_ref_.resize(odb_.nlocs());
    slant_pressure_ref_.resize(odb_.nlocs());
    odb_.get_db("MetOfficeHofX", "slant_path_location", slant_path_location_ref_);
    odb_.get_db("MetOfficeHofX", "slant_pressure", slant_pressure_ref_);
  }

  void ObsProfileAverageData::compareAuxiliaryReferenceVariables
  (const std::vector<std::size_t> & locsExtended,
   const std::vector<std::size_t> & slant_path_location,
   const std::vector<float> & slant_pressure) const {
    std::vector<int> slant_path_location_ref_profile;
    std::vector<float> slant_pressure_ref_profile;
    for (const std::size_t loc : locsExtended) {
      slant_path_location_ref_profile.push_back(slant_path_location_ref_[loc]);
      slant_pressure_ref_profile.push_back(slant_pressure_ref_[loc]);
    }
    std::stringstream errmsg;
    for (std::size_t mlev = 0; mlev < locsExtended.size(); ++mlev) {
      if (slant_path_location[mlev] != slant_path_location_ref_profile[mlev]) {
        errmsg << "Mismatch for slant_path_location, level = " << mlev
               << " (this code, OPS): "
               << slant_path_location[mlev] << ", "
               << slant_path_location_ref_profile[mlev];
        throw eckit::BadValue(errmsg.str(), Here());
      }
      if (!oops::is_close_relative(slant_pressure[mlev],
                                   slant_pressure_ref_profile[mlev],
                                   1e-5f)) {
        errmsg << "Mismatch for slant_pressure, level = " << mlev
               << " (this code, OPS): "
               << slant_pressure[mlev] << ", "
               << slant_pressure_ref_profile[mlev];
        throw eckit::BadValue(errmsg.str(), Here());
      }
    }
  }

  void ObsProfileAverageData::print(std::ostream & os) const {
    os << "ObsProfileAverage operator" << std::endl;
    os << "config = " << options_ << std::endl;
  }
}  // namespace ufo

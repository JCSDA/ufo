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
                                               const ObsProfileAverageParameters & parameters,
                                               const VariableNameMap & nameMap)
    : odb_(odb),
      options_(parameters),
      modelVerticalCoord_(oops::Variable{options_.modelVerticalCoordinate})
  {
    // Ensure observations have been grouped into profiles.
    if (odb_.obs_group_vars().empty())
      throw eckit::UserError("Group variables configuration is empty", Here());

    // Ensure observations have been sorted by air pressure in descending order.
    if (options_.requireDescendingPressureSort) {
      if (odb_.obs_sort_var() != parameters.pressureCoord.value())
        throw eckit::UserError(std::string("Sort variable must be ") +
                               parameters.pressureCoord.value(),
                               Here());
      if (odb_.obs_sort_order() != "descending")
        throw eckit::UserError("Profiles must be sorted in descending order", Here());
    }

    // Check the ObsSpace has been extended. If this is not the case
    // then it will not be possible to access profiles in the original and
    // extended sections of the ObsSpace.
    if (!odb_.has("MetaData", "extendedObsSpace"))
      throw eckit::UserError("The extended obs space has not been produced", Here());

    // Add model air pressure to the list of variables used in this operator.
    // This GeoVaL is used to determine the slant path locations.
    requiredVars_.push_back(modelVerticalCoord_);
    geovalsObsSameDir_ = options_.geovalsObsSameDir;

    // Add any simulated variables to the list of variables used in this operator.
    getOperatorVariables(parameters.variables.value(), odb_.assimvariables(),
                         operatorVars_, operatorVarIndices_);
    requiredVars_ += nameMap.convertName(operatorVars_);

    // If required, set up vectors for OPS comparison.
    if (options_.compareWithOPS.value())
      this->setUpAuxiliaryReferenceVariables();
  }

  const oops::Variables & ObsProfileAverageData::requiredVars() const {
    return requiredVars_;
  }

  const oops::ObsVariables & ObsProfileAverageData::simulatedVars() const {
    return operatorVars_;
  }

  const std::vector<int> & ObsProfileAverageData::operatorVarIndices() const {
    return operatorVarIndices_;
  }

  void ObsProfileAverageData::cacheGeoVaLs(const GeoVaLs & gv) const {
    // Only perform the caching once.
    if (!cachedGeoVaLs_)
      cachedGeoVaLs_.reset(new GeoVaLs(gv));
  }

  std::vector<std::size_t> ObsProfileAverageData::getSlantPathLocations
  (const std::vector<std::size_t> & locsOriginal,
   const std::vector<std::size_t> & locsExtended) const
  {
    const std::vector<std::size_t> slant_path_location =
      ufo::getSlantPathLocations(odb_,
                                 *cachedGeoVaLs_,
                                 locsOriginal,
                                 options_.pressureGroup.value() +
                                 std::string("/") +
                                 options_.pressureCoord.value(),
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
        slant_pressure.push_back(pressure_gv[nlevs_p - 1 - mlev]);
      }
      this->compareAuxiliaryReferenceVariables(locsExtended,
                                               slant_path_location,
                                               slant_pressure);
    }

    return slant_path_location;
  }

  void ObsProfileAverageData::setUpAuxiliaryReferenceVariables() {
    if (!(odb_.has("MetOfficeHofX", "slantPathLocation") &&
          odb_.has("MetOfficeHofX", "slantPressure")))
      throw eckit::UserError("At least one reference variable is not present", Here());
    // Get reference values of the slant path locations and pressures.
    slant_path_location_ref_.resize(odb_.nlocs());
    slant_pressure_ref_.resize(odb_.nlocs());
    odb_.get_db("MetOfficeHofX", "slantPathLocation", slant_path_location_ref_);
    odb_.get_db("MetOfficeHofX", "slantPressure", slant_pressure_ref_);
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
        errmsg << "Mismatch for slantPathLocation, level = " << mlev
               << " (this code, OPS): "
               << slant_path_location[mlev] << ", "
               << slant_path_location_ref_profile[mlev];
        throw eckit::BadValue(errmsg.str(), Here());
      }
      if (!oops::is_close_relative(slant_pressure[mlev],
                                   slant_pressure_ref_profile[mlev],
                                   1e-5f)) {
        errmsg << "Mismatch for slantPressure, level = " << mlev
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

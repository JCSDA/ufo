/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PROFILE_OBSPROFILEAVERAGEDATA_H_
#define UFO_PROFILE_OBSPROFILEAVERAGEDATA_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/profile/ObsProfileAverageParameters.h"
#include "ufo/utils/VariableNameMap.h"

namespace ufo {

  /// \brief Data handling class for the ProfileAverage observation operator and TL/AD code.
  class ObsProfileAverageData {
   public:
    ObsProfileAverageData(const ioda::ObsSpace & odb,
                          const ObsProfileAverageParameters & parameters,
                          const VariableNameMap & nameMap = VariableNameMap(boost::none));

    /// Return required variables for the operator.
    const oops::Variables & requiredVars() const;

    /// Return simulated variables for the operator.
    const oops::ObsVariables & simulatedVars() const;

    /// Return operator variable indices for the operator.
    const std::vector<int> & operatorVarIndices() const;

    /// Cache the initial values of the GeoVaLs.
    void cacheGeoVaLs(const GeoVaLs & gv) const;

    /// Get slant path locations. This determines, for each model level, the location that
    /// corresponds to the intersection of the observed profile with that level.
    std::vector<std::size_t> getSlantPathLocations
      (const std::vector<std::size_t> & locsOriginal,
       const std::vector<std::size_t> & locsExtended) const;

    /// Print operator configuration options.
    void print(std::ostream & os) const;

    /// Get the name of the model vertical coordinate.
    oops::Variable getModelVerticalCoord() const {return modelVerticalCoord_;}

    /// Get the switch indicating if geovals are in the same direction as observations.
    /// (Beware, default 'false'!)
    bool geovalsObsSameDir() const {return geovalsObsSameDir_;}

   private:
    /// Set up auxiliary reference variables that are used for comparison with OPS.
    /// These reference variables are called MetOfficeHofX/slant_path_location and
    /// MetOfficeHofX/slant_pressure. If a comparison with OPS is to be performed
    /// then these variables must be present in the input data set.
    void setUpAuxiliaryReferenceVariables();

    /// Compare auxiliary reference variables with those obtained in OPS.
    void compareAuxiliaryReferenceVariables(const std::vector<std::size_t> & locsExtended,
                                            const std::vector<std::size_t> & slant_path_location,
                                            const std::vector<float> & slant_pressure) const;

   private:
    /// ObsSpace.
    const ioda::ObsSpace & odb_;

    /// Options for this observation operator.
    ObsProfileAverageParameters options_;

    /// Name of model vertical coordinate.
    oops::Variable modelVerticalCoord_;

    /// Required variables.
    oops::Variables requiredVars_;

    /// Operator variables.
    oops::ObsVariables operatorVars_;

    /// Indices of operator variables.
    std::vector<int> operatorVarIndices_;

    /// Cached GeoVaLs.
    mutable std::unique_ptr<GeoVaLs> cachedGeoVaLs_;

    /// Reference values of slant path locations.
    std::vector<int> slant_path_location_ref_;

    /// Reference values of slant path pressures.
    std::vector<float> slant_pressure_ref_;

    /// Whether geovals and observations are in the same direction. (Beware, default 'false'!)
    bool geovalsObsSameDir_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_OBSPROFILEAVERAGEDATA_H_

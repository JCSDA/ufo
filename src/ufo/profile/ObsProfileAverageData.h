/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PROFILE_OBSPROFILEAVERAGEDATA_H_
#define UFO_PROFILE_OBSPROFILEAVERAGEDATA_H_

#include <ostream>
#include <vector>

#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/profile/ObsProfileAverageParameters.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

  /// \brief Data handling class for the ProfileAverage observation operator and TL/AD code.
  class ObsProfileAverageData {
   public:
    ObsProfileAverageData(const ioda::ObsSpace & odb,
                          const eckit::Configuration & config);

    /// Return required variables for the operator.
    const oops::Variables & requiredVars() const;

    /// Get slant path locations. This determines, for each model level, the location that
    /// corresponds to the intersection of the observed profile with that level.
    std::vector<std::size_t> getSlantPathLocations(const std::vector<std::size_t> & locsOriginal,
                                                   const std::vector<std::size_t> & locsExtended,
                                                   const GeoVaLs & gv) const;

    /// Print operator configuration options.
    void print(std::ostream & os) const;

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

    /// Required variables.
    oops::Variables requiredVars_;

    /// Reference values of slant path locations.
    std::vector<int> slant_path_location_ref_;

    /// Reference values of slant path pressures.
    std::vector<float> slant_pressure_ref_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_OBSPROFILEAVERAGEDATA_H_

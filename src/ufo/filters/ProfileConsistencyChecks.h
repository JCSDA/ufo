/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PROFILECONSISTENCYCHECKS_H_
#define UFO_FILTERS_PROFILECONSISTENCYCHECKS_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

#include "ufo/utils/Constants.h"
#include "ufo/utils/Flags.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {
  class ProfileConsistencyCheckParameters;
}

namespace ufo {

  /// \brief Check profile temperature/height data for internal consistency.
  ///
  /// See OSDP 5 section 3. The checks made are as follows:
  ///    -# Basic checks of pressure
  ///    -# Check for different temperatures at the same pressure
  ///    -# Temperature sign change test using model background
  ///    -# Check for superadiabatic layers more than 50 hPa above surface
  ///    -# Check of standard level temperatures against values interpolated from
  ///       significant levels
  ///    -# Check for hydrostatic consistency of standard levels (+ surface)
  ///       -# calculation of thickness residuals
  ///       -# decision making algorithm used where there are large residuals
  ///    -# Junk whole report if > 8 errors
  ///
  /// The sign, standard level and hydrostatic checks are largely based on methods described in the
  /// WMO Guide on the Global Data-Procesing System (1993).
  /// In general the tolerances have been relaxed slightly so that fewer
  /// 'marginal' cases are flagged.
  /// In the hydrostatic check the effect of moisture on the thickness is
  /// neglected and no special account is taken of the tropopause.
  /// Corrections are suggested for apparently simple temperature or height
  /// errors.  If a height correction is within 10 m of a multiple of 100 m then
  /// the rounded correction is applied, temperature corrections are not applied.
  /// In some cases interpolation flags are switched off if the hydrostatic
  /// check is fairly clear-cut.
  ///
  /// Inputs:
  /// - Arrays containing observed and background values
  ///   - Should be in order of increasing height/decreasing pressure
  ///   - interleaved standard and significant levels
  ///     (if input is standard levels only then superadiabatic and
  ///     hydrostatic checks will be performed)
  ///
  /// Outputs:
  /// - QC flags
  /// - Some corrections to reported values (not applied by default):
  ///   - corrections to sign of temperature
  ///   - corrections to height (multiples of 100 m)

  class ProfileConsistencyChecks : public FilterBase,
    private util::ObjectCounter<ProfileConsistencyChecks> {
   public:
      static const std::string classname() {return "ufo::ProfileConsistencyChecks";}

      ProfileConsistencyChecks(ioda::ObsSpace &, const eckit::Configuration &,
                               boost::shared_ptr<ioda::ObsDataVector<int> >,
                               boost::shared_ptr<ioda::ObsDataVector<float> >);
      ~ProfileConsistencyChecks();

      /// Return the number of mismatches between values produced by the checking routines
      /// and the equivalents produced in the OPS code.
      /// The values checked are: the QC flags for each observation,
      /// intermediate values used in various calculations for each observation,
      /// and error counters for each profile.
      std::vector <int> getMismatches() const {return nMismatches_;}

   private:
      void print(std::ostream &) const override;
      void applyFilter(const std::vector<bool> &, const Variables &,
                       std::vector<std::vector<bool>> &) const override;
      int qcFlag() const override {return QCflags::profile;}

      /// Configurable options
      std::unique_ptr <ProfileConsistencyCheckParameters> options_;

      /// Number of mismatches between values produced in this code
      /// and their OPS equivalents (used for validation)
      mutable std::vector <int> nMismatches_;
  };

}  // namespace ufo

#endif  // UFO_FILTERS_PROFILECONSISTENCYCHECKS_H_

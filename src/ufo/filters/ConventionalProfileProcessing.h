/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_CONVENTIONALPROFILEPROCESSING_H_
#define UFO_FILTERS_CONVENTIONALPROFILEPROCESSING_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"

#include "ufo/filters/ConventionalProfileProcessingParameters.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

#include "ufo/profile/EntireSampleDataHandler.h"
#include "ufo/profile/ProfileChecker.h"
#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/VariableNames.h"

#include "ufo/utils/Constants.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

  /// \brief Conventional profile processing.
  ///
  /// This filter applies a variety of QC checks to profile data.
  /// The filter also averages profiles onto model levels prior to their use in data assimilation.
  ///
  /// The temperature consistency checks available are as follows:
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
  //
  /// There are also consistency checks for the interpolated wind speed and the relative humidity.
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

  class ConventionalProfileProcessing : public FilterBase,
    private util::ObjectCounter<ConventionalProfileProcessing> {
   public:
      typedef ConventionalProfileProcessingParameters Parameters_;

      static const std::string classname() {return "ufo::ConventionalProfileProcessing";}

      ConventionalProfileProcessing(ioda::ObsSpace &, const Parameters_ &,
                                    std::shared_ptr<ioda::ObsDataVector<int> >,
                                    std::shared_ptr<ioda::ObsDataVector<float> >);
      ~ConventionalProfileProcessing();

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

      /// Run checks on individual profiles sequentially.
      void individualProfileChecks(ProfileDataHandler &profileDataHandler,
                                   ProfileCheckValidator &profileCheckValidator,
                                   ProfileChecker &profileChecker,
                                   const CheckSubgroup &subGroupChecks) const;

      /// Run checks that use all of the profiles at once.
      void entireSampleChecks(ProfileDataHandler &profileDataHandler,
                              ProfileCheckValidator &profileCheckValidator,
                              ProfileChecker &profileChecker,
                              const CheckSubgroup &subGroupChecks) const;

      int qcFlag() const override {return QCflags::profile;}

      /// Configurable options
      ConventionalProfileProcessingParameters options_;

      /// Number of mismatches between values produced in this code
      /// and their OPS equivalents (used for validation)
      mutable std::vector <int> nMismatches_;
  };

}  // namespace ufo

#endif  // UFO_FILTERS_CONVENTIONALPROFILEPROCESSING_H_

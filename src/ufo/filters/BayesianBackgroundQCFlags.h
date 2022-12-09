/*
 * (C) Crown Copyright 2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_BAYESIANBACKGROUNDQCFLAGS_H_
#define UFO_FILTERS_BAYESIANBACKGROUNDQCFLAGS_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/ProbabilityOfGrossErrorParameters.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the BayesianBackgroundQCFlags filter.
class BayesianBackgroundQCFlagsParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(BayesianBackgroundQCFlagsParameters, FilterParametersBase)

 public:
  /// This parameter is a collection of PGE variable name substitutions.
  /// The PGE of the second variable in each pair is used to set the QC flags
  /// of the first variable; this happens for (e.g.) wind u and v components.
  /// todo(UKMO): might need to add the variables known as
  /// u10, v10 and uSuperob, vSuperob in OPS.
  oops::Parameter<std::map<std::string, std::string>> PGEsubstituteNames
    {"PGE variable name substitutions", {{"windNorthward", "windEastward"}}, this};

  /// Parameters related to PGE calculations. The value of \p PGECrit
  /// is obtained from here.
  ProbabilityOfGrossErrorParameters PGEParameters{this};
};

/// BayesianBackgroundQCFlags: apply QC flags based on values of probability of
/// gross error (PGE). If the PGE is larger than the threshold \p PGEcrit then
/// the observation is rejected.
/// If the Bayesian background or buddy checks were applied,
/// use the PGE that was obtained from those checks.
/// Sometimes the PGE of one variable is used to set the QC flags of another;
/// this happens for (e.g.) wind u and v components.
///
/// todo(UKMO): deal with Pstar/Pmsl and u10AmbWind/v10AmbWind (as they are known in OPS).
/// These are each treated slightly differently in the OPS code.
class BayesianBackgroundQCFlags : public FilterBase,
  private util::ObjectCounter<BayesianBackgroundQCFlags> {
 public:
  typedef BayesianBackgroundQCFlagsParameters Parameters_;

  static const std::string classname() {return "ufo::BayesianBackgroundQCFlags";}

  BayesianBackgroundQCFlags(ioda::ObsSpace &, const Parameters_ &,
                            std::shared_ptr<ioda::ObsDataVector<int> >,
                            std::shared_ptr<ioda::ObsDataVector<float> >);
  ~BayesianBackgroundQCFlags();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::bayesianQC;}

  /// Get the name of the variable whose PGE is tested in order to
  /// set the QC flags for the variable \p varname.
  /// By default, \p varname is returned by this routine; any substitutions
  /// are listed in the \p PGEsubstituteNames parameter.
  std::string getPGEsubstituteName(const std::string& varname) const;

  /// Set flags for the variable \p varname given the \p apply vector.
  /// Set both integer bitmap flags and an overall filter flag (bayesianQC).
  void setFlags(const std::string& varname,
                const std::vector<bool>& apply,
                std::vector<bool>& flagged) const;

  /// Parameters used in this filter.
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_BAYESIANBACKGROUNDQCFLAGS_H_

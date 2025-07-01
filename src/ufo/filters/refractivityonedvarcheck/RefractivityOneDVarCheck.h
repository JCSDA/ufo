/*
 * (C) Copyright 2025 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_REFRACTIVITYONEDVARCHECK_REFRACTIVITYONEDVARCHECK_H_
#define UFO_FILTERS_REFRACTIVITYONEDVARCHECK_REFRACTIVITYONEDVARCHECK_H_

#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/refractivityonedvarcheck/RefractivityOneDVarCheck.interface.h"
#include "ufo/filters/refractivityonedvarcheck/RefractivityOneDVarCheckParameters.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

//! \brief RefractivityOneDVarCheck
//!
//! \details RefractivityOneDVarCheck performs a 1D-Var minimization for satellite using the Met
//! Office's GNSS-RO forward operator.  If a profile does not converge then all observations
//! in this profile are flagged.  The code is based on the Met Office 1D-Var scheme and thus
//! is predominently in Fortran.
//!
//! \author Met Office
//!
//! \date 13/05/2025
//!

class RefractivityOneDVarCheck : public FilterBase,
                     private util::ObjectCounter<RefractivityOneDVarCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef RefractivityOneDVarCheckParameters Parameters_;

  static const std::string classname() {return "ufo::RefractivityOneDVarCheck";}

  RefractivityOneDVarCheck(
    ioda::ObsSpace &,
    const Parameters_ &,
    std::shared_ptr<ioda::ObsDataVector<int> >,
    std::shared_ptr<ioda::ObsDataVector<float> >);
  ~RefractivityOneDVarCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &,
                   const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::onedvar;}

  F90onedvarcheck key_;
  RefractivityOneDVarCheckParameters parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_REFRACTIVITYONEDVARCHECK_REFRACTIVITYONEDVARCHECK_H_

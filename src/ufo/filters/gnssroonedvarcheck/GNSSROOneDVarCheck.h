/*
 * (C) Copyright 2017-2020 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_GNSSROONEDVARCHECK_GNSSROONEDVARCHECK_H_
#define UFO_FILTERS_GNSSROONEDVARCHECK_GNSSROONEDVARCHECK_H_

#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/gnssroonedvarcheck/GNSSROOneDVarCheck.interface.h"
#include "ufo/filters/gnssroonedvarcheck/GNSSROOneDVarCheckParameters.h"
#include "ufo/filters/QCflags.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

//! \brief GNSSROOneDVarCheck
//!
//! \details GNSSROOneDVarCheck performs a 1D-Var minimization for satellite using the Met
//! Office's GNSS-RO forward operator.  If a profile does not converge the all observations
//! in this profile are flagged.  The code is based on the Met Office 1D-Var scheme and thus
//! is predominently in Fortran.
//!
//! \author Met Office
//!
//! \date 16/11/2020
//!

class GNSSROOneDVarCheck : public FilterBase,
                     private util::ObjectCounter<GNSSROOneDVarCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef GNSSROOneDVarCheckParameters Parameters_;

  static const std::string classname() {return "ufo::GNSSROOneDVarCheck";}

  GNSSROOneDVarCheck(ioda::ObsSpace &,
                     const Parameters_ &,
                     std::shared_ptr<ioda::ObsDataVector<int> >,
                     std::shared_ptr<ioda::ObsDataVector<float> >);
  ~GNSSROOneDVarCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &,
                   const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::onedvar;}

  F90onedvarcheck key_;
  GNSSROOneDVarCheckParameters parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_GNSSROONEDVARCHECK_GNSSROONEDVARCHECK_H_

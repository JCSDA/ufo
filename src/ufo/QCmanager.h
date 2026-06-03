/*
 * (C) Copyright 2026 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_QCMANAGER_H_
#define UFO_QCMANAGER_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsVector;
  class ObsSpace;
}

namespace ufo {

/// \brief QCManager encodes basic QC logic for missing data and passivated observations.
/// It runs whenever one or more filters are requested, so that the basic logic doesn't have
/// to be duplicated in every filter set.
///
/// The preSetQc() function sets the QC flag to `missing` at all
/// locations with missing or invalid obs values.
/// The postSetQc() function sets the QC flag to `Hfailed` if it was previously set to
/// `pass`, but the obs operator failed to produce a valid HofX value.
/// The finalSetQc() function sets the QC flag to `processed` if
/// the observation is not to be assimilated and also sets the QC flag to `missing`
/// if it is not already rejected but has missing error estimates.
class QCmanager : public util::Printable,
                  private util::ObjectCounter<QCmanager> {
 public:
  static const std::string classname() {return "ufo::QCmanager";}
  QCmanager(ioda::ObsSpace &,
            ioda::ObsDataVector<int> &,
            ioda::ObsDataVector<float> &);
  ~QCmanager();

  void preSetQc();
  // No QC logic that depends on GeoVaLs is needed at this time,
  // but priorSetQc function is provided for future use if needed.
  void priorSetQc() {}
  void postSetQc(const ioda::ObsVector &);
  void finalSetQc();

 private:
  void print(std::ostream &) const;

  ioda::ObsSpace & obsdb_;
  ioda::ObsDataVector<int> & flags_;
  ioda::ObsDataVector<float> & obserr_;
};

}  // namespace ufo

#endif  // UFO_QCMANAGER_H_

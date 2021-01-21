/*
 * (C) Copyright 2017-2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSPROCESSORBASE_H_
#define UFO_FILTERS_OBSPROCESSORBASE_H_

#include <memory>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Base class for UFO observation processors (including QC filters).
///
/// Observation processors only need to implement the constructor and the doFilter method;
/// the base class takes care of applying the processor at the pre, prior or post stage.

class ObsProcessorBase : public util::Printable {
 public:
  ObsProcessorBase(ioda::ObsSpace &, bool deferToPost,
                   std::shared_ptr<ioda::ObsDataVector<int> >,
                   std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsProcessorBase();

  void preProcess();
  void priorFilter(const GeoVaLs &);
  void postFilter(const ioda::ObsVector &, const ObsDiagnostics &);

  oops::Variables requiredVars() const {
    return allvars_.allFromGroup("GeoVaLs").toOopsVariables();}
  oops::Variables requiredHdiagnostics() const {
    return allvars_.allFromGroup("ObsDiag").toOopsVariables();}

 protected:
  ioda::ObsSpace & obsdb_;
  std::shared_ptr<ioda::ObsDataVector<int>> flags_;
  std::shared_ptr<ioda::ObsDataVector<float>> obserr_;
  ufo::Variables allvars_;
  ObsFilterData data_;

 private:
  virtual void doFilter() const = 0;

  bool prior_;
  bool post_;

  // Variables extracted from the filter parameters.
  bool deferToPost_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSPROCESSORBASE_H_

/*
 * (C) Copyright 2017-2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSPROCESSORBASE_H_
#define UFO_FILTERS_OBSPROCESSORBASE_H_

#include <memory>

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"
#include "ufo/ObsTraits.h"

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

class ObsProcessorBase : public oops::interface::ObsFilterBase<ObsTraits> {
 public:
  ObsProcessorBase(ioda::ObsSpace &, bool deferToPost,
                   std::shared_ptr<ioda::ObsDataVector<int> >,
                   std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsProcessorBase();

  void preProcess() override;
  void priorFilter(const GeoVaLs &) override;
  void postFilter(const GeoVaLs &,
                  const ioda::ObsVector &,
                  const ioda::ObsVector &,
                  const ObsDiagnostics &) override;
  void checkFilterData(const oops::FilterStage filterStage) override;

  oops::Variables requiredVars() const override {
    return allvars_.allFromGroup("GeoVaLs").toOopsVariables();}
  oops::ObsVariables requiredHdiagnostics() const override {
    return allvars_.allFromGroup("ObsDiag").toOopsObsVariables();}

 protected:
  ioda::ObsSpace & obsdb_;
  std::shared_ptr<ioda::ObsDataVector<int>> flags_;
  std::shared_ptr<ioda::ObsDataVector<float>> obserr_;
  ufo::Variables allvars_;
  ObsFilterData data_;
  bool prior_;
  bool post_;

 private:
  virtual void doFilter() = 0;

  // Variables extracted from the filter parameters.
  bool deferToPost_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSPROCESSORBASE_H_

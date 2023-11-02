/*
 * (C) Crown Copyright 2023 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_METOFFICEDUPLICATECHECK_H_
#define UFO_FILTERS_METOFFICEDUPLICATECHECK_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"

namespace ufo {

class MetOfficeDuplicateCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(MetOfficeDuplicateCheckParameters, FilterParametersBase)

 public:
  oops::RequiredParameter <std::string> priorityName
    {"priority name",
     "Name of the priority variable in the MetaData group. "
     "Observations with higher priority are preferentially retained.",
     this};
  oops::RequiredParameter <float> latBandWidth
    {"latitude band width",
     "Width of latitude bands used in the duplicate check.",
     this};
  oops::RequiredParameter <float> latBinHalfWidth
    {"latitude bin half-width",
     "Half-width (degrees) of latitude bin used to define duplicate pairs of observations.",
     this};
  oops::RequiredParameter <float> lonBinHalfWidth
    {"longitude bin half-width",
     "Half-width (degrees) of longitude bin used to define duplicate pairs of observations.",
     this};
  oops::RequiredParameter <int64_t> timeBinHalfWidth
    {"time bin half-width",
     "Half-width (s) of time bin used to define duplicate pairs of observations.",
     this};
  oops::RequiredParameter <float> pBinHalfWidth
    {"pressure bin half-width",
     "Half-width (Pa) of pressure bin used to define duplicate pairs of observations.",
     this};
  oops::OptionalParameter <std::string> categoryVariableName
    {"category variable name",
     "Name of a category variable in the MetaData group. "
     "If defined, whenever an observation O is thinned, all other observations with the same value "
     "of the category variable as O are also thinned.",
     this};
  oops::Parameter <bool> preSortByRecordID
    {"pre-sort by record ID",
     "If true, sort the valid observation IDs before the filter by recordID to guarantee the check "
     "is independent of the number of MPI processes. The default value is false.",
     false,
     this};
  oops::Parameter <Variable> verticalCoordinateVariableName
    {"vertical coordinate variable",
     "If set, use this variable as the vertical coordinate instead of the pressure."
     " The default value is MetaData/pressure.",
     Variable("MetaData/pressure"),
     this};
};

}  // namespace ufo


namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

class ObsAccessor;

/// \brief Thins observations according to the Met Office Duplicate Check algorithm.
/// This algorithm divides the globe into latitude bands of width `latitude band width`.
/// Observations in each latitude band are sorted on longitude from lowest to highest;
/// for simplicity, the discontinuity at the edges is not considered.
/// Starting at the band nearest to the North Pole, and moving downwards, the check
/// determines whether any pairs of observations are colocated inside a volume defined by the
/// parameters `latitude bin half-width`, `longitude bin half-width`, `time bin half-width`
/// and `pressure bin half-width`. If such a pair is found, the observation with the highest
/// priority  variable (which is in the MetaData group with name equal to `priority name`)
/// is retained.
/// In the event of a tie the observation with the lower value of longitude is retained.
/// For a given latitude band the algorithm searches in that band and the adjacent one.
/// his filter is designed to reproduce the Met Office OPS code so the sorting by longitude
/// is performed using the Met Office sorting algorithm.

class MetOfficeDuplicateCheck : public FilterBase,
  private util::ObjectCounter<MetOfficeDuplicateCheck> {
 public:
  typedef MetOfficeDuplicateCheckParameters Parameters_;

  static const std::string classname() {return "ufo::MetOfficeDuplicateCheck";}

  MetOfficeDuplicateCheck(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                      std::shared_ptr<ioda::ObsDataVector<int> > flags,
                      std::shared_ptr<ioda::ObsDataVector<float> > obserr);

  ~MetOfficeDuplicateCheck() override;

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::thinned;}
  ObsAccessor createObsAccessor() const;

 private:
  Parameters_ options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_METOFFICEDUPLICATECHECK_H_

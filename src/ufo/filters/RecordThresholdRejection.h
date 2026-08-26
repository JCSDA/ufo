/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_RECORDTHRESHOLDREJECTION_H_
#define UFO_FILTERS_RECORDTHRESHOLDREJECTION_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/FilterBase.h"
#include "ufo/filters/FilterParametersBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

enum class RecordThresholdRejectionType {
  LT, LTE, GT, GTE
};

enum class RecordThresholdRejectionDataOrder {
  Ascending, Descending
};

struct RecordThresholdRejectionTypeParameterTraitsHelper {
  typedef RecordThresholdRejectionType EnumType;
  static constexpr char enumTypeName[] = "RecordThresholdRejectionType";
  static constexpr util::NamedEnumerator<RecordThresholdRejectionType> namedValues[] = {
    { RecordThresholdRejectionType::LT, "less than" },
    { RecordThresholdRejectionType::LTE, "less than or equal to" },
    { RecordThresholdRejectionType::GT, "greater than" },
    { RecordThresholdRejectionType::GTE, "greater than or equal to" }
  };
};

struct RecordThresholdRejectionDataOrderParameterTraitsHelper {
  typedef RecordThresholdRejectionDataOrder EnumType;
  static constexpr char enumTypeName[] = "RecordThresholdRejectionDataOrder";
  static constexpr util::NamedEnumerator<RecordThresholdRejectionDataOrder> namedValues[] = {
    { RecordThresholdRejectionDataOrder::Ascending, "ascending" },
    { RecordThresholdRejectionDataOrder::Descending, "descending" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::RecordThresholdRejectionType> :
    public EnumParameterTraits<ufo::RecordThresholdRejectionTypeParameterTraitsHelper>
{};

template <>
struct ParameterTraits<ufo::RecordThresholdRejectionDataOrder> :
    public EnumParameterTraits<ufo::RecordThresholdRejectionDataOrderParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// \brief Parameters controlling the operation of the RecordThresholdRejection filter.
class RecordThresholdRejectionParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(RecordThresholdRejectionParameters, FilterParametersBase)

 public:
  /// Name of the threshold value, can be a real number or the name of a variable
  oops::RequiredParameter<std::string> threshold_value{"threshold value", this};
  /// Name of the threshold variable
  oops::RequiredParameter<std::string> threshold_variable{"threshold variable", this};
  /// Rejection type: less than, less than or equal to, greater than, greater than or equal to
  oops::RequiredParameter<RecordThresholdRejectionType> rejection_type{"rejection type", this};

  oops::Parameter<RecordThresholdRejectionDataOrder> data_order{"data order",
    RecordThresholdRejectionDataOrder::Ascending, this};
};

// -----------------------------------------------------------------------------
/*! \brief Filter to perform actions for all entries within a record after a
 * certain value is reached (in ascending or descending data ordering).
 *
 * \details The specified record variable is compared to the specified threshold
 * value for each entry within the record. If the rejection type is less than 
 * (greater than), then the filter flags all entries before (or after,
 * depending on the data order value) this entry within the record.
 *
 * Example for relative humidity:
 * \code{.unparsed}
 *  obs filters:
 *  - filter: Record Threshold Rejection
 *    filter variables: [relativeHumidity]
 *    threshold value: MetaData/AvgRH_InstrTThresholds
 *    # alternatively, could put a scalar value (e.g. -40.0) here to apply the same threshold for all instruments
 *    threshold variable: DerivedObsValue/airTemperature
 *    rejection type: less than # could be 'less than', 'less than or equal to', 'greater than', or 'greater than or equal to'
 *    data order: ascending  # optional, default is ascending; can be ascending or descending
 *    action:
 *      name: set
 *      flag: PermanentStationRejection
 *    where:
 *    - variable: MetaData/extendedObsSpace
 *      is_in: 1
 * \endcode
 *
 * \author S.Cameron (Met Office)
 *
 * \date 11/06/2026: Created
 */
class RecordThresholdRejection : public FilterBase,
                       private util::ObjectCounter<RecordThresholdRejection> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef RecordThresholdRejectionParameters Parameters_;

  static const std::string classname() {return "ufo::RecordThresholdRejection";}

  RecordThresholdRejection(ioda::ObsSpace &, const Parameters_ &,
                   ioda::ObsDataVector<int> &,
                   ioda::ObsDataVector<float> &);
  ~RecordThresholdRejection();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::recordthreshold;}
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_RECORDTHRESHOLDREJECTION_H_

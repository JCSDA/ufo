/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_CAL_QNHTOQFEPRESSURE_H_
#define UFO_VARIABLETRANSFORMS_CAL_QNHTOQFEPRESSURE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

/*!
* \brief Invert conversion of QFE (Query: Field Elevation; mean sea level pressure corrected
* for temperature, adjusted for a specific altitude) to QNH (Query: Nautical Height;
* The pressure measured at mean sea level). Also adds a bias correction to reports that have
* been rounded down to the nearest whole hPa. Applicable to Metars only.
*
* \details  For Metars invert conversion of QFE to QNH following the equations given in ICAO (2011)
* Manual on Automatic Meterological Observing Systems at Aerodromes Doc 9837 AN/454 - chapter 9
* on pressure. Also adds a bias correction of 0.5hPa to reports that have been rounded
* down to the nearest whole hPa; the observation error variance is increased to take account of
* this. Note: QNH is stored as Pmsl and QFE as Pstation
*
*/

class Cal_QNHtoQFEpressure : public TransformBase {
 public:
  Cal_QNHtoQFEpressure(const GenericVariableTransformParameters &options,
                       const ObsFilterData &data,
                       const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                       const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
  // Run variable conversion
  void runTransform(const std::vector<bool> &apply) override;
};
}  // namespace ufo
#endif  // UFO_VARIABLETRANSFORMS_CAL_QNHTOQFEPRESSURE_H_

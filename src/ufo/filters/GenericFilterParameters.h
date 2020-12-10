/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_GENERICFILTERPARAMETERS_H_
#define UFO_FILTERS_GENERICFILTERPARAMETERS_H_

#include "oops/util/parameters/ConfigurationParameter.h"
#include "ufo/filters/FilterParametersBase.h"

namespace ufo {

/// \brief A subclass of FilterParametersBase storing the values of all filter options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// \note **Do not use this class if you're writing a new filter** because the
/// ConfigurationParameter class does not perform any parameter validation. Instead, define a new
/// subclass of FilterParametersBase storing each parameter accepted by your filter in a separate
/// (Optional/Required)Parameter object. Existing filters using GenericFilterParameters should also
/// be refactored in this way when time allows.
class GenericFilterParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericFilterParameters, FilterParametersBase)
 public:
  oops::ConfigurationParameter config{this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_GENERICFILTERPARAMETERS_H_

/*
 * (C) Copyright 2021- UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/sfcpcorrected/ObsSfcPCorrectedParameters.h"

namespace ufo {

// -----------------------------------------------------------------------------
// Required for enum Parameters
constexpr char SfcPCorrectionTypeParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<SfcPCorrectionType>
          SfcPCorrectionTypeParameterTraitsHelper::namedValues[];

}  // namespace ufo

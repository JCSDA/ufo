/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <utility>

#include "ufo/obslocalization/ObsLocParameters.h"

namespace ufo {

constexpr char DistanceTypeParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<DistanceType> DistanceTypeParameterTraitsHelper::namedValues[];

constexpr char SearchMethodParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<SearchMethod> SearchMethodParameterTraitsHelper::namedValues[];

constexpr double ObsLocParameters::radius_earth;

}  // namespace ufo

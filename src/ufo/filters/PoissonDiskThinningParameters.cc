/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/PoissonDiskThinningParameters.h"

#include <utility>

namespace ufo {
constexpr char ExclusionVolumeShapeParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<ExclusionVolumeShape>
  ExclusionVolumeShapeParameterTraitsHelper::namedValues[];
}  // namespace ufo

/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "ufo/profile/ProfileDataBase.h"

#include "ufo/utils/Flags.h"

namespace ufo {
  ProfileDataBase::ProfileDataBase(ioda::ObsSpace &obsdb,
                                   const ProfileConsistencyCheckParameters &options,
                                   const ProfileIndices &profileIndices)
    : obsdb_(obsdb),
      options_(options),
      profileIndices_(profileIndices)
  {}
}  // namespace ufo

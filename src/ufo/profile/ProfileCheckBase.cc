/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileCheckBase.h"

namespace ufo {
  ProfileCheckBase::ProfileCheckBase(const ProfileConsistencyCheckParameters &options,
                                     const ProfileIndices &profileIndices,
                                     const ProfileData &profileData,
                                     ProfileFlags &profileFlags,
                                     ProfileCheckValidator &profileCheckValidator)
    : options_(options),
      profileIndices_(profileIndices),
      profileData_(profileData),
      profileFlags_(profileFlags),
      profileCheckValidator_(profileCheckValidator)
  {}
}  // namespace ufo


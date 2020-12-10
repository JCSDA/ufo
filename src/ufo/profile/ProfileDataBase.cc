/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/profile/ProfileDataBase.h"

namespace ufo {
  ProfileDataBase::ProfileDataBase(ProfileDataHandler &profileDataHandler)
    : profileDataHandler_(profileDataHandler)
  {}
}

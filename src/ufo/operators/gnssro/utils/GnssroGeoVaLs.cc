/*
 * (C) Copyright 2025 Space Sciences and Engineering, LLC (dba PlanetiQ).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * author Steve Marshall (smarshall@planetiq.com)
 */
#include <utility>
#include "oops/base/Variable.h"
#include "oops/util/Logger.h"
#include "ufo/operators/gnssro/utils/GnssroGeoVaLs.h"

namespace ufo {

// oops variables used in GeoVaLs interface.
static const oops::Variable sTempVar(GnssroGeoVaLs::VAR_TEMP);
static const oops::Variable sSphumVar(GnssroGeoVaLs::VAR_SPHUM);
static const oops::Variable sPresVar(GnssroGeoVaLs::VAR_PRES);
static const oops::Variable sGphVar(GnssroGeoVaLs::VAR_GPH);

// -----------------------------------------------------------------------------
GnssroGeoVaLs::GnssroGeoVaLs(const GeoVaLs & cgv)
  : nprofiles_(validateProfiles(cgv))
  , nlevs_(validateLevels(cgv))
  , temp_()
  , sphum_()
  , pres_()
  , gph_()
{
  temp_.reserve(nprofiles_);
  sphum_.reserve(nprofiles_);
  pres_.reserve(nprofiles_);
  gph_.reserve(nprofiles_);

  for (std::size_t profIdx = 0; profIdx < nprofiles_; ++profIdx)
  {
    getProfile(sTempVar, profIdx, temp_, cgv);
    getProfile(sSphumVar, profIdx, sphum_, cgv);
    getProfile(sPresVar, profIdx, pres_, cgv);
    getProfile(sGphVar, profIdx, gph_, cgv);
  }
}
// -----------------------------------------------------------------------------

GnssroGeoVaLs::~GnssroGeoVaLs()
{
}
// -----------------------------------------------------------------------------

void GnssroGeoVaLs::saveTo(GeoVaLs & gv)
{
  putProfiles(sTempVar, temp_, gv);
  putProfiles(sSphumVar, sphum_, gv);
  putProfiles(sPresVar, pres_, gv);
  putProfiles(sGphVar, gph_, gv);
}
// -----------------------------------------------------------------------------

std::size_t GnssroGeoVaLs::validateProfiles(const GeoVaLs& cgv)
{
  std::size_t nprofTemp = cgv.nprofiles(sTempVar);
  std::size_t nprofSphum = cgv.nprofiles(sSphumVar);
  std::size_t nprofPres = cgv.nprofiles(sPresVar);
  std::size_t nprofGph = cgv.nprofiles(sGphVar);
  if (nprofTemp != nprofSphum || nprofTemp != nprofPres || nprofTemp != nprofGph)
  {
    throw eckit::BadValue("GnssroGeoVaLs: number of profiles for " + sTempVar.name()
            + ", " + sSphumVar.name() + ", " + sPresVar.name() + ", and " + sGphVar.name()
            + " must be the same", Here());
  }
  return nprofTemp;
}
// -----------------------------------------------------------------------------

std::size_t GnssroGeoVaLs::validateLevels(const GeoVaLs& cgv)
{
  std::size_t nlevsTemp = cgv.nlevs(sTempVar);
  std::size_t nlevsSphum = cgv.nlevs(sSphumVar);
  std::size_t nlevsPres = cgv.nlevs(sPresVar);
  std::size_t nlevsGph = cgv.nlevs(sGphVar);
  if (nlevsTemp != nlevsSphum || nlevsTemp != nlevsPres || nlevsTemp != nlevsGph)
  {
    throw eckit::BadValue("GnssroGeoVaLs: number of vertical levels for " + sTempVar.name()
            + ", " + sSphumVar.name() + ", " + sPresVar.name() + ", and " + sGphVar.name()
            + " must be the same", Here());
  }
  return nlevsTemp;
}
// -----------------------------------------------------------------------------

void GnssroGeoVaLs::getProfile(
        const oops::Variable& var,
        std::size_t profIdx,
        Profiles_ & profiles,
        const GeoVaLs & cgv)
{
  ProfilePtr_ profPtr = std::make_unique<Profile_>(nlevs_, 0.0);
  cgv.getProfile(*profPtr, var, profIdx);
  profiles.push_back(std::move(profPtr));
}
// -----------------------------------------------------------------------------

void GnssroGeoVaLs::putProfiles(const oops::Variable& var, const Profiles_& profiles, GeoVaLs& gv)
{
  for (std::size_t profIdx = 0; profIdx < nprofiles_; ++profIdx)
  {
    const ProfilePtr_& profPtr = profiles.at(profIdx);
    gv.putProfile(*profPtr, var, profIdx);
  }
}
// -----------------------------------------------------------------------------

}  // namespace ufo


/*
 * CC0 - No Copyright (Public Domain)
 * 
 * The person who associated a work with this deed has dedicated the work to the public domain by waiving all of his or her rights to the work worldwide under copyright law, including all related and neighboring rights, to the extent allowed by law.
 *
 * You can copy, modify, distribute and perform the work, even for commercial purposes, all without asking permission. See Other Information below.
 *
 * Other Information
 * In no way are the patent or trademark rights of any person affected by CC0, nor are the rights that other persons may have in the work or in how the work is used, such as publicity or privacy rights.
 *
 * The person who associated a work with this deed makes no warranties about the work, and disclaims liability for all uses of the work, to the fullest extent permitted by applicable law.
 *
 * When using or citing the work, you should not imply endorsement by the author or the affirmer.
 */

#include "ufo/filters/ObsPolygonCheck.h"

#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include <boost/geometry.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
  // Read the entire contents of a text file and return it as a string.
  std::string readWholeFile(const std::string &filename) {
    return (std::ostringstream() << std::ifstream(filename).rdbuf()).str();
  }
}

namespace ufo {

// -----------------------------------------------------------------------------

ObsPolygonCheck::ObsPolygonCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "ObsPolygonCheck constructor" << std::endl;
  oops::Log::debug() << "ObsPolygonCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsPolygonCheck::~ObsPolygonCheck() {
  oops::Log::trace() << "ObsPolygonCheck destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPolygonCheck::applyFilter(const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ObsPolygonCheck applyFilter start" << std::endl;

  namespace bg = boost::geometry;

  // point_t = a boost::geometry point as a longitude-latitude pair in degrees
  using point_t = bg::model::point<double, 2, bg::cs::geographic<bg::degree>>;

  // polygon_t = a closed polygon of those points
  using polygon_t = bg::model::polygon<point_t>;

  // Get the observation locations.
  ioda::ObsDataVector<float> lons(obsdb_, "longitude", "MetaData");
  ioda::ObsDataVector<float> lats(obsdb_, "latitude", "MetaData");

  // Read the polygon from a wkt-format file.
  polygon_t poly = bg::from_wkt<polygon_t>(readWholeFile(filename));

  // Ask boost::geometry to correct common problems in the polygon definition.
  bg::correct(poly);

  // Scan for obvious errors that bg::correct couldn't correct.
  if(std::string reason; !bg::is_valid(poly, reason)) {
    oops::Log::error() << "ObsPolygonCheck: unable to correct invalid polygon ("
                       << reason << "): " << bg::wkt(poly) << std::endl;
    return;
  }

  // Figure out which side is the inside by checking a point that is known to be inside.
  point_t insidePoint = point_t(parmeters_.insideLon, parameters_.insideLat);
  bool useThisSide = bg::within(insidePoint, poly);

  // Find all points that are on the opposite side from the "insidePoint"
  vector<bool> notInside(obsdb_.nlocs(), true);
  for (size_t iloc = 0; iloc < obsdb_.nlocs(); iloc++)
    try {
      bool inside = useThisSide == bg::within(point_t(lons[iloc], lats[iloc]), poly);
      notInside[iloc] = not inside;
    } catch(const bg::exception &ex) {
      // We only catch boost geometry exceptions here; anything else is passed through.
      oops::Log::error() << "ObsPolygonCheck: boost::geometry error: " << ex.what() << std::endl;
      // The default value for all points is "not inside" so any erroring points will be rejected.
    }

  for (auto &vec : flagged)
    vec.assign(notInside.begin(), notInside.end());

  oops::Log::trace() << "ObsPolygonCheck applyFilter complete" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPolygonCheck::print(std::ostream & os) const {
  os << "ObsPolygonCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo

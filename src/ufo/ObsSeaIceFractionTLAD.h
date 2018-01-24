/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSSEAICEFRACTIONTLAD_H_
#define UFO_OBSSEAICEFRACTIONTLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperBase.h"
#include "ObsSpace.h"
#include "UfoTrait.h"
#include "util/ObjectCounter.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVector;

  // -----------------------------------------------------------------------------
  /// Sea-ice fraction observation for  model.

  class ObsSeaIceFractionTLAD : public oops::LinearObsOperBase<UfoTrait>, 
                                private util::ObjectCounter<ObsSeaIceFractionTLAD> {
  public:
    static const std::string classname() {return "ufo::ObsSeaIceFractionTLAD";}

    ObsSeaIceFractionTLAD(const ObsSpace &, const eckit::Configuration &);    
    virtual ~ObsSeaIceFractionTLAD();

    // Obs Operators
    void setTrajectory(const GeoVaLs &, const ObsBias &);
    void obsEquivTL(const GeoVaLs &, ObsVector &, const ObsBiasIncrement &) const;
    void obsEquivAD(GeoVaLs &, const ObsVector &, ObsBiasIncrement &) const;

    // Other
    const oops::Variables & variables() const {return *varin_;}

    int & toFortran() {return keyOperSeaIceFraction_;}
    const int & toFortran() const {return keyOperSeaIceFraction_;}

  private:
    void print(std::ostream &) const;
    F90hop keyOperSeaIceFraction_;
    boost::scoped_ptr<const oops::Variables> varin_;
  };
  // -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSSEAICEFRACTIONTLAD_H_

/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSSTERICHEIGHTTLAD_H_
#define UFO_OBSSTERICHEIGHTTLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "ufo/LinearObsOperatorBase.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/marine/FortranMarine.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;

  // -----------------------------------------------------------------------------
  /// Simulated Steric height for  model.
    class ObsStericHeightTLAD : public LinearObsOperatorBase, 
                                private util::ObjectCounter<ObsStericHeightTLAD> {
  public:
      static const std::string classname() {return "ufo::ObsStericHeightTLAD";}

      ObsStericHeightTLAD(const ioda::ObsSpace &, const eckit::Configuration &);    
      virtual ~ObsStericHeightTLAD();

      // Obs Operators
      void setTrajectory(const GeoVaLs &, const ObsBias &);
      void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
      void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

      // Other
      const oops::Variables & variables() const {return *varin_;}
      int & toFortran() {return keyOperStericHeight_;}
      const int & toFortran() const {return keyOperStericHeight_;}
      //const GeoVaLs * traj_;
  
  private:
      void print(std::ostream &) const;
      F90hop keyOperStericHeight_;
      boost::scoped_ptr<const GeoVaLs> traj_;
      boost::scoped_ptr<const oops::Variables> varin_;
    };

  // -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSSTERICHEIGHTTLAD_H_

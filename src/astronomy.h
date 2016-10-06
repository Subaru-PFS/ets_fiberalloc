#ifndef ETS_ASTRONOMY_H
#define ETS_ASTRONOMY_H

#include <string>
#include "rotmatrix.h"
#include "pointing.h"

class eq2hor
  {
  private:
    const double j2000= 2451545.0;
    const double sec2rad=degr2rad/3600.;

    rotmatrix npmat;
    double lat, lon, gast, eps, jd, sunlon, altitude;

  public:
    eq2hor (double lat_obs, double lon_obs, double alt_obs, const std::string &time);
    pointing radec2altaz (const pointing &radec) const;
  };

#endif

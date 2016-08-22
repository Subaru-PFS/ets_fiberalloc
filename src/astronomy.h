#ifndef ETS_ASTRONOMY_H
#define ETS_ASTRONOMY_H

#include <string>

void eq2hor_subaru (double ra, double decl, const std::string &time,
  double &alt, double &az);

#endif

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

using namespace std;

vector<complex<double>> cconv (const vector<double> &ra, const vector<double> &dec,
  double tel_ra, double tel_dec, double psi, const string &time);

PYBIND11_MODULE(pycconv,m)
  {
  namespace py = pybind11;
  using namespace pybind11::literals;
  m.doc() = "Python interface for coordinate transforms needed by ETS";

  m.def("cconv", &cconv,
    "converts Ra/Dec to PFS focal plane coordinates (very crudely)"
    "Args:\n"
    "  ra      : Target rectascensions (in deg)\n"
    "  dec     : Target declinations (in deg)\n",
    "  tel_ra  : rectascension of the telescope optical axis (deg)\n"
    "  tel_dec : declination of the telescope optical axis (deg)\n"
    "  psi     : position angle of the telescope (deg)\n"
    "  time    : observation time in UTC-compliant form,\n"
    "            e.g. '2016-04-03T08:00:00Z'\n"
    "Returns:\n"
    "  a list of complex numbers representing x/y focal plane coordinates\n"
    "  (in mm)",
    "ra"_a, "dec"_a, "tel_ra"_a, "tel_dec"_a, "psi"_a, "time"_a);
  }

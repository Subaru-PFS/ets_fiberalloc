#include <iostream>
#include <vector>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ets_tools.h"

using namespace std;

namespace py = pybind11;

class ETShelper
  {
  private:
    vector<vector<size_t>> f2t,t2f;
    vector<Target> tgt;
    pointing center;

  public:
    ETShelper(const vector<string> &t_id,
              const vector<double> &t_ra,
              const vector<double> &t_dec,
              double tel_ra, double tel_dec, double posang, string time)
      {
      center=radec2ptg(tel_ra, tel_dec);
      planck_assert(t_ra.size()==t_dec.size(),"array size mismatch");
      planck_assert(t_ra.size()==t_id.size(),"array size mismatch");
      for (size_t i = 0; i < t_ra.size(); i++)
        tgt.emplace_back(radec2ptg(t_ra[i],t_dec[i]),1.,t_id[i],0);
      eq2hor eqtest(obs_lat, obs_lon, obs_height, time);
      center = eqtest.radec2altaz(center);
      for (auto &t:tgt)
        t.altaz=eqtest.radec2altaz(t.radec);
      targetToPFI(tgt, center, posang);
      fpraster raster=tgt2raster(tgt,100,100);
      calcMappings(tgt,raster,f2t,t2f);
      }
    vector<vector<double>> getFiberpos() const
      {
      vector<vector<double>> res;
      for (size_t i=0; i<nfiber; ++i)
        {
        vec2 fpos=id2fiberpos(i);
        res.push_back({fpos.x,fpos.y});
        }
      return res;
      }
    map<string,vector<double>> getTgtpos() const
      {
      map<string,vector<double>> res;
      for (auto t:tgt)
        res[t.id]={t.pos.x,t.pos.y};
      return res;
      }
    map<string,vector<size_t>> getVis() const
      {
      map<string,vector<size_t>> res;
      for (size_t i=0; i<t2f.size(); ++i)
        res[tgt[i].id]=t2f[i];
      return res;
      }
  };

PYBIND11_PLUGIN(pyETS)
  {
  using namespace pybind11::literals;
  py::module m("pyETS",
    "Python interface for some of the ETS C++ functionality");

  py::class_<ETShelper>(m, "ETShelper")
    .def(py::init<const vector<string>&,const vector<double>&,
      const vector<double> &, double, double, double, string >(),
      "Args:\n"
      "  t_id   : Target IDs\n"
      "  t_ra   : Target rectascensions (in degrees)\n"
      "  t_dec  : target declinations (in degrees)\n"
      "  tel_ra : telescope pointing rectascension (in degrees)\n"
      "  tel_dec: telescope declination (in degrees)\n"
      "  posang : telescope position angle (in degrees)\n"
      "  time   : obervation time in the format 'yyyy-mm-ddThh:mm:ssZ'",
      "t_id"_a,"t_ra"_a,"t_dec"_a,"tel_ra"_a,"tel_dec"_a,"posang"_a,"time"_a)
    .def("getFiberpos", &ETShelper::getFiberpos,
      "returns the positions (in mm on the PFI plane) of the Cobra centers")
    .def("getTgtpos", &ETShelper::getTgtpos,
      "returns the positions (in mm on the PFI plane) of the visible targets")
    .def("getVis", &ETShelper::getVis,
      "returns a list of the visible targets and the fibers that can observe them")
;
  return m.ptr();
  }

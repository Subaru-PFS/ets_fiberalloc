#include <iostream>
#include <vector>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ets_tools.h"
#include "ets_assigners.h"

using namespace std;

namespace {

namespace py = pybind11;

/*! Discard targets that are too far away from the PFS. */
vector<size_t> select_observable (const vector<Target> &tgt, double safety)
  {
  vector<vec2> fpos(nfiber);
  for (size_t i=0; i<fpos.size(); ++i) fpos[i]=id2fiberpos(i);
  fpraster raster(fpos,100,100);
  vector<size_t> res;
  for (size_t i=0; i<tgt.size(); ++i)
    if (raster.anyIn(tgt[i].pos,rmax+safety))
      res.push_back(i);
  return res;
  }

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
              const vector<double> &t_time,
              const vector<int> &t_pri,
              double tel_ra, double tel_dec, double posang, string time)
      {
      center=radec2ptg(tel_ra, tel_dec);
      planck_assert(t_ra.size()==t_dec.size(),"array size mismatch");
      planck_assert(t_ra.size()==t_id.size(),"array size mismatch");
      for (size_t i=0; i<t_ra.size(); i++)
        tgt.emplace_back
          (radec2ptg(t_ra[i],t_dec[i]),t_time[i],t_id[i],t_pri[i]);
      eq2hor eqtest(obs_lat, obs_lon, obs_height, time);
      center = eqtest.radec2altaz(center);
      for (auto &t:tgt)
        t.altaz=eqtest.radec2altaz(t.radec);
      targetToPFI(tgt, center, posang);
      fpraster raster=tgt2raster(tgt,100,100);
      calcMappings(tgt,raster,f2t,t2f);
      }
    map<string,vector<double>> getFiberpos() const
      {
      map<string,vector<double>> res;
      for (size_t i=0; i<nfiber; ++i)
        {
        vec2 fpos=id2fiberpos(i);
        res[dataToString(i+1)] = {fpos.x,fpos.y};
        }
      return res;
      }
    map<string,vector<double>> getDotpos() const
      {
      map<string,vector<double>> res;
      for (size_t i=0; i<nfiber; ++i)
        {
        vec2 fpos=id2dotpos(i);
        res[dataToString(i+1)] = {fpos.x,fpos.y};
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
    map<string,vector<string>> getVis() const
      {
      map<string,vector<string>> res;
      for (size_t i=0; i<t2f.size(); ++i)
        {
        vector<string> tmp;
        for (auto j : t2f[i])
          tmp.push_back(dataToString(j+1));
        res[tgt[i].id]=tmp;
        }
      return res;
      }
    map<string,string> getObservation(string assigner) const
      {
      map<string,string> res;
      vector<size_t> tid, fid;
      vector<Target> tgt1(tgt);
      vector<size_t> idx = select_observable (tgt1, r_kernel);
      vector<Target> tgt2;
      for (auto i:idx)
        tgt2.push_back(tgt1[i]);
      if (!tgt2.empty())
        {
        unique_ptr<FiberAssigner> pass=make_assigner(assigner);
        pass->assign(tgt2,tid,fid);
        }
      for (size_t i=0; i<tid.size(); ++i)
        res[tgt[idx[tid[i]]].id] = dataToString(fid[i]+1);
      return res;
      }
  };

} // unnamed namespace

PYBIND11_PLUGIN(pyETS)
  {
  using namespace pybind11::literals;
  py::module m("pyETS",
    "Python interface for some of the ETS C++ functionality");

  py::class_<ETShelper>(m, "ETShelper")
    .def(py::init<const vector<string>&,const vector<double>&,
      const vector<double> &, const vector<double> &, const vector<int> &,
      double, double, double, string>(),
      "Args:\n"
      "  t_id   : Target IDs\n"
      "  t_ra   : Target rectascensions (in degrees)\n"
      "  t_dec  : target declinations (in degrees)\n"
      "  t_time : requested target observation times (in seconds) (unused)\n"
      "  t_pri  : target priorities\n"
      "  tel_ra : telescope pointing rectascension (in degrees)\n"
      "  tel_dec: telescope pointing declination (in degrees)\n"
      "  posang : telescope position angle (in degrees)\n"
      "  time   : obervation time in the format 'yyyy-mm-ddThh:mm:ssZ'",
      "t_id"_a,"t_ra"_a,"t_dec"_a,"t_time"_a,"t_pri"_a,"tel_ra"_a,"tel_dec"_a,
      "posang"_a,"time"_a)
    .def("getFiberpos", &ETShelper::getFiberpos,
      "returns the positions (in mm on the PFI plane) of the Cobra centers")
    .def("getDotpos", &ETShelper::getDotpos,
      "returns the positions (in mm on the PFI plane) of the dots")
    .def("getTgtpos", &ETShelper::getTgtpos,
      "returns the positions (in mm on the PFI plane) of the visible targets")
    .def("getVis", &ETShelper::getVis,
      "returns a list of the visible targets and the fibers that can observe them")
    .def("getObservation", &ETShelper::getObservation,
      "returns the targets that would be chosen by the provided assigner "
      "algorithm, as well as the fibers observing them","assigner"_a)
;
  return m.ptr();
  }

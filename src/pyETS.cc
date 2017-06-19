#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "string_utils.h"
#include "ets.h"
#include "ets_helpers.h"

using namespace std;

namespace {

namespace py = pybind11;

map<string,vector<double>> getAllCobras()
  {
  map<string,vector<double>> res;
  auto c=makeCobras();
  for (size_t i=0; i<c.size(); ++i)
    res[dataToString(i+1)] = {c[i].center.x(),c[i].center.y(),c[i].l1,c[i].l2,
                              c[i].dotpos.x(),c[i].dotpos.y(),c[i].rdot};
  return res;
  }

class ETShelper
  {
  private:
    vector<Target> tgt;
    vector<Cobra> cobras;

  public:
    ETShelper(const vector<string> &t_id,
              const vector<double> &t_x,
              const vector<double> &t_y,
              const vector<double> &t_time,
              const vector<int> &t_pri,
              const map<string,vector<double>> &cbr)
      {
      planck_assert((t_id.size()==t_x.size())
                  &&(t_id.size()==t_y.size())
                  &&(t_id.size()==t_time.size())
                  &&(t_id.size()==t_pri.size()), "vector length mismatch");
      for (size_t i=0; i< t_id.size(); ++i)
        tgt.emplace_back(vec2(t_x[i],t_y[i]),t_time[i],t_id[i],t_pri[i]);

      if (cbr.empty())
        cobras=makeCobras();
      else
        for (auto it=cbr.begin(); it!=cbr.end(); ++it)
          {
          auto d(it->second);
          cobras.emplace_back(vec2(d[0],d[1]),d[2],d[3],vec2(d[4],d[5]),d[6]);
          }
      }
    map<string,vector<string>> getVis() const
      {
      auto t2f = getT2F(tgt,cobras);
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
    map<string,string> getObs(const string &assigner) const
      {
      map<string,string> res;
      vector<size_t> tid, fid;
      if (!tgt.empty())
        getObservation(tgt,cobras,assigner,tid,fid);
      for (size_t i=0; i<tid.size(); ++i)
        res[tgt[tid[i]].id] = dataToString(fid[i]+1);
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
      const map<string,vector<double>> &>(),
      "Args:\n"
      "  t_id   : Target IDs\n"
      "  t_x    : Target x coordinates on the focal plane (in mm)\n"
      "  t_y    : Target y coordinates on the focal plane (in mm)\n"
      "  t_time : requested target observation times (in seconds) (unused)\n"
      "  t_pri  : target priorities\n"
      "  cbr    : dictionary of cobras (if empty, a default set is used\n"
      "t_id"_a,"t_x"_a,"t_y"_a,"t_time"_a,"t_pri"_a,"cbr"_a)
    .def("getVis", &ETShelper::getVis,
      "returns a list of the visible targets and the fibers that can observe them")
    .def("getObs", &ETShelper::getObs,
      "returns the targets that would be chosen by the provided assigner "
      "algorithm, as well as the fibers observing them","assigner"_a)
;
  m.def("getAllCobras", &getAllCobras,
    "returns a dictionary containing all cobras. A cobra is defined by a\n"
    "6-tuple of real numbers (unit is mm):\n"
    "  x position of the cobra center\n"
    "  y position pf the cobra center\n"
    "  length l1\n"
    "  length l2\n"
    "  x position of the dot\n"
    "  y position of the dot\n"
    "  dot radius");
  return m.ptr();
  }

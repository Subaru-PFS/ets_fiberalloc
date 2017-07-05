#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

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

map<size_t,vector<size_t>> getVis(const vector<complex<double>> &t_pos,
                                  const map<string,vector<double>> &cbr)
  {
  vector<Target> tgt;
  for (size_t i=0; i<t_pos.size(); ++i)
    tgt.emplace_back(t_pos[i],1.,1);

  vector<Cobra> cobras;
  if (cbr.empty())
    cobras=makeCobras();
  else
    for (auto it=cbr.begin(); it!=cbr.end(); ++it)
      {
      auto d(it->second);
      cobras.emplace_back(vec2(d[0],d[1]),d[2],d[3],vec2(d[4],d[5]),d[6]);
      }
  auto tmp = getT2F(tgt,cobras);
  map<size_t,vector<size_t>> res;
  for (size_t i=0; i< tmp.size(); ++i)
    if (tmp[i].size()>0) res[i]=tmp[i];
  return res;
  }

map<size_t,size_t> getObs(const vector<complex<double>> &t_pos,
                          const vector<double> &t_time,
                          const vector<int> &t_pri,
                          const map<string,vector<double>> &cbr,
                          const string &assigner)
  {
  planck_assert((t_pos.size()==t_time.size())
              &&(t_pos.size()==t_pri.size()), "vector length mismatch");
  vector<Target> tgt;
  for (size_t i=0; i<t_pos.size(); ++i)
    tgt.emplace_back(t_pos[i],t_time[i],t_pri[i]);

  vector<Cobra> cobras;
  if (cbr.empty())
    cobras=makeCobras();
  else
    for (auto it=cbr.begin(); it!=cbr.end(); ++it)
      {
      auto d(it->second);
      cobras.emplace_back(vec2(d[0],d[1]),d[2],d[3],vec2(d[4],d[5]),d[6]);
      }

  vector<size_t> tid, fid;
  if (!tgt.empty())
    getObservation(tgt,cobras,assigner,tid,fid);
  map<size_t,size_t> res;
  for (size_t i=0; i<tid.size(); ++i)
    res[tid[i]] = fid[i];
  return res;
  }

} // unnamed namespace

PYBIND11_PLUGIN(pyETS)
  {
  using namespace pybind11::literals;
  py::module m("pyETS",
    "Python interface for some of the ETS C++ functionality");

  m.def("getVis", &getVis,
    "returns a list of the visible targets and the fibers that can observe them.\n"
    "Args:\n"
    "  t_pos  : Target x/y coordinates on the focal plane (in mm)\n"
    "  cbr    : dictionary of cobras (if empty, a default set is used)\n",
    "t_pos"_a, "cbr"_a);
  m.def("getObs", &getObs,
    "performs an assignment step and returns a dictionary containing the\n"
    "observed target numbers and the assigned cobra numbers\n"
    "Args:\n"
    "  t_pos   : Target x/y coordinates on the focal plane (in mm)\n"
    "  t_time  : requested target observation times (in seconds) (unused)\n"
    "  t_pri   : target priorities\n"
    "  cbr     : dictionary of cobras (if empty, a default set is used)\n"
    "  assigner: algorithm to do the assignment. Must be one of 'naive',"
    "            'draining', 'draining_closest' or 'new'\n",
    "t_pos"_a, "t_time"_a, "t_pri"_a, "cbr"_a, "assigner"_a);
  m.def("getAllCobras", &getAllCobras,
    "returns a dictionary containing all cobras of an idealized instrument\n"
    "configuration. A cobra is defined by a 6-tuple of real numbers (unit is\n"
    "mm):\n"
    "  x position of the cobra center\n"
    "  y position pf the cobra center\n"
    "  length l1 (link between center of cobra and 'elbow')\n"
    "  length l2 (link between 'elbow' and tip\n"
    "  x position of the dot\n"
    "  y position of the dot\n"
    "  dot radius");
  return m.ptr();
  }

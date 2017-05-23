/*
 *  This file is part of ets_fiber_assigner.
 *
 *  ets_fiber_assigner is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  ets_fiber_assigner is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ets_fiber_assigner; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  ets_fiber_assigner is being developed at the Max-Planck-Institut fuer
 *  Astrophysik.
 */

/*! \file ets_demo.cc
 *  Copyright (C) 2016-2017 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <regex>

#include "ets_tools.h"
#include "ets_assigners.h"

using namespace std;

namespace {

/*! Discard targets that are too far away from the PFS. */
vector<size_t> select_observable (const vector<Target> &tgt,
  const vector<Cobra> &cobras, double safety)
  {
  vector<vec2> fpos(cobras.size());
  double rmax=0;
  for (size_t i=0; i<fpos.size(); ++i)
    {
    fpos[i]=cobras[i].center;
    rmax=max(rmax,cobras[i].rmax);
    }
  fpraster raster(fpos,100,100);
  vector<size_t> res;
  for (size_t i=0; i<tgt.size(); ++i)
    if (raster.anyIn(tgt[i].pos,rmax+safety))
      res.push_back(i);
  return res;
  }

void single_exposure(const vector<Target> &tgt, const vector<Cobra> &cobras,
  const pointing &center, double posang, const FiberAssigner &ass,
  vector<size_t> &tid, vector<size_t> &fid)
  {
  tid.clear(); fid.clear();
  vector<Target> tgt1(tgt);
  targetToPFI(tgt1, center, posang);
#if 1
  vector<size_t> idx = select_observable (tgt1, cobras, r_kernel);
  vector<Target> tgt2;
  for (auto i:idx)
    tgt2.push_back(tgt1[i]);
  if (!tgt2.empty()) ass.assign(tgt2,cobras,tid,fid);
  for (size_t i=0; i<tid.size(); ++i)
    tid[i]=idx[tid[i]];
#else
  if (!tgt1.empty()) ass.assign(tgt1,tid,fid);
#endif
  }

} // unnamed namespace

void optimal_exposure(const vector<Target> &tgt, const vector<Cobra> &cobras,
  pointing &center, double dptg, int nptg, double &posang, double dposang,
  int nposang, const FiberAssigner &ass, vector<size_t> &tid,
  vector<size_t> &fid)
  {
  double posang0=posang;
  tid.clear(); fid.clear();
  vec3 vcenter(center);
  vec3 vdx=crossprod(vcenter,vec3(0,0,1));
  if (vdx.SquaredLength()==0.) // center lies at a pole
    vdx=vec3(1,0,0);
  else
    vdx.Normalize();
  vec3 vdy=crossprod(vcenter,vdx);
  //FIXME: make this user-definable!
  for (int idx=0; idx<nptg; ++idx)
    for (int idy=0; idy<nptg; ++idy)
      for (int ida=0; ida<nposang; ++ida)
        {
        double dx=-dptg+2*dptg*(idx+0.5)/nptg;
        double dy=-dptg+2*dptg*(idy+0.5)/nptg;
        double da=-dposang+2*dposang*(ida+0.5)/nposang;
        pointing newcenter(vcenter+(vdx*dx+vdy*dy));
        double newposang=posang0+da;
        vector<size_t> tid2,fid2;
        single_exposure (tgt, cobras, newcenter, newposang, ass, tid2, fid2);
        if (tid2.size()>tid.size())
          { tid=tid2; fid=fid2; center=newcenter; posang=newposang; }
        }
  }

namespace {

void strip (vector<Target> &tgt, const vector<size_t> &remove, double time)
  {
  vector<bool> obs(tgt.size(),false);
  for (size_t i=0; i<remove.size(); ++i)
    obs[remove[i]]=true;
  vector<Target> t2;
  for (size_t i=0; i<obs.size(); ++i)
    if (!obs[i])
      t2.push_back(tgt[i]);
    else
      if (tgt[i].time>time+1e-7)
        {
        t2.push_back(tgt[i]);
        t2.back().time-=time;
        }
  tgt.swap(t2);
  }

template<typename T> string toString(const T&val, int w)
  {
  ostringstream o;
  o<<setw(w)<<val;
  return o.str();
  }
template<typename T> string toString(const T&val, int w, int p)
  {
  ostringstream o;
  o<<fixed<<setw(w)<<setprecision(p)<<val;
  return o.str();
  }

void subprocess (const vector<Target> &tgt, const vector<Cobra> &cobras,
  const pointing &center0,
  double dptg, int nptg, double posang0, double dposang, int nposang,
  int n_exposures, ofstream &fout, const FiberAssigner &ass)
  {
  vector<Target> tgt1=tgt;
  double ttime=0., acc=0., time2=0.;
  for (const auto &t: tgt)
    ttime += t.time;
  cout << endl << "Total observation time: " << ttime << endl;
  size_t cnt=0;
  cout << endl << "tile # | fiber allocation fraction | "
                  "total observation fraction | time"
       << endl;
  for (int expo=0; expo<n_exposures; ++expo)
    {
    pointing center(center0);
    double posang(posang0);
    vector<size_t> tidmax, fidmax;
    optimal_exposure(tgt1, cobras, center, dptg, nptg, posang, dposang, nposang,
      ass, tidmax, fidmax);
    if (tidmax.empty()) break; // stop if no more fibers could be assigned
    double time=tgt1[tidmax[0]].time;
    for (const auto i: tidmax)
      if (tgt1[i].time<time) time=tgt1[i].time;
    time2+=time;
    acc+=tidmax.size()*time;
    targetToPFI(tgt1,center,posang);
    if (fout.is_open())
      {
      fout << "# Exposure " << cnt << ": duration " << time << "s, "
        "AZ: " << rad2degr*center.phi << ", ALT " << 90-rad2degr*center.theta
        << " PA: " << rad2degr*posang << endl
        << "# Target    Fiber          X          Y         "
           "RA        DEC     exp" << endl;
      //FIXME: add PFI coordinates
      for (size_t i=0; i<tidmax.size(); ++i)
        fout << tgt1[tidmax[i]].id << toString(fidmax[i]+1,10)
        << " " << toString(tgt1[tidmax[i]].pos.x,10,5)
        << " " << toString(tgt1[tidmax[i]].pos.y,10,5)
        << " " << toString(rad2degr*tgt1[tidmax[i]].radec.phi,10,5)
        << " " << toString(90-rad2degr*tgt1[tidmax[i]].radec.theta,10,5)
        << " " << toString(cnt,7)
        << endl;
      }
    cout << toString(cnt++,6)
         << toString(tidmax.size()/double(cobras.size()),18,5)
         << toString(acc/ttime,28,5)
         << toString(time2,20,0) << endl;
    strip (tgt1,tidmax,time);
    }
  }

/*! Reads targets from the ASCII file \a name and returns them in a \a vector.
    The returned coordinates are RA/DEC in degrees. */
vector<Target> readTargets (const string &name, const string &time)
  {
  int lineno=0;
  vector<Target> res;
  ifstream inp(name);
  planck_assert (inp,"Could not open target file '"+name+"'.");
  while (inp)
    {
    string line;
    getline(inp, line);
    ++lineno;
    // remove potential carriage returns at the end of the line
    line=line.substr(0,line.find("\r"));
    line=line.substr(0,line.find("#"));
    line=trim(line);
    if (line.size()>0)
      {
      istringstream iss(line);
      double x,y,time;
      string id;
      int pri;
      iss >> id >> x >> y >> time >> pri;
      if (iss)
        res.emplace_back(radec2ptg(x,y),time,id,pri);
      else
        cerr << "Warning: unrecognized format in '" << name << "', line "
             << lineno << ":\n" << line << endl;
      }
    }
  eq2hor eqtest(obs_lat, obs_lon, obs_height, time);
  for (auto &t:res)
    t.altaz=eqtest.radec2altaz(t.radec);

  return res;
  }

void process(const string &name, int n_exposures,
  const pointing &center, double dptg, int nptg, double posang, double dposang,
  int nposang, const string &out, const string &time,
  const FiberAssigner &ass)
  {
  vector<Target> tgt=readTargets(name,time);
  vector<Cobra> cobras(makeCobras());
  eq2hor eqtest(obs_lat, obs_lon, obs_height, time);
  pointing center_altaz(eqtest.radec2altaz(center));

#if 0
  {
  vector<Target> tmp(tgt), tgt2;
  targetToPFI(tmp, center_altaz, posang);
  for (size_t i=0; i<tmp.size(); ++i)
    if (tmp[i].pos.dsq(vec2(0,0))<190*190)
      tgt2.push_back(tgt[i]);
  tgt.swap(tgt2);
  }
#endif

  ofstream fout;
  if (out!="")
    { fout.open(out); planck_assert(fout,"error opening output file"); }
  subprocess (tgt, cobras, center_altaz, dptg, nptg, posang, dposang, nposang,
    n_exposures, fout, ass);
  }

/*! Finds the smallest circle enclosing all locations in \a tgt and returns
    its center. Used to find a telescope pointing that hits the given target
    list. Only for temporary use. */
vec3 getCenter(const vector<Target> &tgt)
  {
  vector<vec3> tmp;
  for (auto t:tgt)
    tmp.push_back(vec3(t.radec));
  double dummy;
  vec3 res;
  find_enclosing_circle(tmp,res,dummy);
  pointing pcnt(res);
  cout << "center of data set: RA " << rad2degr*pcnt.phi
       << ", DEC " << 90-rad2degr*pcnt.theta << endl;
  return res;
  }

} // unnamed namespace

int main(int argc, const char ** argv)
  {
  map<string,string> paramdict;
  parse_cmdline_equalsign (argc, argv, paramdict);
  paramfile params (paramdict);

  unique_ptr<FiberAssigner> pass=make_assigner(params.find<string>("assigner"));
  pointing center;
  if (params.param_present("ra")||params.param_present("dec"))
    center=radec2ptg (params.find<double>("ra"), params.find<double>("dec"));
  else
    center=pointing(getCenter(readTargets(params.find<string>("input"),
           params.find<string>("time"))));

  double posang=degr2rad*params.find<double>("posang",0.);
  double dposang=degr2rad*params.find<double>("dposang",4.);
  int nposang=params.find<int>("nposang",5);
  double dptg=degr2rad*params.find<double>("dptg",4./320.);// should roughly correspond to 4mm in PFI plane
  int nptg=params.find<int>("nptg",5);
  process (params.find<string>("input"),
    params.find<int>("n_exposures",1),center,dptg,nptg,posang,dposang,nposang,
    params.find<string>("output",""),params.find<string>("time"),*pass);
  }

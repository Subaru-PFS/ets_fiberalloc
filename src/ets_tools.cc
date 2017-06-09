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

/*! \file ets_tools.cc
 *  Copyright (C) 2017 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "ets_tools.h"

using namespace std;

namespace {

double colldist=2;

vec2 elbow_pos(const Cobra &c, const vec2 &tip)
  {
  vec2 pt(tip-c.center);
  double apt=abs(pt);
  auto rot = pt/apt;
  const double l1=2.375, l2=2.375; // fixed for the moment
  double x=(l2*l2-l1*l1-apt*apt)/(-2*apt);
  double y=-sqrt(l1*l1-x*x);
  return vec2(x,y)*rot + c.center;
  }

bool line_segment_collision (const vec2 &x1, const vec2 &x2, const vec2 &y,
  double dist)
  {
  // interpret x1 as origin
  auto p2 = x2-x1;
  auto q = y-x1;
  double ap2=abs(p2);
  // rotate p2 to lie on positive real axis
  auto rot = conj(p2)/ap2;
  q*=rot;
  if (q.real()<=0) return norm(q)<=dist*dist;
  if (q.real()>=ap2) return norm(q-ap2)<=dist*dist;
  return abs(q.imag())<=dist;
  }

} // unnamed namespace

void setCollisionDistance(double dist)
  {colldist=dist; }

void targetToPFI(vector<Target> &tgt, const pointing &los, double psi)
  {
  // altitude and azimuth of North celestial pole:
  pointing altaz_ncp (halfpi-obs_lat,0.);
  vec3 z{los}, skypole(altaz_ncp);
  vec3 x=(skypole-z*dotprod(z,skypole)).Norm();
  vec3 y=crossprod(z,x);
  double cpsi=cos(psi),spsi=sin(psi);
  const double a0=0., a1=-3.2e2, a2=-1.37e1, a3=-7.45e0;
  for (auto&& t:tgt)
    {
    vec3 pos=t.altaz;
    vec3 xp=pos-y*dotprod(pos,y);
    vec3 yp=pos-x*dotprod(pos,x);
    vec2 pnew (atan2(dotprod(xp,x),dotprod(xp,z))*rad2degr,
               atan2(dotprod(yp,y),dotprod(yp,z))*rad2degr);
    rotate (pnew,cpsi,spsi);
    double rsq=norm(pnew);
    t.pos = vec2 ((a3*rsq*rsq+a2*rsq+a1)*pnew.x()+a0,
                  (-a3*rsq*rsq-a2*rsq-a1)*pnew.y()+a0);
    }
  }

namespace {

/*! Computes the central fiber position in PFI coordinates, given the fiber ID.
    Fiber ID is zero-based throughout this code, i.e. ranging from 0 to 2393. */
vec2 id2fiberpos(int id)
  {
  int field=id/(57*14);
  id-=field*57*14;
  int module=id/57;
  int cobra=id-module*57;
  const double vspace=sqrt(0.75); // cos(30deg)
  vec2 res (-vspace*(1.+2*module+(cobra&1)),
            0.5+module-0.5*cobra);
  if (field==1) rotate(res,-vspace,-0.5);
  if (field==2) rotate(res,vspace,-0.5);
  res*=8;
  return res;
  }

/*! Computes the position of a dot center in PFI coordinates, given a fiber ID.
    Fiber ID is zero-based throughout this code, i.e. ranging from 0 to 2393. */
vec2 id2dotpos(int id) // id is assumed to be in [0; 2394[
  {
  const double dot_shift_y = 2.35;
  vec2 res=id2fiberpos(id);
  res+=vec2(0.,dot_shift_y);
  return res;
  }

}

vector<Cobra> makeCobras()
  {
  constexpr size_t nfiber=3*57*14;
  constexpr double rmax=4.75; // maximum radius of a fiber patrol area
  constexpr double dotdist=1.375; // radius of the dot blocking area
  vector<Cobra> res;
  for (size_t i=0; i<nfiber; ++i)
    res.emplace_back(id2fiberpos(i),rmax,id2dotpos(i),dotdist);
  return res;
  }

void calcMappings (const vector<Target> &tgt, const vector<Cobra> &cobras,
  const fpraster &raster,
  vector<vector<size_t>> &f2t, vector<vector<size_t>> &t2f)
  {
  f2t=vector<vector<size_t>>(cobras.size());
  for (size_t i=0; i<cobras.size(); ++i)
    {
    const auto &c(cobras[i]);
    vector<size_t> tmp=raster.query(c.center,c.rmax);
    for (auto j : tmp)
      if (c.dotpos.dsq(tgt[j].pos)>=c.rdot*c.rdot) f2t[i].push_back(j);
    }
  t2f=vector<vector<size_t>>(tgt.size());
  for (size_t i=0; i<f2t.size(); ++i)
    for (auto t : f2t[i])
      t2f[t].push_back(i);

 // for (auto &v: f2t) sort(v.begin(), v.end());
 // for (auto &v: t2f) sort(v.begin(), v.end());
  }

void cleanup (const vector<Target> &tgt, const vector<Cobra> &cobras,
  const fpraster &raster, vector<vector<size_t>> &f2t,
  vector<vector<size_t>> &t2f, int fiber, int itgt)
  {
  // remove everything related to the selected fiber
  for (auto curtgt : f2t[fiber]) stripout(t2f[curtgt],fiber);
  f2t[fiber].clear();
  // remove target
  for (auto j : t2f[itgt]) stripout(f2t[j],itgt);
  t2f[itgt].clear();
  // remove everything in "lower arm" area of the assigned cobra
  if (colldist>0.)
    {
    vec2 tippos (tgt[itgt].pos),
         elbowpos(elbow_pos(cobras[fiber], tippos));
    vector<size_t> tmp=raster.query(0.5*(tippos+elbowpos),colldist+2.375/2.); //FIXME
// classical version:
//    vector<size_t> tmp=raster.query(tippos,colldist); //FIXME
    for (auto i : tmp)
      if (line_segment_collision (elbowpos, tippos, tgt[i].pos, colldist))
        {
        for (auto j : t2f[i]) stripout(f2t[j],i);
        t2f[i].clear();
        }
    }
//  checkMappings(tgt,f2t,t2f);
  }

fpraster tgt2raster (const vector<Target> &tgt, int nx, int ny)
  {
  vector<vec2> tpos(tgt.size());
  for (size_t i=0; i<tgt.size(); ++i) tpos[i]=tgt[i].pos;
  return fpraster (tpos,nx,ny);
  }

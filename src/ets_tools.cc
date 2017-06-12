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

fpraster tgt2raster (const vector<Target> &tgt, int nx, int ny)
  {
  vector<vec2> tpos(tgt.size());
  for (size_t i=0; i<tgt.size(); ++i) tpos[i]=tgt[i].pos;
  return fpraster (tpos,nx,ny);
  }

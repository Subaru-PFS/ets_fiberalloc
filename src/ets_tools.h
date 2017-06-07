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

/*! \file ets_tools.h
 *  Copyright (C) 2017 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef ETS_TOOLS_H
#define ETS_TOOLS_H

#include "error_handling.h"
#include "string_utils.h"
#include "geom_utils.h"
#include "lsconstants.h"
#include "pointing.h"
#include "paramfile.h"
#include "rotmatrix.h"
#include "astronomy.h"

/*! Simple class for storing a position in a 2D plane. */
class vec2
  {
  public:
    double x,y;

    vec2() {}
    vec2(double x_, double y_) : x(x_), y(y_) {}
    double dsq(const vec2 &b) const
      { return (x-b.x)*(x-b.x) + (y-b.y)*(y-b.y); }
  };

/*! Simple class containing all relevant properties of a PFS observation
    target. */
class Target
  {
  public:
    vec2 pos;
    pointing radec, altaz;
    double time;
    int pri;
    std::string id;

    Target (const pointing &radec_,double time_,const std::string &id_,int pri_)
      : radec(radec_), time(time_), pri(pri_), id(id_) {}
  };

class Cobra
  {
  public:
    vec2 center;
    double rmax;
    vec2 dotpos;
    double rdot;

    Cobra(const vec2 &center_, double rmax_, const vec2 &dotpos_, double rdot_)
      : center(center_), rmax(rmax_), dotpos(dotpos_), rdot(rdot_) {}
  };

std::vector<Cobra> makeCobras();

void setCollisionDistance(double dist);

constexpr double r_kernel=4.75; // radius of the priority function kernel
// Latitude and longitude of the obseratory (radian)
const double obs_lat=(19+49/60.+32/3600.)*degr2rad,
                 obs_lon=-(155+28/60.+34/3600.)*degr2rad;
// Height of the observatory (meters above sea level)
const double obs_height=4139.;

/*! Class providing efficient queries for locations on a 2D plane */
class fpraster
  {
  private:
    double x0, y0, x1, y1, idx, idy;
    size_t nx, ny;
    std::vector<std::vector<size_t>> data;
    std::vector<vec2> loc;

    size_t indexx (double x) const
      { return size_t(std::max(0,std::min(int(nx)-1,int((x-x0)*idx)))); }
    size_t indexy (double y) const
      { return size_t(std::max(0,std::min(int(ny)-1,int((y-y0)*idy)))); }
    size_t index (const vec2 &pos) const
      { return indexx(pos.x) + nx*indexy(pos.y); }

  public:
    /*! Constructs an \a fpraster with \a nx_ bins in x direction and
        \a ny_ bins in y direction, and sorts the entries in \a loc_ into this
        structure. */
    fpraster (const std::vector<vec2> &loc_, size_t nx_, size_t ny_) :
      nx(nx_),ny(ny_),data(nx*ny), loc(loc_)
      {
      planck_assert ((nx>0) && (ny>0), "bad array sizes");
      planck_assert(loc.size()>0,"input array too small");
      x0=x1=loc[0].x;
      y0=y1=loc[0].y;
      for (size_t i=1; i<loc.size(); ++i)
        {
        x0=std::min(x0,loc[i].x); x1=std::max(x1,loc[i].x);
        y0=std::min(y0,loc[i].y); y1=std::max(y1,loc[i].y);
        }
      if (x0==x1) x1+=1e-9;
      if (y0==y1) y1+=1e-9;
      idx=nx/(x1-x0);
      idy=ny/(y1-y0);
      for (size_t i=0; i<loc.size(); ++i)
        data[index(loc[i])].push_back(i);
      }
    /*! Returns the indices of all \a loc entries that lie within a circle of
        radius \a rad around \a center. */
    std::vector<size_t> query(const vec2 &center, double rad) const
      {
      std::vector<size_t> res;
      if ((center.x<x0-rad)||(center.x>x1+rad)
        ||(center.y<y0-rad)||(center.y>y1+rad))
        return res;
      double rsq=rad*rad;
      size_t i0=indexx(center.x-rad), i1=indexx(center.x+rad),
            j0=indexy(center.y-rad), j1=indexy(center.y+rad);
      for (size_t j=j0; j<=j1; ++j)
        for (size_t i=i0; i<=i1; ++i)
          for (auto k : data[i+nx*j])
            if (center.dsq(loc[k])<=rsq) res.push_back(k);
      return res;
      }
    bool anyIn (const vec2 &center, double rad) const
      {
      if ((center.x<x0-rad)||(center.x>x1+rad)
        ||(center.y<y0-rad)||(center.y>y1+rad))
        return false;
      double rsq=rad*rad;
      size_t i0=indexx(center.x-rad), i1=indexx(center.x+rad),
             j0=indexy(center.y-rad), j1=indexy(center.y+rad);
      for (size_t j=j0; j<=j1; ++j)
        for (size_t i=i0; i<=i1; ++i)
          for (auto k : data[i+nx*j])
            if (center.dsq(loc[k])<=rsq) return true;
      return false;
      }
  };

/*! Converts RA/DEC in degrees to colatitude/longitude in radians. */
inline pointing radec2ptg (double ra, double dec)
  { return pointing((90-dec)*degr2rad,ra*degr2rad); }

inline void rotate (vec2 &pos, double sa, double ca)
  {
  vec2 t{pos};
  pos.x = ca*t.x - sa*t.y;
  pos.y = sa*t.x + ca*t.y;
  }

/*! Converts target coordinates from alt/az to PFI coordinates in
    millimeters, given a telescope pointing and orientation.
    \note This is still very preliminary, incomplete and approximate! */
void targetToPFI(std::vector<Target> &tgt, const pointing &los, double psi);

/*! Remove a given value from a vector of integers. Assert that exactly one
    value was removed. */
inline void stripout (std::vector<size_t> &v, size_t val)
  {
  size_t s1=v.size();
  v.erase(remove(v.begin(),v.end(),val));
  planck_assert(v.size()+1==s1,"oops");
  }

fpraster tgt2raster (const std::vector<Target> &tgt, int nx, int ny);

#if 0
/*! Diagnostic function to check for inconsistencies in the fiber->target and
    target->fiber mappings. */
void checkMappings (const vector<Target> &tgt,
  const vector<vector<size_t>> &/*f2t*/, const vector<vector<size_t>> &t2f)
  {
  vector<bool> tmp(tgt.size(), false);
  for (const auto &i : t2f)
    {
    planck_assert(i.size()<=3,"too many fibers");
    set<size_t> tmp2;
    for (auto j:i)
      {
      planck_assert(tmp2.find(j)==tmp2.end(),"multiply defined fiber");
      tmp2.insert(j);
      }
    }
  }
#endif

/*! Computes the fiber->target and target->fiber mappings. */
void calcMappings (const std::vector<Target> &tgt,
  const std::vector<Cobra> &cobras, const fpraster &raster,
  std::vector<std::vector<size_t>> &f2t, std::vector<std::vector<size_t>> &t2f);

/*! Given a target index \a itgt and a fiber index \a fiber observing this
    target, remove all references to \a itgt from the mappings and also remove
    all targets that lie in the blocking area around \a itgt and all targets
    exclusively visible from \a fiber. */
void cleanup (const std::vector<Target> &tgt, const std::vector<Cobra> &cobras,
  const fpraster &raster,
  std::vector<std::vector<size_t>> &f2t, std::vector<std::vector<size_t>> &t2f,
  int fiber, int itgt);

inline double kernelfunc(double rsq)
  {
  // simple parabola - quick but probably not optimal
  return std::max(0.,r_kernel*r_kernel-rsq);
//  return sqrt(max(0.,r_kernel*r_kernel-rsq)); // linear decrease
//  return exp(-9*rsq/(r_kernel*r_kernel)); // Gaussian kernel
  }

#endif

/*
 *  This file is part of ets_fiber_assigner.
 *
 *  ets_finer_assigner is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  ets_finer_assigner is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ets_finer_assigner; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  ets_finer_assigner is being developed at the Max-Planck-Institut fuer
 *  Astrophysik.
 */

/*! \file ets_demo.cc */

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <regex>

#include "error_handling.h"
#include "string_utils.h"
#include "geom_utils.h"
#include "lsconstants.h"
#include "pointing.h"
#include "paramfile.h"
#include "rotmatrix.h"
#include "astronomy.h"

using namespace std;

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
    int id;

    Target (const pointing &radec_, double time_, int id_, int pri_)
      : radec(radec_), time(time_), pri(pri_), id(id_) {}
  };

namespace {

/*! Priority queue that allows changing the priority of its entries after
    its creation
    Originally developed for Gadget 4. */
template <typename T, typename Compare=std::less<T>> class pqueue
  {
  private:
    Compare comp;
    struct node_t
      {
      T pri;
      size_t pos;
      };
    vector<node_t> nodes;
    vector<size_t> idx;

    void sift_up (size_t i)
      {
      size_t moving_node = idx[i];
      T moving_pri = nodes[moving_node].pri;

      for (size_t parent_node=i>>1;
           (i>1) && (comp(nodes[idx[parent_node]].pri,moving_pri));
           i=parent_node, parent_node=i>>1)
        {
        idx[i] = idx[parent_node];
        nodes[idx[i]].pos=i;
        }

      idx[i] = moving_node;
      nodes[idx[i]].pos=i;
      }

    size_t maxchild(size_t i) const
      {
      size_t child_node = i<<1;

      if (child_node>=idx.size())
        return 0;

      if (((child_node+1)<idx.size())
          && (comp(nodes[idx[child_node]].pri,nodes[idx[child_node+1]].pri)))
        child_node++; /* use right child instead of left */

      return child_node;
      }

    void sift_down(size_t i)
      {
      size_t moving_node = idx[i];
      T moving_pri = nodes[moving_node].pri;

      size_t child_node;
      while ((child_node = maxchild(i))
             && comp(moving_pri,nodes[idx[child_node]].pri))
        {
        idx[i] = idx[child_node];
        nodes[idx[i]].pos=i;
        i = child_node;
        }

      idx[i] = moving_node;
      nodes[idx[i]].pos=i;
      }

    /*! Rearranges the internal data structure to ensure the heap property. */
    void heapify()
      {
      size_t startnode=idx.size()>>1;
      for (size_t i=startnode; i>=1; --i)
        sift_down(i);
      }

  public:
    /*! Constructs a \a pqueue of size \a n. All priorities are set to zero. */
    pqueue (size_t n) : nodes(n), idx(n+1)
      {
      idx[0]=0;
      for (size_t i=0; i<n; ++i)
        {
        nodes[i]= {0.,i+1};
        idx[i+1]=i;
        }
      }
    /*! Constructs a \a pqueue with priorities taken from \a pri. */
    pqueue (const vector<T> &pri) : nodes(pri.size()), idx(pri.size()+1)
      {
      idx[0]=0;
      for (size_t i=0; i<pri.size(); ++i)
        {
        nodes[i]= {pri[i],i+1};
        idx[i+1]=i;
        }
      heapify();
      }

    /*! Sets the priority of the entry \a pos to \a new_pri. The heap is rebuilt
        automatically. */
    void set_priority(T new_pri, size_t pos)
      {
      T old_pri = nodes[pos].pri;
      nodes[pos].pri=new_pri;
      size_t posn = nodes[pos].pos;
      comp(old_pri,new_pri) ? sift_up(posn) : sift_down(posn);
      }

    /*! Returns the priority of the entry \a pos. */
    T priority(size_t pos) const
      { return nodes[pos].pri; }
    /*! Returns the lowest priority in the queue. */
    T top_priority() const
      { return nodes[idx[1]].pri; }
    /*! Returns entry with the lowest priority in the queue. */
    size_t top() const
      { return idx[1]; }
  };
constexpr size_t nfiber=3*57*14;
constexpr double rmax=4.75; // maximum radius of a fiber patrol area
constexpr double r_kernel=4.75; // radius of the priority function kernel
constexpr double dotdist=1.375; // radius of the dot blocking area
constexpr double colldist=2; // minimum distance between fiber positioners

/*! Class providing efficient queries for locations on a 2D plane */
class fpraster
  {
  private:
    double x0, y0, x1, y1, idx, idy;
    size_t nx, ny;
    vector<vector<size_t>> data;
    vector<vec2> loc;

    size_t indexx (double x) const
      { return size_t(max(0,min(int(nx)-1,int((x-x0)*idx)))); }
    size_t indexy (double y) const
      { return size_t(max(0,min(int(ny)-1,int((y-y0)*idy)))); }
    size_t index (const vec2 &pos) const
      { return indexx(pos.x) + nx*indexy(pos.y); }

  public:
    /*! Constructs an \a fpraster with \a nx_ bins in x direction and
        \a ny_ bins in y direction, and sorts the entries in \a loc_ into this
        structure. */
    fpraster (const vector<vec2> &loc_, size_t nx_, size_t ny_) :
      nx(nx_),ny(ny_),data(nx*ny), loc(loc_)
      {
      planck_assert ((nx>0) && (ny>0), "bad array sizes");
      planck_assert(loc.size()>0,"input array too small");
      x0=x1=loc[0].x;
      y0=y1=loc[0].y;
      for (size_t i=1; i<loc.size(); ++i)
        {
        x0=min(x0,loc[i].x); x1=max(x1,loc[i].x);
        y0=min(y0,loc[i].y); y1=max(y1,loc[i].y);
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
    vector<size_t> query(const vec2 &center, double rad) const
      {
      vector<size_t> res;
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

void rotate (vec2 &pos, double sa, double ca)
  {
  vec2 t{pos};
  pos.x = ca*t.x - sa*t.y;
  pos.y = sa*t.x + ca*t.y;
  }

/*! Converts target coordinates from alt/az to PFI coordinates in
    millimeters, given a telescope pointing and orientation.
    \note This is still very preliminary, incomplete and approximate! */
void targetToPFI(vector<Target> &tgt, const pointing &los, double psi)
  {
  vec3 z{los}, sky{0,0,1};
  vec3 x=(sky-z*dotprod(z,sky)).Norm();
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
    double rsq=pnew.x*pnew.x+pnew.y*pnew.y;
    t.pos.x= (a3*rsq*rsq+a2*rsq+a1)*pnew.x+a0;
    t.pos.y= (-a3*rsq*rsq-a2*rsq-a1)*pnew.y+a0;
    }
  }

/*! Computes the central fiber position in PFI coordinates, given the fiber ID.
    Fiber ID is zero-based throughout this code, i.e. ranging from 0 to 2393. */
vec2 id2fiberpos(int id)
  {
  int field=id/(57*14);
  id-=field*57*14;
  int module=id/57;
  int cobra=id-module*57;
  const double vspace=sqrt(0.75); // cos(30deg)
  vec2 res;
  res.y=0.5+module-0.5*cobra;
  res.x=-vspace*(1.+2*module+(cobra&1));
  if (field==1) rotate(res,-vspace,-0.5);
  if (field==2) rotate(res,vspace,-0.5);
  res.x*=8; res.y*=8;
  return res;
  }

/*! Computes the position of a dot center in PFI coordinates, given a fiber ID.
    Fiber ID is zero-based throughout this code, i.e. ranging from 0 to 2393. */
vec2 id2dotpos(int id) // id is assumed to be in [0; 2394[
  {
  const double dot_shift_y = 2.35;
  vec2 res=id2fiberpos(id);
  res.y+=dot_shift_y;
  return res;
  }

/*! Remove a given value from a vector of integers. Assert that exactly one
    value was removed. */
inline void stripout (vector<size_t> &v, size_t val)
  {
  size_t s1=v.size();
  v.erase(remove(v.begin(),v.end(),val));
  planck_assert(v.size()+1==s1,"oops");
  }
fpraster tgt2raster (const vector<Target> &tgt, int nx, int ny)
  {
  vector<vec2> tpos(tgt.size());
  for (size_t i=0; i<tgt.size(); ++i) tpos[i]=tgt[i].pos;
  return fpraster (tpos,nx,ny);
  }

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
void calcMappings (const vector<Target> &tgt, const fpraster &raster,
  vector<vector<size_t>> &f2t, vector<vector<size_t>> &t2f)
  {
  f2t=vector<vector<size_t>>(nfiber);
  for (size_t i=0; i<nfiber; ++i)
    {
    vec2 fp=id2fiberpos(i), dp=id2dotpos(i);
    vector<size_t> tmp=raster.query(fp,rmax);
    for (auto j : tmp)
      if (dp.dsq(tgt[j].pos)>=dotdist*dotdist) f2t[i].push_back(j);
    }
  t2f=vector<vector<size_t>>(tgt.size());
  for (size_t i=0; i<f2t.size(); ++i)
    for (auto t : f2t[i])
      t2f[t].push_back(i);

 // for (auto &v: f2t) sort(v.begin(), v.end());
 // for (auto &v: t2f) sort(v.begin(), v.end());
  }

/*! Given a target index \a itgt and a fiber index \a fiber observing this
    target, remove all references to \a itgt from the mappings and also remove
    all targets that lie in the blocking area around \a itgt and all targets
    exclusively visible from \a fiber. */
void cleanup (const vector<Target> &tgt, const fpraster &raster,
  vector<vector<size_t>> &f2t, vector<vector<size_t>> &t2f, int fiber, int itgt)
  {
  // remove everything related to the selected fiber
  for (auto curtgt : f2t[fiber]) stripout(t2f[curtgt],fiber);
  f2t[fiber].clear();
  // remove target and everything in blocking area
  vector<size_t> tmp=raster.query(tgt[itgt].pos,colldist);
  for (auto i : tmp)
    {
    for (auto j : t2f[i]) stripout(f2t[j],i);
    t2f[i].clear();
    }
//  checkMappings(tgt,f2t,t2f);
  }

inline double kernelfunc(double rsq)
  {
  // simple parabola - quick but probably not optimal
  return max(0.,r_kernel*r_kernel-rsq);
//  return sqrt(max(0.,r_kernel*r_kernel-rsq)); // linear decrease
//  return exp(-9*rsq/(r_kernel*r_kernel)); // Gaussian kernel
  }

struct pq_entry
  {
  double prox;
  int pri;

  pq_entry() : prox(0.), pri(0) {}
  pq_entry(double prox_,int pri_) : prox(prox_), pri(pri_) {}

  bool operator< (const pq_entry &other) const
    {
    if (pri!=other.pri) return pri>other.pri;
    return prox<other.prox;
    }
  };

pqueue<pq_entry> calc_pri(const vector<Target> &tgt,
  const vector<vector<size_t>> &t2f, const fpraster &raster)
  {
  vector<pq_entry> pri(tgt.size());
  for (size_t i=0; i<tgt.size(); ++i)
    {
    if (t2f[i].size()>0)
      {
      vector<size_t> ngb = raster.query(tgt[i].pos,r_kernel);

      for (auto j : ngb)
        {
        if (i==j)
          pri[i].prox+=tgt[i].time*tgt[i].time*kernelfunc(0.);
        if (i<j)
          {
          double tmp=tgt[i].time*tgt[j].time
                    *kernelfunc(tgt[i].pos.dsq(tgt[j].pos));
          pri[i].prox+=tmp;
          pri[j].prox+=tmp;
          }
        }
      }
    }
  for (size_t i=0; i<tgt.size(); ++i)
    pri[i].pri=tgt[i].pri;
  pqueue<pq_entry> res(pri);
  return res;
  }

void fix_priority(const vector<Target> &tgt, const vector<vector<size_t>> &t2f,
  const fpraster &raster, size_t itgt, pqueue<pq_entry> &pri)
  {
  vector<size_t> ngb = raster.query(tgt[itgt].pos,r_kernel);
  for (auto j : ngb)
    if ((!t2f[j].empty())||(pri.priority(j).prox!=0.))
      {
      pq_entry tpri=pri.priority(j);
      tpri.prox-=tgt[j].time*tgt[itgt].time
                *kernelfunc(tgt[itgt].pos.dsq(tgt[j].pos));
      pri.set_priority(tpri,j);
      }
  }

class FiberAssigner
  {
  public:
    /*! Assign target from \a tgt to fibers according to a given strategy.
        On return, \a tid and \a fid contain target resp. fiber IDs for the
        assigned targets. Target IDs range from 0 to \a tgt.size()-1, fiber
        IDs from 0 to 2393. */
    virtual void assign (const vector<Target> &tgt,
      vector<size_t> &tid, vector<size_t> &fid) const = 0;

    virtual ~FiberAssigner() {}
  };

int maxpri_in_fiber (size_t fiber, const vector<Target> &tgt,
  const vector<vector<size_t>> &f2t)
  {
  planck_assert(!f2t[fiber].empty(), "searching in empty fiber");
  size_t idx=0;
  int maxpri = tgt[f2t[fiber][idx]].pri;
  for (size_t j=1; j<f2t[fiber].size(); ++j)
    if (tgt[f2t[fiber][idx]].pri<maxpri)
      { maxpri=tgt[f2t[fiber][idx]].pri; idx=j; }
  return f2t[fiber][idx];
  }

class NaiveAssigner: public FiberAssigner
  {
  /*! Naive assignment algorithm: iterate over all fibers, and if a fiber has
      targets in its patrol area, assign the target with the highest priority
      to it. */
  virtual void assign (const vector<Target> &tgt,
    vector<size_t> &tid, vector<size_t> &fid) const
    {
    tid.clear(); fid.clear();
    fpraster raster=tgt2raster(tgt,100,100);
    vector<vector<size_t>> f2t,t2f;
    calcMappings(tgt,raster,f2t,t2f);

    for (size_t fiber=0; fiber<f2t.size(); ++fiber)
      {
      if (f2t[fiber].empty()) continue;
      int itgt = maxpri_in_fiber(fiber,tgt,f2t);
      tid.push_back(itgt);
      fid.push_back(fiber);
      cleanup (tgt, raster, f2t, t2f, fiber, itgt);
      }
    }
  };

class DrainingAssigner: public FiberAssigner
  {
  /*! Assignment strategy modeled after Morales et al. 2012: MNRAS 419, 1187
      find the fiber(s) with the smallest number of observable targets >0;
      for the first of the returned fibers, assign the target with highest
      priority to it; repeat until no more targets are observable. */
  virtual void assign (const vector<Target> &tgt,
    vector<size_t> &tid, vector<size_t> &fid) const
    {
    tid.clear(); fid.clear();
    fpraster raster=tgt2raster(tgt,100,100);
    vector<vector<size_t>> f2t,t2f;
    calcMappings(tgt,raster,f2t,t2f);

    size_t maxtgt=0;
    for (const auto &f:f2t)
      maxtgt=max(maxtgt,f.size());

    while (true)
      {
      int fiber=-1;
      size_t mintgt=maxtgt+1;
      for (size_t i=0; i<f2t.size(); ++i)
        if ((f2t[i].size()<mintgt)&&(f2t[i].size()>0))
          { fiber=i; mintgt=f2t[i].size(); }
      if (fiber==-1) break; // assignment done
      int itgt = maxpri_in_fiber(fiber,tgt,f2t);
      tid.push_back(itgt);
      fid.push_back(fiber);
      cleanup(tgt,raster,f2t,t2f,fiber,itgt);
      }
    }
  };

class NewAssigner: public FiberAssigner
  {
  /*! Assignment strategy with the goal of reducing inhomogeneity in the
      target distribution: assign a priority to each target that depends on
      the distance of all other targets in its close vicinity; process targets
      in order of decreasing priority and assign them to fibers, if possible.
      After each assignment, update the priority of the remaining targets. */
  virtual void assign (const vector<Target> &tgt,
    vector<size_t> &tid, vector<size_t> &fid) const
    {
    tid.clear(); fid.clear();
    fpraster raster=tgt2raster(tgt,100,100);
    vector<vector<size_t>> f2t,t2f;
    calcMappings(tgt,raster,f2t,t2f);
    pqueue<pq_entry> pri=calc_pri(tgt,t2f,raster);

    while (true)
      {
      if (pri.top_priority().pri==(1<<30)) break;
      size_t itgt=pri.top();
      if (t2f[itgt].empty())
        { pri.set_priority(pq_entry(0.,(1<<30)),itgt); continue; }
      size_t ifib=0, mintgt=f2t[t2f[itgt][ifib]].size();
      for (size_t i=1; i<t2f[itgt].size(); ++i)
        if (f2t[t2f[itgt][i]].size()<mintgt)
          { ifib=i; mintgt=f2t[t2f[itgt][i]].size(); }
      int fiber=t2f[itgt][ifib];
      tid.push_back(itgt);
      fid.push_back(fiber);
      cleanup(tgt,raster,f2t,t2f,fiber,itgt);
      fix_priority(tgt,t2f,raster,itgt,pri);
      }
    }
};

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

void single_exposure(const vector<Target> &tgt, const pointing &center,
  double posang, const FiberAssigner &ass,
  vector<size_t> &tid, vector<size_t> &fid)
  {
  tid.clear(); fid.clear();
  vector<Target> tgt1(tgt);
  targetToPFI(tgt1, center, posang);
#if 1
  vector<size_t> idx = select_observable (tgt1, r_kernel);
  vector<Target> tgt2;
  for (auto i:idx)
    tgt2.push_back(tgt1[i]);
  if (!tgt2.empty()) ass.assign(tgt2,tid,fid);
  for (size_t i=0; i<tid.size(); ++i)
    tid[i]=idx[tid[i]];
#else
  if (!tgt1.empty()) ass.assign(tgt1,tid,fid);
#endif
  }

} // unnamed namespace

void optimal_exposure(const vector<Target> &tgt, pointing &center, double dptg,
  int nptg, double &posang, double dposang, int nposang,
  const FiberAssigner &ass, vector<size_t> &tid, vector<size_t> &fid)
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
        single_exposure (tgt, newcenter, newposang, ass, tid2, fid2);
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

void subprocess (const vector<Target> &tgt, const pointing &center0,
  double dptg, int nptg, double posang0, double dposang, int nposang,
  double fract, ofstream &fout, const FiberAssigner &ass)
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
  while (true)
    {
    pointing center(center0);
    double posang(posang0);
    vector<size_t> tidmax, fidmax;
    optimal_exposure(tgt1, center, dptg, nptg, posang, dposang, nposang,
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
      fout << "Exposure " << cnt << ": duration " << time << "s, "
        "AZ: " << rad2degr*center.phi << ", ALT " << 90-rad2degr*center.theta
        << " PA: " << rad2degr*posang << endl
        << "  Target     Fiber          X          Y         "
           "RA        DEC" << endl;
      //FIXME: add PFI coordinates
      for (size_t i=0; i<tidmax.size(); ++i)
        fout << toString(tgt1[tidmax[i]].id,8) << toString(fidmax[i]+1,10)
        << " " << toString(tgt1[tidmax[i]].pos.x,10,5)
        << " " << toString(tgt1[tidmax[i]].pos.y,10,5)
        << " " << toString(rad2degr*tgt1[tidmax[i]].radec.phi,10,5)
        << " " << toString(90-rad2degr*tgt1[tidmax[i]].radec.theta,10,5)
        << endl;
      }
    cout << toString(cnt++,6)
         << toString(tidmax.size()/double(nfiber),18,5)
         << toString(acc/ttime,28,5)
         << toString(time2,20,0) << endl;
    if (acc/ttime>fract) break;
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
      string id0;
      int pri;
      iss >> id0 >> x >> y >> time >> pri;
      if (iss)
        {
        planck_assert((id0.length()>2) && (id0.substr(0,2)=="ID"),
          "identifier not starting with 'ID'");
        int id=stringToData<int>(id0.substr(2));
        res.emplace_back(radec2ptg(x,y),time,id,pri);
        }
      else
        cerr << "Warning: unrecognized format in '" << name << "', line "
             << lineno << ":\n" << line << endl;
      }
    }
  eq2hor eqtest((19+49/60.+32/3600.)*degr2rad, -(155+28/60.+34/3600.)*degr2rad, 4139.,time);
  for (auto &t:res)
    t.altaz=eqtest.radec2altaz(t.radec);

  return res;
  }

void process(const string &name, double fract,
  const pointing &center, double dptg, int nptg, double posang, double dposang,
  int nposang, const string &out, const string &time,
  const FiberAssigner &ass)
  {
  vector<Target> tgt=readTargets(name,time);
  eq2hor eqtest((19+49/60.+32/3600.)*degr2rad, -(155+28/60.+34/3600.)*degr2rad, 4139.,time);
  pointing center_altaz(eqtest.radec2altaz(center));
cout << "center altaz: "<<center_altaz << endl;
  {
  vector<Target> tmp(tgt), tgt2;
  targetToPFI(tmp, center_altaz, posang);
  for (size_t i=0; i<tmp.size(); ++i)
    if (tmp[i].pos.dsq(vec2(0,0))<190*190)
      tgt2.push_back(tgt[i]);
  tgt.swap(tgt2);
  }
  ofstream fout;
  if (out!="")
    { fout.open(out); planck_assert(fout,"error opening output file"); }
  subprocess (tgt, center_altaz, dptg, nptg, posang, dposang, nposang,
    fract, fout, ass);
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

  unique_ptr<FiberAssigner> pass;
  string assignerName=params.find<string>("assigner");
  if (assignerName=="naive")
    pass=make_unique<NaiveAssigner>();
  else if (assignerName=="draining")
    pass=make_unique<DrainingAssigner>();
  else if (assignerName=="new")
    pass=make_unique<NewAssigner>();
  else
    planck_fail("unknown assigner");
  pointing center;
  if (params.param_present("ra")||params.param_present("dec"))
    center=radec2ptg (params.find<double>("ra"), params.find<double>("dec"));
  else
    center=pointing(getCenter(readTargets(params.find<string>("input"),
           params.find<string>("time"))));

cout << "center radec: "<<center << endl;

  double posang=degr2rad*params.find<double>("posang",0.);
  double dposang=degr2rad*params.find<double>("dposang",4.);
  int nposang=params.find<int>("nposang",5);
  double dptg=degr2rad*params.find<double>("dptg",4./320.);// should roughly correspond to 4mm in PFI plane
  int nptg=params.find<int>("nptg",5);
  process (params.find<string>("input"),
    params.find<double>("fract"),center,dptg,nptg,posang,dposang,nposang,
    params.find<string>("output",""),params.find<string>("time"),*pass);
  }

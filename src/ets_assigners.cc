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

/*! \file ets_assigners.cc
 *  Copyright (C) 2017 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <random>
#include "ets_assigners.h"

using namespace std;

namespace {

constexpr int randseed=42;

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
    std::vector<node_t> nodes;
    std::vector<size_t> idx;

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
    pqueue (const std::vector<T> &pri) : nodes(pri.size()), idx(pri.size()+1)
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

} // unnamed namespace

void setCollisionDistance(double dist)
  {colldist=dist; }

    /*! Assign target from \a tgt to fibers according to a given strategy.
        On return, \a tid and \a fid contain target resp. fiber IDs for the
        assigned targets. Target IDs range from 0 to \a tgt.size()-1, fiber
        IDs from 0 to 2393. */
class FiberAssigner
  {
  protected:
    const vector<Target> &tgt;
    const vector<Cobra> &cobras;
    vector<size_t> &tid;
    vector<size_t> &fid;

    fpraster raster;
    vector<vector<size_t>> f2t,t2f;

    /*! Computes the fiber->target and target->fiber mappings. */
    void calcMappings ()
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

    int maxpri_in_fiber (size_t fiber) const
      {
      using namespace std;
      static std::mt19937 engine{randseed};

      planck_assert(!f2t[fiber].empty(), "searching in empty fiber");
      vector<size_t> tmp;
      int maxpri = tgt[f2t[fiber][0]].pri;
      tmp.push_back(0);
      for (size_t j=1; j<f2t[fiber].size(); ++j)
        {
        if (tgt[f2t[fiber][j]].pri==maxpri)
          tmp.push_back(j);
        else if (tgt[f2t[fiber][j]].pri<maxpri)
          {
          maxpri=tgt[f2t[fiber][j]].pri;
          tmp.clear();
          tmp.push_back(j);
          }
        }
      std::uniform_int_distribution<size_t> dist(0, tmp.size() - 1);
      return f2t[fiber][tmp[dist(engine)]];
      }

    int maxpri_in_fiber_closest (size_t fiber) const
      {
      using namespace std;
      vec2 fpos=cobras[fiber].center;
      int maxpri = tgt[f2t[fiber][0]].pri;
      double mindsq=fpos.dsq(tgt[f2t[fiber][0]].pos);
      size_t idx=0;
      for (size_t j=1; j<f2t[fiber].size(); ++j)
        {
        if (tgt[f2t[fiber][j]].pri==maxpri)
          {
          if (fpos.dsq(tgt[f2t[fiber][j]].pos)<mindsq)
            { idx=j; mindsq=fpos.dsq(tgt[f2t[fiber][j]].pos); }
          }
        else if (tgt[f2t[fiber][j]].pri<maxpri)
          {
          maxpri=tgt[f2t[fiber][j]].pri;
          idx=j;
          mindsq=fpos.dsq(tgt[f2t[fiber][j]].pos);
          }
        }
      return f2t[fiber][idx];
      }

/*! Given a target index \a itgt and a fiber index \a fiber observing this
    target, remove all references to \a itgt from the mappings and also remove
    all targets that lie in the blocking area around \a itgt and all targets
    exclusively visible from \a fiber. */
    void cleanup (int fiber, int itgt)
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

  public:
    FiberAssigner(const vector<Target> &tgt_, const std::vector<Cobra> &cobras_,
    vector<size_t> &tid_, vector<size_t> &fid_)
      : tgt(tgt_), cobras(cobras_), tid(tid_), fid(fid_),
        raster(tgt2raster(tgt,100,100))
      {
      tid.clear(); fid.clear();
      calcMappings();
      }
  };

/*! Naive assignment algorithm: iterate over all fibers, and if a fiber has
    targets in its patrol area, assign the target with the highest priority
    to it. */
class NaiveAssigner: public FiberAssigner
  {
  public:
    NaiveAssigner (const vector<Target> &tgt_, const std::vector<Cobra> &cobras_,
      vector<size_t> &tid_, vector<size_t> &fid_)
      : FiberAssigner(tgt_, cobras_, tid_, fid_)
      {
      for (size_t fiber=0; fiber<f2t.size(); ++fiber)
        {
        if (f2t[fiber].empty()) continue;
        int itgt = maxpri_in_fiber(fiber);
   //     try_to_assign(itgt, fiber, f2t, t2f, tgt, cobras, tid, fid, larms);
        tid.push_back(itgt);
        fid.push_back(fiber);
        cleanup (fiber, itgt);
        }
      }
  };

  /*! Assignment strategy modeled after Morales et al. 2012: MNRAS 419, 1187
      find the fiber(s) with the smallest number of observable targets >0;
      for the first of the returned fibers, assign the target with highest
      priority to it; repeat until no more targets are observable. */
class DrainingAssigner: public FiberAssigner
  {
  public:
    DrainingAssigner (const vector<Target> &tgt_, const std::vector<Cobra> &cobras_,
      vector<size_t> &tid_, vector<size_t> &fid_)
      : FiberAssigner(tgt_, cobras_, tid_, fid_)
      {
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
        int itgt = maxpri_in_fiber(fiber);
        tid.push_back(itgt);
        fid.push_back(fiber);
        cleanup(fiber,itgt);
        }
      }
  };

  /*! Assignment strategy modeled after Morales et al. 2012: MNRAS 419, 1187
      find the fiber(s) with the smallest number of observable targets >0;
      for the first of the returned fibers, assign the target with highest
      priority to it; if there is more than one of those, choose the one closest
      to the center of the patrol area;
      repeat until no more targets are observable. */
class DrainingClosestAssigner: public FiberAssigner
  {
  public:
    DrainingClosestAssigner (const vector<Target> &tgt_, const std::vector<Cobra> &cobras_,
      vector<size_t> &tid_, vector<size_t> &fid_)
      : FiberAssigner(tgt_, cobras_, tid_, fid_)
      {
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
        int itgt = maxpri_in_fiber_closest(fiber);
        tid.push_back(itgt);
        fid.push_back(fiber);
        cleanup(fiber,itgt);
        }
      }
  };

  /*! Assignment strategy with the goal of reducing inhomogeneity in the
      target distribution: assign a priority to each target that depends on
      the distance of all other targets in its close vicinity; process targets
      in order of decreasing priority and assign them to fibers, if possible.
      After each assignment, update the priority of the remaining targets. */
class NewAssigner: public FiberAssigner
  {
  private:
    pqueue<pq_entry> calc_pri() const
      {
      std::vector<pq_entry> pri(tgt.size());
      for (size_t i=0; i<tgt.size(); ++i)
        {
        if (t2f[i].size()>0)
          {
          std::vector<size_t> ngb = raster.query(tgt[i].pos,r_kernel);

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

    void fix_priority(size_t itgt, pqueue<pq_entry> &pri)
      {
      std::vector<size_t> ngb = raster.query(tgt[itgt].pos,r_kernel);
      for (auto j : ngb)
        if ((!t2f[j].empty())||(pri.priority(j).prox!=0.))
          {
          pq_entry tpri=pri.priority(j);
          tpri.prox-=tgt[j].time*tgt[itgt].time
                    *kernelfunc(tgt[itgt].pos.dsq(tgt[j].pos));
          pri.set_priority(tpri,j);
          }
      }

  public:
    NewAssigner (const vector<Target> &tgt_, const std::vector<Cobra> &cobras_,
      vector<size_t> &tid_, vector<size_t> &fid_)
      : FiberAssigner(tgt_, cobras_, tid_, fid_)
      {
      pqueue<pq_entry> pri=calc_pri();

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
        cleanup(fiber,itgt);
        fix_priority(itgt,pri);
        }
      }
  };

void ets_assign (const std::string &name, const std::vector<Target> &tgt,
      const std::vector<Cobra> &cobras,
      std::vector<size_t> &tid, std::vector<size_t> &fid)
  {
  if (name=="naive")
    NaiveAssigner tmp(tgt, cobras, tid, fid);
  else if (name=="draining")
    DrainingAssigner tmp(tgt, cobras, tid, fid);
  else if (name=="draining_closest")
    DrainingClosestAssigner tmp(tgt, cobras, tid, fid);
  else if (name=="new")
    NewAssigner tmp(tgt, cobras, tid, fid);
  else
    planck_fail("unknown assigner");
  }

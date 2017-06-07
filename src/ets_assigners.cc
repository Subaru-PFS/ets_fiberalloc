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

#ifdef HAVE_ORTOOLS
#include "graph/min_cost_flow.h"
#endif
#include "ets_assigners.h"

using namespace std;

namespace {

constexpr int randseed=42;

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

pqueue<pq_entry> calc_pri(const std::vector<Target> &tgt,
  const std::vector<std::vector<size_t>> &t2f, const fpraster &raster)
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

void fix_priority(const std::vector<Target> &tgt,
  const std::vector<std::vector<size_t>> &t2f,
  const fpraster &raster, size_t itgt, pqueue<pq_entry> &pri)
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


int maxpri_in_fiber (size_t fiber, const std::vector<Target> &tgt,
  const std::vector<std::vector<size_t>> &f2t)
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

int maxpri_in_fiber_closest (size_t fiber, const std::vector<Target> &tgt,
  const std::vector<Cobra> &cobras, const std::vector<std::vector<size_t>> &f2t)
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

} // unnamed namespace

class NaiveAssigner: public FiberAssigner
  {
  /*! Naive assignment algorithm: iterate over all fibers, and if a fiber has
      targets in its patrol area, assign the target with the highest priority
      to it. */
  virtual void assign (const vector<Target> &tgt, const std::vector<Cobra> &cobras,
    vector<size_t> &tid, vector<size_t> &fid) const
    {
    tid.clear(); fid.clear();
    fpraster raster=tgt2raster(tgt,100,100);
    vector<vector<size_t>> f2t,t2f;
    calcMappings(tgt,cobras,raster,f2t,t2f);

    for (size_t fiber=0; fiber<f2t.size(); ++fiber)
      {
      if (f2t[fiber].empty()) continue;
      int itgt = maxpri_in_fiber(fiber,tgt,f2t);
      tid.push_back(itgt);
      fid.push_back(fiber);
      cleanup (tgt, cobras, raster, f2t, t2f, fiber, itgt);
      }
    }
  };

class DrainingAssigner: public FiberAssigner
  {
  /*! Assignment strategy modeled after Morales et al. 2012: MNRAS 419, 1187
      find the fiber(s) with the smallest number of observable targets >0;
      for the first of the returned fibers, assign the target with highest
      priority to it; repeat until no more targets are observable. */
  virtual void assign (const vector<Target> &tgt, const std::vector<Cobra> &cobras,
    vector<size_t> &tid, vector<size_t> &fid) const
    {
    tid.clear(); fid.clear();
    fpraster raster=tgt2raster(tgt,100,100);
    vector<vector<size_t>> f2t,t2f;
    calcMappings(tgt,cobras,raster,f2t,t2f);

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
      cleanup(tgt,cobras,raster,f2t,t2f,fiber,itgt);
      }
    }
  };

class DrainingClosestAssigner: public FiberAssigner
  {
  /*! Assignment strategy modeled after Morales et al. 2012: MNRAS 419, 1187
      find the fiber(s) with the smallest number of observable targets >0;
      for the first of the returned fibers, assign the target with highest
      priority to it; if there is more than one of those, choose the one closest
      to the center of the patrol area;
      repeat until no more targets are observable. */
  virtual void assign (const vector<Target> &tgt, const std::vector<Cobra> &cobras,
    vector<size_t> &tid, vector<size_t> &fid) const
    {
    tid.clear(); fid.clear();
    fpraster raster=tgt2raster(tgt,100,100);
    vector<vector<size_t>> f2t,t2f;
    calcMappings(tgt,cobras,raster,f2t,t2f);

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
      int itgt = maxpri_in_fiber_closest(fiber,tgt,cobras,f2t);
      tid.push_back(itgt);
      fid.push_back(fiber);
      cleanup(tgt,cobras,raster,f2t,t2f,fiber,itgt);
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
  virtual void assign (const vector<Target> &tgt, const std::vector<Cobra> &cobras,
    vector<size_t> &tid, vector<size_t> &fid) const
    {
    tid.clear(); fid.clear();
    fpraster raster=tgt2raster(tgt,100,100);
    vector<vector<size_t>> f2t,t2f;
    calcMappings(tgt,cobras,raster,f2t,t2f);
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
      cleanup(tgt,cobras,raster,f2t,t2f,fiber,itgt);
      fix_priority(tgt,t2f,raster,itgt,pri);
      }
    }
  };

#ifdef HAVE_ORTOOLS
class NetworkAssigner: public FiberAssigner
  {
  /*! Assignment using min-cost network flow algorithm. Currently, collisions
      are ignored. Priorities are ignored as well.
      Performance of the algorithm lies between naive and draining assigner */
  virtual void assign (const vector<Target> &tgt,
    vector<size_t> &tid, vector<size_t> &fid) const
    {
    using namespace operations_research;
    tid.clear(); fid.clear();
    vector<size_t> head,tail,xarc;
    fpraster raster=tgt2raster(tgt,100,100);
    vector<vector<size_t>> f2t,t2f;
    calcMappings(tgt,raster,f2t,t2f);
    size_t nfibers=f2t.size();
    size_t ntgt=t2f.size();
    size_t nnodes=nfibers+ntgt+2;
    size_t narcs=0;
    for (auto x: f2t) narcs+=x.size();
    narcs+=nfibers+ntgt+1;
    StarGraph graph(nnodes, narcs);
    MinCostFlow min_cost_flow(&graph);
    // Source node is 0;
    // fiber nodes are [1; nfibers];
    // target nodes are [nfibers+1; nfibers+ntgt];
    // sink node is nfibers+ntgt+1.

    // Arcs from source to fibers
    for (size_t i=0; i<nfibers; ++i)
      {
      ArcIndex arc = graph.AddArc(0, 1+i);
      min_cost_flow.SetArcUnitCost(arc, 1);
      min_cost_flow.SetArcCapacity(arc, 1);
      }
    // Arcs between fibers and targets
    for (size_t i=0; i<nfibers; ++i)
      for (size_t j=0; j<f2t[i].size(); ++j)
        {
        size_t arc=graph.AddArc(1+i, 1+nfibers+f2t[i][j]);
        xarc.push_back(arc);
        head.push_back(i); tail.push_back(f2t[i][j]);
        min_cost_flow.SetArcUnitCost(arc, 1);
        min_cost_flow.SetArcCapacity(arc, 1);
        }
    // Arcs from targets to sink
    for (size_t i=0; i<ntgt; ++i)
      {
      ArcIndex arc = graph.AddArc(1+nfibers+i, 1+nfibers+ntgt);
      min_cost_flow.SetArcUnitCost(arc, 1);
      min_cost_flow.SetArcCapacity(arc, 1);
      }
    // overflow arc
    {
    ArcIndex arc = graph.AddArc(0, 1+nfibers+ntgt);
    min_cost_flow.SetArcUnitCost(arc, 1000000);
    min_cost_flow.SetArcCapacity(arc, 1000000);
    }
    // source node
    min_cost_flow.SetNodeSupply(0, nfibers);
    // passive nodes
    for (size_t i=0; i<nfibers+ntgt; ++i)
      min_cost_flow.SetNodeSupply(i+1,0);
    // sink node
    min_cost_flow.SetNodeSupply(1+nfibers+ntgt,-nfibers);
    CHECK(min_cost_flow.Solve());
    CHECK_EQ(MinCostFlow::OPTIMAL, min_cost_flow.status());
    // Extract the solution
    for (size_t i=0; i<xarc.size(); ++i)
      if (min_cost_flow.Flow(xarc[i])>0)
        {
        fid.push_back(head[i]);
        tid.push_back(tail[i]);
        }
    CostValue total_flow_cost = min_cost_flow.GetOptimalCost();
    }
  };

#endif

unique_ptr<FiberAssigner> make_assigner(const string &name)
  {
  if (name=="naive")
    return make_unique<NaiveAssigner>();
  else if (name=="draining")
    return make_unique<DrainingAssigner>();
  else if (name=="draining_closest")
    return make_unique<DrainingClosestAssigner>();
  else if (name=="new")
    return make_unique<NewAssigner>();
#ifdef HAVE_ORTOOLS
  else if (name=="network")
    return make_unique<NetworkAssigner>();
#endif
  else
    planck_fail("unknown assigner");
  }

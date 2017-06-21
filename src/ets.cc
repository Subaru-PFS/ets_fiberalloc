#include <random>
#include "ets.h"
#include "ets_helpers.h"
#include "error_handling.h"

using namespace std;

namespace {

inline vec2 elbow_pos(const Cobra &c, const vec2 &tip)
  {
  vec2 pt(tip-c.center);
  double apt=abs(pt);
  auto rot = pt/apt;
  double x=(c.l2*c.l2-c.l1*c.l1-apt*apt)/(-2*apt);
  double y=-sqrt(c.l1*c.l1-x*x);
  return vec2(x,y)*rot + c.center;
  }

inline bool line_segment_collision (const vec2 &x1, const vec2 &x2,
  const vec2 &y, double dist)
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

/*! Class providing efficient queries for locations on a 2D plane */
class fpraster
  {
  protected:
    double x0, y0, x1, y1, idx, idy;
    size_t nx, ny;
    std::vector<std::vector<size_t>> data;
    std::vector<vec2> loc;

    size_t indexx (double x) const
      { return size_t(std::max(0,std::min(int(nx)-1,int((x-x0)*idx)))); }
    size_t indexy (double y) const
      { return size_t(std::max(0,std::min(int(ny)-1,int((y-y0)*idy)))); }
    size_t index (const vec2 &pos) const
      { return indexx(pos.x()) + nx*indexy(pos.y()); }

  public:
    fpraster (const vec2 &pmin, const vec2 &pmax, size_t nx_, size_t ny_)
      : nx(nx_),ny(ny_),data(nx*ny)
      {
      planck_assert ((nx>0) && (ny>0), "bad array sizes");
      x0=pmin.x(); y0=pmin.y();
      x1=pmax.x(); y1=pmax.y();
      idx=nx/(x1-x0);
      idy=ny/(y1-y0);
      }
    /*! Constructs an \a fpraster with \a nx_ bins in x direction and
        \a ny_ bins in y direction, and sorts the entries in \a loc_ into this
        structure. */
    fpraster (const std::vector<vec2> &loc_, size_t nx_, size_t ny_)
      : nx(nx_),ny(ny_),data(nx*ny), loc(loc_)
      {
      planck_assert ((nx>0) && (ny>0), "bad array sizes");
      planck_assert(loc.size()>0,"input array too small");
      x0=x1=loc[0].x();
      y0=y1=loc[0].y();
      for (size_t i=1; i<loc.size(); ++i)
        {
        x0=std::min(x0,loc[i].x()); x1=std::max(x1,loc[i].x());
        y0=std::min(y0,loc[i].y()); y1=std::max(y1,loc[i].y());
        }
      if (x0==x1) x1+=1e-9;
      if (y0==y1) y1+=1e-9;
      idx=nx/(x1-x0);
      idy=ny/(y1-y0);
      for (size_t i=0; i<loc.size(); ++i)
        data[index(loc[i])].push_back(i);
      }
    void add (const vec2 &pos)
      {
      loc.push_back(pos);
      data[index(pos)].push_back(loc.size()-1);
      }
    /*! Returns the indices of all \a loc entries that lie within a circle of
        radius \a rad around \a center. */
    std::vector<size_t> query(const vec2 &center, double rad) const
      {
      std::vector<size_t> res;
      if ((center.x()<x0-rad)||(center.x()>x1+rad)
        ||(center.y()<y0-rad)||(center.y()>y1+rad))
        return res;
      double rsq=rad*rad;
      size_t i0=indexx(center.x()-rad), i1=indexx(center.x()+rad),
            j0=indexy(center.y()-rad), j1=indexy(center.y()+rad);
      for (size_t j=j0; j<=j1; ++j)
        for (size_t i=i0; i<=i1; ++i)
          for (auto k : data[i+nx*j])
            if (std::norm(center-loc[k])<=rsq) res.push_back(k);
      return res;
      }
    bool anyIn (const vec2 &center, double rad) const
      {
      if ((center.x()<x0-rad)||(center.x()>x1+rad)
        ||(center.y()<y0-rad)||(center.y()>y1+rad))
        return false;
      double rsq=rad*rad;
      size_t i0=indexx(center.x()-rad), i1=indexx(center.x()+rad),
             j0=indexy(center.y()-rad), j1=indexy(center.y()+rad);
      for (size_t j=j0; j<=j1; ++j)
        for (size_t i=i0; i<=i1; ++i)
          for (auto k : data[i+nx*j])
            if (std::norm(center-loc[k])<=rsq) return true;
      return false;
      }
  };

/*! Remove a given value from a vector of integers. Assert that exactly one
    value was removed. */
inline void stripout (std::vector<size_t> &v, size_t val)
  {
  size_t s1=v.size();
  v.erase(std::remove(v.begin(),v.end(),val));
  planck_assert(v.size()+1==s1,"oops");
  }

fpraster tgt2raster (const vector<Target> &tgt, int nx, int ny)
  {
  vector<vec2> tpos(tgt.size());
  for (size_t i=0; i<tgt.size(); ++i) tpos[i]=tgt[i].pos;
  return fpraster (tpos,nx,ny);
  }

fpraster cbr2raster (const vector<Cobra> &cbr, int nx, int ny)
  {
  vector<vec2> cpos(cbr.size());
  for (size_t i=0; i<cbr.size(); ++i) cpos[i]=cbr[i].center;
  return fpraster (cpos,nx,ny);
  }

constexpr double r_kernel=4.75; // radius of the priority function kernel
inline double kernelfunc(double rsq)
  {
  // simple parabola - quick but probably not optimal
  return std::max(0.,r_kernel*r_kernel-rsq);
//  return sqrt(max(0.,r_kernel*r_kernel-rsq)); // linear decrease
//  return exp(-9*rsq/(r_kernel*r_kernel)); // Gaussian kernel
  }

class ETS_data
  {
  protected:
    const std::vector<Target> &tgt;
    const std::vector<Cobra> &cbr;
    fpraster rtgt, rcbr;
    std::vector<std::vector<size_t>> f2t,t2f;
    double colldist, rmax;

    /*! Computes the fiber->target and target->fiber mappings. */
    void calcMappings()
      {
      f2t=std::vector<std::vector<size_t>>(cbr.size());
      for (size_t i=0; i<cbr.size(); ++i)
        {
        const auto &c(cbr[i]);
        std::vector<size_t> tmp=rtgt.query(c.center,c.l1+c.l2);
        for (auto j : tmp)
          if ((std::norm(c.dotpos-tgt[j].pos)>=c.rdot*c.rdot)
            &&(std::norm(c.center-tgt[j].pos)>=(c.l1-c.l2)*(c.l1-c.l2)))
            f2t[i].push_back(j);
        }
      t2f=std::vector<std::vector<size_t>>(tgt.size());
      for (size_t i=0; i<f2t.size(); ++i)
        for (auto t : f2t[i])
          t2f[t].push_back(i);
      }

    int maxpri_in_fiber (size_t fiber) const
      {
      using namespace std;
      constexpr auto randseed=42;
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
      vec2 fpos=cbr[fiber].center;
      int maxpri = tgt[f2t[fiber][0]].pri;
      double mindsq=std::norm(fpos-tgt[f2t[fiber][0]].pos);
      size_t idx=0;
      for (size_t j=1; j<f2t[fiber].size(); ++j)
        {
        if (tgt[f2t[fiber][j]].pri==maxpri)
          {
          if (std::norm(fpos-tgt[f2t[fiber][j]].pos)<mindsq)
            { idx=j; mindsq=std::norm(fpos-tgt[f2t[fiber][j]].pos); }
          }
        else if (tgt[f2t[fiber][j]].pri<maxpri)
          {
          maxpri=tgt[f2t[fiber][j]].pri;
          idx=j;
          mindsq=std::norm(fpos-tgt[f2t[fiber][j]].pos);
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
             elbowpos(elbow_pos(cbr[fiber], tippos));
        std::vector<size_t> tmp=rtgt.query(0.5*(tippos+elbowpos),
             colldist+0.5*cbr[fiber].l2);
        for (auto i : tmp)
          if (line_segment_collision (elbowpos, tippos, tgt[i].pos, colldist))
            {
            for (auto j : t2f[i]) stripout(f2t[j],i);
            t2f[i].clear();
            }
        // remove all other target-fiber combinations that would leave the elbow
        // within the lower arm region of this cobra.
        auto tmp2 = rcbr.query(0.5*(tippos+elbowpos),
                    colldist+rmax+0.5*cbr[fiber].l2);//FIXME!!
        for (auto i:tmp2)
          {
          auto cpy(f2t[i]);
          for (auto j:cpy)
            {
            vec2 tippos2 (tgt[j].pos),
                 elbowpos2(elbow_pos(cbr[i], tippos2));
            if (line_segment_collision (elbowpos, tippos, elbowpos2, colldist))
              {
              stripout(f2t[i],j);
              stripout(t2f[j],i);
              }
            }
          }
        }
    //  checkMappings(tgt,f2t,t2f);
      }
  public:
    ETS_data(const std::vector<Target> &tgt_, const std::vector<Cobra> &cbr_,
      double colldist_=2.)
      : tgt(tgt_), cbr(cbr_), rtgt(tgt2raster(tgt,100,100)),
        rcbr(cbr2raster(cbr,100,100)),colldist(colldist_)
      {
      calcMappings();
      rmax=0;
      for (auto c : cbr)
        rmax=max(rmax,c.l1+c.l2);
      }
    std::vector<std::vector<size_t>> F2T() const { return f2t; }
    std::vector<std::vector<size_t>> T2F() const { return t2f; }
  };

inline void rotate (vec2 &pos, double sa, double ca)
  {
  pos*=std::complex<double>(ca,sa);
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

/*! Naive assignment algorithm: iterate over all fibers, and if a fiber has
    targets in its patrol area, assign the target with the highest priority
    to it. */
class NaiveAssigner: public ETS_data
  {
  public:
    NaiveAssigner (const vector<Target> &tgt_, const std::vector<Cobra> &cobras_,
      vector<size_t> &tid, vector<size_t> &cid)
      : ETS_data(tgt_, cobras_)
      {
      tid.clear(); cid.clear();
      for (size_t fiber=0; fiber<f2t.size(); ++fiber)
        {
        if (f2t[fiber].empty()) continue;
        int itgt = maxpri_in_fiber(fiber);
        tid.push_back(itgt);
        cid.push_back(fiber);
        cleanup (fiber, itgt);
        }
      }
  };

  /*! Assignment strategy modeled after Morales et al. 2012: MNRAS 419, 1187
      find the fiber(s) with the smallest number of observable targets >0;
      for the first of the returned fibers, assign the target with highest
      priority to it; repeat until no more targets are observable. */
class DrainingAssigner: public ETS_data
  {
  public:
    DrainingAssigner (const vector<Target> &tgt_,
      const std::vector<Cobra> &cobras_, vector<size_t> &tid,
      vector<size_t> &cid)
      : ETS_data(tgt_, cobras_)
      {
      tid.clear(); cid.clear();
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
        cid.push_back(fiber);
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
class DrainingClosestAssigner: public ETS_data
  {
  public:
    DrainingClosestAssigner (const vector<Target> &tgt_,
      const std::vector<Cobra> &cobras_,
      vector<size_t> &tid, vector<size_t> &cid)
      : ETS_data(tgt_, cobras_)
      {
      tid.clear(); cid.clear();
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
        cid.push_back(fiber);
        cleanup(fiber,itgt);
        }
      }
  };

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


  /*! Assignment strategy with the goal of reducing inhomogeneity in the
      target distribution: assign a priority to each target that depends on
      the distance of all other targets in its close vicinity; process targets
      in order of decreasing priority and assign them to fibers, if possible.
      After each assignment, update the priority of the remaining targets. */
class NewAssigner: public ETS_data
  {
  private:
    pqueue<pq_entry> calc_pri() const
      {
      std::vector<pq_entry> pri(tgt.size());
      for (size_t i=0; i<tgt.size(); ++i)
        {
        if (t2f[i].size()>0)
          {
          std::vector<size_t> ngb = rtgt.query(tgt[i].pos,r_kernel);

          for (auto j : ngb)
            {
            if (i==j)
              pri[i].prox+=tgt[i].time*tgt[i].time*kernelfunc(0.);
            if (i<j)
              {
              double tmp=tgt[i].time*tgt[j].time
                        *kernelfunc(std::norm(tgt[i].pos-tgt[j].pos));
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
      std::vector<size_t> ngb = rtgt.query(tgt[itgt].pos,r_kernel);
      for (auto j : ngb)
        if ((!t2f[j].empty())||(pri.priority(j).prox!=0.))
          {
          pq_entry tpri=pri.priority(j);
          tpri.prox-=tgt[j].time*tgt[itgt].time
                    *kernelfunc(std::norm(tgt[itgt].pos-tgt[j].pos));
          pri.set_priority(tpri,j);
          }
      }

  public:
    NewAssigner (const vector<Target> &tgt_, const std::vector<Cobra> &cobras_,
      vector<size_t> &tid, vector<size_t> &cid)
      : ETS_data(tgt_, cobras_)
      {
      tid.clear(); cid.clear();
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
        cid.push_back(fiber);
        cleanup(fiber,itgt);
        fix_priority(itgt,pri);
        }
      }
  };

} // unnamed namespace

std::vector<Cobra> makeCobras()
  {
  constexpr size_t nfiber=3*57*14;
  constexpr double armlen=2.375; // length of upper/lower arms
  constexpr double dotdist=1.375; // radius of the dot blocking area
  std::vector<Cobra> res;
  for (size_t i=0; i<nfiber; ++i)
    res.emplace_back(id2fiberpos(i),armlen,armlen,id2dotpos(i),dotdist);
  return res;
  }

std::vector<std::vector<size_t>> getT2F (const std::vector<Target> &targets,
  const std::vector<Cobra> &cobras)
  {
  ETS_data tmp(targets,cobras);
  return tmp.T2F();
  }

void getObservation (const std::vector<Target> &targets,
  const std::vector<Cobra> &cobras, const std::string &algorithm,
  std::vector<size_t> &tid, std::vector<size_t> &cid)
  {
  if (algorithm=="naive")
    NaiveAssigner dummy(targets,cobras,tid,cid);
  else if (algorithm=="draining")
    DrainingAssigner dummy(targets,cobras,tid,cid);
  else if (algorithm=="draining_closest")
    DrainingClosestAssigner dummy(targets,cobras,tid,cid);
  else if (algorithm=="new")
    DrainingClosestAssigner dummy(targets,cobras,tid,cid);
  else
    planck_fail("unknown assignment algorithm");
  }

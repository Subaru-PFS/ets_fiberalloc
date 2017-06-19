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

#ifndef ETS_H
#define ETS_H

#include <complex>
#include <vector>
#include <algorithm>
#include <string>

/*! Simple class for storing a position in a 2D plane. */
class vec2: public std::complex<double>
  {
  public:
    vec2() {}
    vec2(double x_, double y_) : std::complex<double>(x_,y_) {}
    vec2(const std::complex<double> &other) : std::complex<double>(other) {}
    double x() const { return real(); };
    double y() const { return imag(); };
  };

/*! Simple class containing all relevant properties of a PFS observation
    target. */
class Target
  {
  public:
    vec2 pos;
    double time;
    int pri;

    Target (const vec2 &pos_,double time_,int pri_)
      : pos(pos_), time(time_), pri(pri_) {}
  };

class Cobra
  {
  public:
    vec2 center;
    double l1,l2;
    vec2 dotpos;
    double rdot;

    Cobra(const vec2 &center_, double l1_, double l2_, const vec2 &dotpos_,
      double rdot_)
      : center(center_), l1(l1_), l2(l2_), dotpos(dotpos_), rdot(rdot_) {}
  };


std::vector<std::vector<size_t>> getT2F (const std::vector<Target> &targets,
  const std::vector<Cobra> &cobras);

void getObservation (const std::vector<Target> &targets,
  const std::vector<Cobra> &cobras, const std::string &algorithm,
  std::vector<size_t> &tid, std::vector<size_t> &cid);

#endif

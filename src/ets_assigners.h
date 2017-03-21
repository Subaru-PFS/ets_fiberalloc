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

/*! \file ets_assigners.h
 *  Copyright (C) 2017 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef ETS_ASSIGNERS_H
#define ETS_ASSIGNERS_H

#include <vector>
#include <memory>
#include "ets_tools.h"

class FiberAssigner
  {
  public:
    /*! Assign target from \a tgt to fibers according to a given strategy.
        On return, \a tid and \a fid contain target resp. fiber IDs for the
        assigned targets. Target IDs range from 0 to \a tgt.size()-1, fiber
        IDs from 0 to 2393. */
    virtual void assign (const std::vector<Target> &tgt,
      std::vector<size_t> &tid, std::vector<size_t> &fid) const = 0;

    virtual ~FiberAssigner() {}
  };

std::unique_ptr<FiberAssigner> make_assigner(const std::string &name);

#endif

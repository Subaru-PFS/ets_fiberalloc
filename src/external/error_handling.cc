/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Utilities for error reporting
 *
 *  Copyright (C) 2003-2016 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include <cstdlib>

#include "error_handling.h"
#ifdef USE_MPI
#include "mpi_support.h"
#endif

#ifdef LEVELS_STACKTRACE

#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {

void print_trace()
  {
  char pid_buf[30];
  sprintf(pid_buf, "%d", getpid());
  char name_buf[512];
  name_buf[readlink("/proc/self/exe", name_buf, 511)]=0;
  int child_pid = fork();
  if (!child_pid) {           
    dup2(2,1); // redirect output to stderr
    fprintf(stdout,"stack trace for %s pid=%s\n",name_buf,pid_buf);
    execlp("gdb", "gdb", "--batch", "-n", "-ex", "thread", "-ex", "bt",
      name_buf, pid_buf, NULL);
    abort(); /* If gdb failed to start */
    }
  else
    {
    waitpid(child_pid,NULL,0);
    }
  }

}

#endif


using namespace std;

PlanckError::PlanckError(const string &message) : msg (message) {}
PlanckError::PlanckError(const char *message) : msg (message) {}

//virtual
PlanckError::~PlanckError() {}

void planck_failure__(const char *file, int line, const char *func,
  const string &msg)
  {
  cerr << "Error encountered at " << file << ", line " << line << endl;
  if (func) cerr << "(function " << func << ")" << endl;
  if (msg!="") cerr << endl << msg << endl;
  cerr << endl;
#ifdef LEVELS_STACKTRACE
  print_trace();
#endif
  }

void planck_failure__(const char *file, int line, const char *func,
  const char *msg)
  { planck_failure__ (file,line,func,string(msg)); }

void killjob__()
  {
#ifdef USE_MPI
  mpiMgr.abort();
#else
  exit(1);
#endif
  }

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Christina Payne (Vanderbilt U)
                        Stan Moore (Sandia) for dipole terms
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include "fix_bfield.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "region.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixBfield::FixBfield(LAMMPS *lmp, int narg, char **arg) :
		          Fix(lmp, narg, arg), xstr(NULL), ystr(NULL), zstr(NULL),
		          estr(NULL), idregion(NULL), bfield(NULL)
{
  if (narg < 6) error->all(FLERR,"Illegal fix bfield command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  scalar_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  extscalar = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
  virial_flag = 1;

  qe2f = force->qe2f;
  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else {
    bx = qe2f * force->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }

  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else {
    by = qe2f * force->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }

  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else {
    bz = qe2f * force->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }

  // optional args

  iregion = -1;
  idregion = NULL;
  estr = NULL;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bfield command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix bfield does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bfield command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        estr = new char[n];
        strcpy(estr,&arg[iarg+1][2]);
      } else error->all(FLERR,"Illegal fix bfield command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix bfield command");
  }

  force_flag = 0;
  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;

  maxatom = atom->nmax;
  memory->create(bfield,maxatom,4,"bfield:bfield");
}

/* ---------------------------------------------------------------------- */

FixBfield::~FixBfield()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] estr;
  delete [] idregion;
  memory->destroy(bfield);
}

/* ---------------------------------------------------------------------- */

int FixBfield::setmask()
{
  int mask = 0;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBfield::init()
{
  mumflag = 0;
  if (atom->mum_flag && atom->torque_flag) mumflag = 1;
  if (!mumflag)
    error->all(FLERR,"Fix bfield requires atom attribute q or mu");

  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix bfield does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix bfield is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix bfield does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix bfield is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix bfield does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix bfield is invalid style");
  }
  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0)
      error->all(FLERR,"Variable name for fix bfield does not exist");
    if (input->variable->atomstyle(evar)) estyle = ATOM;
    else error->all(FLERR,"Variable for fix bfield is invalid style");
  } else estyle = NONE;


  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix aveforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (mumflag && varflag == ATOM)
    error->all(FLERR,"Fix bfield with dipoles cannot use atom-style variables");

  if (mumflag && update->whichflag == 2 && comm->me == 0)
    error->warning(FLERR,
        "The minimizer does not re-orient dipoles "
        "when using fix bfield");

  if (varflag == CONSTANT && estyle != NONE)
    error->all(FLERR,"Cannot use variable energy with "
        "constant bfield in fix bfield");
  if ((varflag == EQUAL || varflag == ATOM) &&
      update->whichflag == 2 && estyle == NONE)
    error->all(FLERR,"Must use variable energy with fix bfield");

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  //Full Neighbor List
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixBfield::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBfield::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixBfield::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   apply F = qE
------------------------------------------------------------------------- */

void FixBfield::post_force(int vflag)
{
  int i,j,k,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double dx, dy, dz, fx, fy, fz, rsq, r, mir, mjr, mumu, A, K, muR, C, mui;
  double *rad = atom->radius;
  double **x = atom->x;
  double **mum = atom->mum;
  double **omega = atom->omega;
  double **f = atom->f;
  double **torque = atom->torque;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  bigint ntimestep = update->ntimestep;
  double p4 = M_PI*4;
  double u = p4*1e-7;

  // energy and virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // reallocate bfield array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(bfield);
    memory->create(bfield,maxatom,4,"bfield:bfield");
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // fsum[0] = "potential energy" for added force
  // fsum[123] = extra force added to atoms

  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  force_flag = 0;


  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (mask[i] & groupbit) {

      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj<jnum; jj++)  {

        j =jlist[jj];
        j &= NEIGHMASK;
        dx = x[i][0] - x[j][0];
        dy = x[i][1] - x[j][1];
        dz = x[i][2] - x[j][2];
        rsq = dx*dx + dy*dy + dz*dz;
        r = sqrt(rsq);
        //std::cout<<"r = "<<r<<" ";

        K = 3e-7/rsq/rsq;

        dx /= r;
        dy /= r;
        dz /= r;

        mir = mum[i][0]*dx+mum[i][1]*dy+mum[i][2]*dz;
        mjr = mum[j][0]*dx+mum[j][1]*dy+mum[j][2]*dz;
        mumu = mum[i][0]*mum[j][0]+mum[i][1]*mum[j][1]+mum[i][2]*mum[j][2];
        //std::cout<<"fm = "<<fx<<","<<fy<<","<<fz<<" ";

        f[i][0] += K*(mir*mum[j][0]+mjr*mum[i][0]+(mumu-5*mjr*mir)*dx);
        f[i][1] += K*(mir*mum[j][1]+mjr*mum[i][1]+(mumu-5*mjr*mir)*dy);
        f[i][2] += K*(mir*mum[j][2]+mjr*mum[i][2]+(mumu-5*mjr*mir)*dz);

        A = u/p4/r/rsq;
        torque[i][0] += mum[i][1]*A*(3*mjr*dz-mum[j][2])-mum[i][2]*A*(3*mjr*dy-mum[j][1]);
        torque[i][1] += mum[i][2]*A*(3*mjr*dx-mum[j][0])-mum[i][0]*A*(3*mjr*dz-mum[j][2]);
        torque[i][2] += mum[i][0]*A*(3*mjr*dy-mum[j][1])-mum[i][1]*A*(3*mjr*dx-mum[j][0]);

      }
    }
  }

    }

    /* ---------------------------------------------------------------------- */

    void FixBfield::post_force_respa(int vflag, int ilevel, int /*iloop*/)
    {
      if (ilevel == ilevel_respa) post_force(vflag);
    }

    /* ---------------------------------------------------------------------- */

    void FixBfield::min_post_force(int vflag)
    {
      post_force(vflag);
    }

    /* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

    double FixBfield::memory_usage()
    {
      double bytes = 0.0;
      if (varflag == ATOM) bytes = atom->nmax*4 * sizeof(double);
      return bytes;
    }

    /* ----------------------------------------------------------------------
   return energy added by fix
------------------------------------------------------------------------- */

    double FixBfield::compute_scalar(void)
    {
      if (force_flag == 0) {
        MPI_Allreduce(fsum,fsum_all,4,MPI_DOUBLE,MPI_SUM,world);
        force_flag = 1;
      }
      return fsum_all[0];
    }

    /* ----------------------------------------------------------------------
   return total extra force due to fix
------------------------------------------------------------------------- */

    double FixBfield::compute_vector(int n)
    {
      if (force_flag == 0) {
        MPI_Allreduce(fsum,fsum_all,4,MPI_DOUBLE,MPI_SUM,world);
        force_flag = 1;
      }
      return fsum_all[n+1];
    }

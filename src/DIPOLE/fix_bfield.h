/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(bfield,FixBfield)

#else

#ifndef LMP_FIX_BFIELD_H
#define LMP_FIX_BFIELD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBfield : public Fix {
 public:
  FixBfield(class LAMMPS *, int, char **);
  ~FixBfield();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  //void init_list(int, class NeighList *);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double memory_usage();
  double compute_scalar();
  double compute_vector(int);

 private:
  double bx,by,bz;
  int varflag,iregion;
  char *xstr,*ystr,*zstr,*estr,*xpstr,*ypstr,*zpstr;
  char *idregion;
  int xvar,yvar,zvar,evar,xstyle,ystyle,zstyle,estyle;
  int ilevel_respa;
  double qe2f;
  int muflag;
  class NeighList *list;

  int maxatom;
  double **bfield;

  int force_flag;
  double fsum[4],fsum_all[4];
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix bfield does not exist

Self-explanatory.

E: Fix bfield requires atom attribute q or mu

The atom style defined does not have this attribute.

E: Variable name for fix bfield does not exist

Self-explanatory.

E: Variable for fix bfield is invalid style

The variable must be an equal- or atom-style variable.

E: Region ID for fix aveforce does not exist

Self-explanatory.

E: Fix bfield with dipoles cannot use atom-style variables

This option is not supported.

W: The minimizer does not re-orient dipoles when using fix bfield

This means that only the atom coordinates will be minimized,
not the orientation of the dipoles.

E: Cannot use variable energy with constant bfield in fix bfield

LAMMPS computes the energy itself when the B-field is constant.

E: Must use variable energy with fix bfield

You must define an energy when performing a minimization with a
variable B-field.

*/

// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ScalarExpr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-31
/// @brief    Implementation of the ScalarExpr class

#include <cstring>

#include "problem.hpp"

//----------------------------------------------------------------------

ScalarExpr::ScalarExpr
(Param * param) throw()
  : param_(param)
{
  printf ("DEBUG %s:%d %p\n",__FILE__,__LINE__,param_);
}

//----------------------------------------------------------------------

void ScalarExpr::copy_(const ScalarExpr & scalar_expr) throw()
{
  printf ("DEBUG %s:%d %p\n",__FILE__,__LINE__,param_);
  param_ = new Param(*scalar_expr.param_);
}

//----------------------------------------------------------------------

double ScalarExpr::evaluate (double t, double x, double y, double z, 
			     Mask * mask, double deflt) const
{
  double value;
  printf ("DEBUG %s:%d %p\n",__FILE__,__LINE__,param_);
  bool m = mask ? mask->evaluate(t,x,y,z) : true;
  if (m) {
    param_->evaluate_float(1,&value,&x,&y,&z,t);
  } else {
    value = deflt;
  }
  return value;
}

//----------------------------------------------------------------------

void ScalarExpr::evaluate (double * value, double t,
			   int ndx, int nx, double * xv,
			   int ndy, int ny, double * yv,
			   int ndz, int nz, double * zv, 
			   Mask * mask, double * deflt) const
{
  ASSERT6("ScalarExpr::evaluate",
	  "value dimension (%d %d %d) needs to be at least (%d %d %d)",
	  ndx,ndy,ndz,nx,ny,nz,
	  (ndx >= nx) && (ndy >= ny) && (ndz >= nz));

  double * x = new double [ nx*ny*nz ];
  double * y = new double [ nx*ny*nz ];
  double * z = new double [ nx*ny*nz ];
  double * value_temp = new double [nx*ny*nz];
  bool * mv = 0;
  if (mask) {
    mv = new bool [ nx*ny*nz ];
    mask->evaluate(mv, t, ndx,nx,xv, ndy,ny,yv, ndz,nz,zv);
  }
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + nx*(iy + ny*iz);
	x[i] = xv[ix];
	y[i] = yv[iy];
	z[i] = zv[iz];
	value_temp[i] = 0.0;
      }
    }
  }

  int n=nx*ny*nz;

  param_->evaluate_float(n, value_temp, x,y,z,t);

  if (mask) {
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  int i=ix + nx*(iy + ny*iz);
	  int id=ix + ndx*(iy + ndy*iz);
	  value[id] = mv[i] ? value_temp[i] : deflt[id];
	}
      }
    }
  } else {
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	for (int iz=0; iz<nz; iz++) {
	  int i=ix + nx*(iy + ny*iz);
	  int id=ix + ndx*(iy + ndy*iz);
	  value[id] = (double) value_temp[i];
	}
      }
    }
  }

  delete [] value_temp;
  delete [] z;
  delete [] y;
  delete [] x;

}
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPmDeposit.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPmDeposit class
///
/// The EnzoMethodPmDeposit method computes a "density_total" field,
/// which includes the "density" field plus mass from gravitating
/// particles (particles in the "mass" group, e.g. "dark" matter
/// particles)

#include "cello.hpp"
#include "enzo.hpp"

// #define DEBUG_PM_DEPOSIT
// #define DEBUG_FIELD

#ifdef DEBUG_PM_DEPOSIT
#  define TRACE_PM(MESSAGE)				\
  CkPrintf ("%s:%d %s\n", __FILE__,__LINE__,MESSAGE);				
#else
#  define TRACE_PM(MESSAGE) /* ... */
#endif

#define FORTRAN_NAME(NAME) NAME##_

extern "C" void  FORTRAN_NAME(dep_grid_cic)
  (enzo_float * de,enzo_float * de_t,enzo_float * temp,
   enzo_float * vx, enzo_float * vy, enzo_float * vz, 
   enzo_float * dt, enzo_float * rfield, int *rank,
   enzo_float * hx, enzo_float * hy, enzo_float * hz,
   int * mx,int * my,int * mz,
   int * gxi,int * gyi,int * gzi,
   int * nxi,int * nyi,int * nzi,
   int * ,int * ,int * ,
   int * nx,int * ny,int * nz,
   int * ,int * ,int * );

//----------------------------------------------------------------------

EnzoMethodPmDeposit::EnzoMethodPmDeposit 
(const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 std::string type,
 double alpha)
  : Method(),
    alpha_(alpha),
    type_(pm_type_unknown)
{
  TRACE_PM("EnzoMethodPmDeposit()");
  if (type == "cic") {
    type_ = pm_type_cic;
  } else if (type == "ngp") {
    type_ = pm_type_ngp;
    WARNING("EnzoMethodPmDeposit()",
	    "PM deposit type 'ngp' is not implemented yet: using 'cic'");
  } else if (type == "tsc") {
    WARNING("EnzoMethodPmDeposit()",
	    "PM deposit type 'tsc' is not implemented yet: using 'cic'");
    type_ = pm_type_tsc;
  } else {
    ERROR1 ("EnzoMethodPmDeposit()",
	    "PM deposit type '%s' is not supported: "
	    "must be 'cic','ngp', or 'tsc'",
	    type.c_str());
  }
  // Initialize default Refresh object

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_fields();

  // PM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPmDeposit::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | alpha_;
  p | type_;
}

//----------------------------------------------------------------------

void EnzoMethodPmDeposit::compute ( Block * block) throw()
{
  TRACE_PM("compute()");

  if (block->is_leaf()) {

    Particle particle (block->data()->particle());
    Field    field    (block->data()->field());

    int rank = block->rank();

    enzo_float  * de_t = (enzo_float *) field.values("density_total");

    int mx,my,mz;
    field.dimensions(0,&mx,&my,&mz);
    int nx,ny,nz;
    field.size(&nx,&ny,&nz);
    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);

    // Initialize "density_total" with gas "density"

    for (int i=0; i<mx*my*mz; i++) de_t[i] = 0;
    // Get block extents and cell widths
    double xm,ym,zm;
    double xp,yp,zp;
    block->lower(&xm,&ym,&zm);
    block->upper(&xp,&yp,&zp);

    double hx,hy,hz;
    hx = (xp-xm)/nx;
    hy = (yp-ym)/ny;
    hz = (zp-zm)/nz;

    // declare particle position arrays

    const int it = particle.type_index ("dark");

    const int ia_mass = particle.constant_index (it,"mass");

    enzo_float dens = *((enzo_float *)(particle.constant_value (it,ia_mass)));

    // Scale by volume if particle value is mass instead of density
    
    // double vol = 1.0;
    // if (rank >= 1) vol *= hx;
    // if (rank >= 2) vol *= hy;
    // if (rank >= 3) vol *= hz;
    //    dens = dens / vol;

    // check precisions match
    
    int ia = particle.attribute_index(it,"x");
    int ba = particle.attribute_bytes(it,ia); // "bytes (actual)"
    int be = sizeof(enzo_float);                // "bytes (expected)"

    ASSERT4 ("EnzoMethodPmUpdate::compute()",
	     "Particle type %s attribute %s defined as %s but expecting %s",
	     particle.type_name(it).c_str(),
	     particle.attribute_name(it,ia).c_str(),
	     ((ba == 4) ? "single" :
	      ((ba == 8) ? "double" : "quadruple")),
	     ((be == 4) ? "single" :
	      ((be == 8) ? "double" : "quadruple")),
	     (ba == be));

    // Accumulate particle density using CIC

    const double dt = alpha_ * block->dt();

    // Accumulated single velocity array for Baryon deposit

    for (int ib=0; ib<particle.num_batches(it); ib++) {

      const int np = particle.num_particles(it,ib);

      if (rank == 1) {

	const int ia_x  = particle.attribute_index(it,"x");
	const int ia_vx = particle.attribute_index(it,"vx");

	enzo_float * xa =  (enzo_float *) particle.attribute_array (it,ia_x,ib);
	enzo_float * vxa = (enzo_float *) particle.attribute_array (it,ia_vx,ib);

	const int dp =  particle.stride(it,ia_x);
	const int dv =  particle.stride(it,ia_vx);

	for (int ip=0; ip<np; ip++) {

	  double x = xa[ip*dp] + vxa[ip*dv]*dt;

	  double tx = nx*(x - xm) / (xp - xm) - 0.5;

	  int ix0 = gx + floor(tx);

	  int ix1 = ix0 + 1;

	  double x0 = 1.0 - (tx - floor(tx));
	  double x1 = 1.0 - x0;

	  de_t[ix0] += dens*x0;
	  de_t[ix1] += dens*x1;

	}

      } else if (rank == 2) {

	const int ia_x  = particle.attribute_index(it,"x");
	const int ia_y  = particle.attribute_index(it,"y");
	const int ia_vx = particle.attribute_index(it,"vx");
	const int ia_vy = particle.attribute_index(it,"vy");

	// Batch arrays
	enzo_float * xa  = (enzo_float *) particle.attribute_array (it,ia_x,ib);
	enzo_float * ya  = (enzo_float *) particle.attribute_array (it,ia_y,ib);

	// Particle batch velocities
	enzo_float * vxa = (enzo_float *) particle.attribute_array (it,ia_vx,ib);
	enzo_float * vya = (enzo_float *) particle.attribute_array (it,ia_vy,ib);

	const int dp =  particle.stride(it,ia_x);
	const int dv =  particle.stride(it,ia_vx);

	for (int ip=0; ip<np; ip++) {

	  double x = xa[ip*dp] + vxa[ip*dv]*dt;
	  double y = ya[ip*dp] + vya[ip*dv]*dt;

	  double tx = nx*(x - xm) / (xp - xm) - 0.5;
	  double ty = ny*(y - ym) / (yp - ym) - 0.5;

	  int ix0 = gx + floor(tx);
	  int iy0 = gy + floor(ty);

	  int ix1 = ix0 + 1;
	  int iy1 = iy0 + 1;

	  double x0 = 1.0 - (tx - floor(tx));
	  double y0 = 1.0 - (ty - floor(ty));

	  double x1 = 1.0 - x0;
	  double y1 = 1.0 - y0;

	  if ( dens < 0.0) {
	    CkPrintf ("%s:%d ERROR: dens = %f\n", __FILE__,__LINE__,dens);
	  }

	  de_t[ix0+mx*iy0] += dens*x0*y0;
	  de_t[ix1+mx*iy0] += dens*x1*y0;
	  de_t[ix0+mx*iy1] += dens*x0*y1;
	  de_t[ix1+mx*iy1] += dens*x1*y1;

	  if ( de_t[ix0+mx*iy0] < 0.0) {
	    CkPrintf ("%s:%d ERROR: de_t %d %d = %f\n",
		      __FILE__,__LINE__,ix0,iy0,dens);
	  }
	  if ( de_t[ix1+mx*iy0] < 0.0) {
	    CkPrintf ("%s:%d ERROR: de_t %d %d = %f\n",
		      __FILE__,__LINE__,ix1,iy0,dens);
	  }
	  if ( de_t[ix0+mx*iy1] < 0.0) {
	    CkPrintf ("%s:%d ERROR: de_t %d %d = %f\n",
		      __FILE__,__LINE__,ix0,iy1,dens);
	  }
	  if ( de_t[ix1+mx*iy1] < 0.0) {
	    CkPrintf ("%s:%d ERROR: de_t %d %d = %f\n",
		      __FILE__,__LINE__,ix1,iy1,dens);
	  }
	}

      } else if (rank == 3) {

	const int ia_x  = particle.attribute_index(it,"x");
	const int ia_y  = particle.attribute_index(it,"y");
	const int ia_z  = particle.attribute_index(it,"z");
	const int ia_vx = particle.attribute_index(it,"vx");
	const int ia_vy = particle.attribute_index(it,"vy");
	const int ia_vz = particle.attribute_index(it,"vz");

	enzo_float * xa  = (enzo_float *) particle.attribute_array (it,ia_x,ib);
	enzo_float * ya  = (enzo_float *) particle.attribute_array (it,ia_y,ib);
	enzo_float * za  = (enzo_float *) particle.attribute_array (it,ia_z,ib);

	// Particle batch velocities
	enzo_float * vxa = (enzo_float *) particle.attribute_array (it,ia_vx,ib);
	enzo_float * vya = (enzo_float *) particle.attribute_array (it,ia_vy,ib);
	enzo_float * vza = (enzo_float *) particle.attribute_array (it,ia_vz,ib);

	const int dp =  particle.stride(it,ia_x);
	const int dv =  particle.stride(it,ia_vx);

	for (int ip=0; ip<np; ip++) {

	  // Copy batch particle velocities to temporary block field velocities
	  
	  double x = xa[ip*dp] + vxa[ip*dv]*dt;
	  double y = ya[ip*dp] + vya[ip*dv]*dt;
	  double z = za[ip*dp] + vza[ip*dv]*dt;

	  double tx = nx*(x - xm) / (xp - xm) - 0.5;
	  double ty = ny*(y - ym) / (yp - ym) - 0.5;
	  double tz = nz*(z - zm) / (zp - zm) - 0.5;

	  int ix0 = gx + floor(tx);
	  int iy0 = gy + floor(ty);
	  int iz0 = gz + floor(tz);

	  int ix1 = ix0 + 1;
	  int iy1 = iy0 + 1;
	  int iz1 = iz0 + 1;

	  double x0 = 1.0 - (tx - floor(tx));
	  double y0 = 1.0 - (ty - floor(ty));
	  double z0 = 1.0 - (tz - floor(tz));

	  double x1 = 1.0 - x0;
	  double y1 = 1.0 - y0;
	  double z1 = 1.0 - z0;

	  de_t[ix0+mx*(iy0+my*iz0)] += dens*x0*y0*z0;
	  de_t[ix1+mx*(iy0+my*iz0)] += dens*x1*y0*z0;
	  de_t[ix0+mx*(iy1+my*iz0)] += dens*x0*y1*z0;
	  de_t[ix1+mx*(iy1+my*iz0)] += dens*x1*y1*z0;
	  de_t[ix0+mx*(iy0+my*iz1)] += dens*x0*y0*z1;
	  de_t[ix1+mx*(iy0+my*iz1)] += dens*x1*y0*z1;
	  de_t[ix0+mx*(iy1+my*iz1)] += dens*x0*y1*z1;
	  de_t[ix1+mx*(iy1+my*iz1)] += dens*x1*y1*z1;

	}
      }
    }

    TRACE_FIELD("density-particle",de_t,1.0);

    enzo_float  * de   = (enzo_float *) field.values("density");
    enzo_float  * de_gas   = (enzo_float *) field.values("density_gas");
    enzo_float * temp = new enzo_float [4*mx*my*mz];
    enzo_float * rfield = new enzo_float[mx*my*mz];
    for (int i=0; i<mx*my*mz; i++) rfield[i] = 0.0;
    for (int i=0; i<4*mx*my*mz; i++) temp[i] = 0.0;

    int gxi=gx;
    int gyi=gy;
    int gzi=gz;
    int nxi=mx-gx-1;
    int nyi=my-gy-1;
    int nzi=mz-gz-1;
    int i0 = 0;
    int i1 = 1;
    enzo_float hxf = hx;
    enzo_float hyf = hy;
    enzo_float hzf = hz;
    enzo_float dtf = 0.0;

#define FIELD_CIC

#ifdef FIELD_CIC    
    enzo_float * de_gas_0 = de_gas + gx + mx*(gy + my*gz);
    for (int i=0; i<mx*my*mz; i++) de_gas[i] = 0.0;
#else
    enzo_float * de_gas_0 = new enzo_float [nx*ny*nz];
    for (int i=0; i<nx*ny*nz; i++) de_gas_0[i] = 0.0;
#endif    
    

    enzo_float * vxf = (enzo_float *) field.values("velocity_x");
    enzo_float * vyf = (enzo_float *) field.values("velocity_y");
    enzo_float * vzf = (enzo_float *) field.values("velocity_z");
    enzo_float * vx = new enzo_float [mx*my*mz];
    enzo_float * vy = new enzo_float [mx*my*mz];
    enzo_float * vz = new enzo_float [mx*my*mz];

    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
	for (int ix=0; ix<mx; ix++) {
	  int k = ix + mx*(iy + my*iz);
	  int i = ix + nx*(iy + ny*iz);
	  vx[k] = vxf[k];
	}
      }
    }

    if (rank >= 2) {
      for (int iz=0; iz<mz; iz++) {
	for (int iy=0; iy<my; iy++) {
	  for (int ix=0; ix<mx; ix++) {
	    int k = ix + mx*(iy + my*iz);
	    int i = ix + nx*(iy + ny*iz);
	    vy[k] = vyf[k];
	  }
	}
      }
    } else {
      for (int i=0; i<mx*my*mz; i++) vy[i] = 0.0;
    }
    
    if (rank >= 3) {
      for (int iz=0; iz<mz; iz++) {
	for (int iy=0; iy<my; iy++) {
	  for (int ix=0; ix<mx; ix++) {
	    int k = ix + mx*(iy + my*iz);
	    int i = ix + nx*(iy + ny*iz);
	    vz[k] = vzf[k];
	  }
	}
      }
    } else {
      for (int i=0; i<mx*my*mz; i++) vz[i] = 0.0;
    }
    
    FORTRAN_NAME(dep_grid_cic)(de,de_gas_0,temp,
			       vx, vy, vz, 
			       &dtf, rfield, &rank,
			       &hxf,&hyf,&hzf,
			       &mx,&my,&mz,
			       &gxi,&gyi,&gzi,
			       &nxi,&nyi,&nzi,
			       &i0,&i0,&i0,
			       &nx,&ny,&nz,
			       &i1,&i1,&i1);

    delete [] rfield;
    delete [] temp;
    delete [] vx;
    delete [] vy;
    delete [] vz;
    TRACE_FIELD("density-gas",de,1.0);
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
    	for (int ix=gx; ix<mx-gx; ix++) {
	  int i = ix + mx*(iy + my*iz);
#ifdef FIELD_CIC
	  int ig = (ix-gx) + nx*((iy-gy) + ny*(iz-gz));
#else
	  int ig = ix + nx*(iy + ny*iz);
#endif	  
    	  de_t[i] += de_gas_0[ig];
    	}
      }
    }

    TRACE_FIELD("density-gas-particle",de_t,1.0);

#ifndef FIELD_CIC
    delete [] de_gas_0;
#endif    

    
  }


  block->compute_done(); 
  
}

//----------------------------------------------------------------------

double EnzoMethodPmDeposit::timestep ( Block * block ) const throw()
{
  TRACE_PM("timestep()");
  double dt = std::numeric_limits<double>::max();

  return dt;
}
// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMethodMultipole.cpp
/// @author   Ryan Golant (ryan.golant@columbia.edu)
/// @date     July 19, 2022
/// @brief    Compute accelerations on gas and particles using FMM


#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/gravity/gravity.hpp"
#include "Enzo/gravity/EnzoEwald.hpp"

//----------------------------------------------------------------------

EnzoMethodMultipole::EnzoMethodMultipole (ParameterGroup p)
  : Method(),
    theta_(p.value_float("theta", 0.0)),
    eps0_(p.value_float("eps0", 0.0)),
    // r0_(p.value_float("r0", 0.0)),
    is_volume_(-1),
    block_volume_(),
    max_volume_(0),
    interp_xpoints_(p.value_integer("interp_xpoints", 32)),
    interp_ypoints_(p.value_integer("interp_ypoints", 32)),
    interp_zpoints_(p.value_integer("interp_zpoints", 32)),
    ewald_order_(p.value_integer("ewald_order", 3)),
    dt_max_(p.value_float("dt_max", 1.0e10)),
    ir_exit_(-1)
{
  // CkPrintf("Have we started?");

  cello::define_field ("density");
  cello::define_field ("acceleration_x");
  cello::define_field ("acceleration_y");
  cello::define_field ("acceleration_z");

  // Initialize default Refresh object
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  // should these be added to the default refresh?
  refresh->add_field("acceleration_x");
  refresh->add_field("acceleration_y");
  refresh->add_field("acceleration_z");

  ir_exit_ = add_refresh_();
  cello::simulation()->refresh_set_name(ir_exit_,name()+":exit");
  Refresh * refresh_exit = cello::refresh(ir_exit_);
  // refresh_exit->set_prolong(index_prolong_);
  refresh_exit->add_field("acceleration_x");
  refresh_exit->add_field("acceleration_y");
  refresh_exit->add_field("acceleration_z");
  refresh_exit->set_callback(CkIndex_EnzoBlock::p_method_multipole_end());


  // Initialize Sync scalars
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  i_sync_restrict_ = scalar_descr_sync->new_value("multipole:restrict");

  // Initialize scalars for receiving messages
  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();
  i_msg_prolong_ = scalar_descr_void->new_value("multipole:msg_prolong");
  for (int ic=0; ic<cello::num_children(); ic++) {
    i_msg_restrict_[ic] = scalar_descr_void->new_value("multipole:msg_restrict");
  }

  // Initialize mass, COM, and quadrupole scalars
  ScalarDescr * scalar_descr_double = cello::scalar_descr_double();
  i_mass_ = scalar_descr_double->new_value("multipole:mass");
  i_com_ = scalar_descr_double->new_value("multipole:com", 3);
  i_quadrupole_ = scalar_descr_double->new_value("multipole:quadrupole", 6);

  // Initialize Taylor coefficient scalars
  i_c1_ = scalar_descr_double->new_value("multipole:c1", 3);
  i_c2_ = scalar_descr_double->new_value("multipole:c2", 6);
  i_c3_ = scalar_descr_double->new_value("multipole:c3", 10);

  
  // Declare long long Block Scalar for volume and save scalar index
  is_volume_ = cello::scalar_descr_long_long()->new_value("solver_fmm_volume");

  ewald_ = nullptr;

  ASSERT1("EnzoMethodMultipole",
          "Method:multipole:ewald_order must be 3 or 4 (it is %d)",
          ewald_order_, (ewald_order_ == 3 || ewald_order_ == 4));

}

//----------------------------------------------------------------------

void EnzoMethodMultipole::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | i_sync_restrict_;
  p | i_msg_prolong_;
  p | ir_exit_;
  PUParray(p, i_msg_restrict_, 8);

  p | i_mass_;
  p | i_com_;
  p | i_quadrupole_;
  p | i_c1_;
  p | i_c2_;
  p | i_c3_;

  p | theta_;
  p | eps0_;
  // p | r0_;
  p | is_volume_;
  p | block_volume_;
  p | max_volume_;
  p | interp_xpoints_;
  p | interp_ypoints_;
  p | interp_zpoints_;
  p | ewald_order_;
  p | dt_max_;

  // bool ewald_is_nullptr = (ewald_ == nullptr);
  // p | ewald_is_nullptr;
  // if (not ewald_is_nullptr) {
  //   if (p.isUnpacking())
  //     ewald_ = new EnzoEwald();
  //   p | *ewald_;
  // }
  
}

//----------------------------------------------------------------------

void EnzoMethodMultipole::compute ( Block * block) throw()
{

  // CkPrintf("have we made it to compute?\n");
  
  Sync * sync_restrict = psync_restrict(block);
  sync_restrict->set_stop(1 + cello::num_children());
  //Sync * sync_prolong = psync_prolong(block);
  //sync_prolong->set_stop(1 + 1);

  double * mass = pmass(block);
  double * com = pcom(block);
  double * quadrupole = pquadrupole(block);
  double * c1 = pc1(block);
  double * c2 = pc2(block);
  double * c3 = pc3(block);
  
  *mass = 0;

  // Is this necessary?
  for (int i = 0; i < 3; i++){
    com[i]  = 0;   
  }
  
  for (int i = 0; i < 6; i++){
    quadrupole[i]  = 0;
  }

  for (int i = 0; i < 3; i++){
      c1[i] = 0;
  }
  
  for (int i = 0; i < 6; i++){
      c2[i] = 0;
    }

  for (int i = 0; i < 10; i++){
      c3[i] = 0;
  }


  // reset cell accelerations to zero
  Field field = block->data()->field();
  
  enzo_float * accel_x = (enzo_float*) field.values ("acceleration_x");
  enzo_float * accel_y = (enzo_float*) field.values ("acceleration_y");
  enzo_float * accel_z = (enzo_float*) field.values ("acceleration_z");

  int mx, my, mz;
  int gx, gy, gz;
  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);

  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {

        int i = ix + mx * (iy + iz * my);
        
        accel_x[i] = 0;
        accel_y[i] = 0;
        accel_z[i] = 0;

      }
    }
  }


  // Allocate and initialize block_volume_[level - min_level] to store
  // weighted volume of blocks in different levels relative to the finest
  // (for the time being, I've set min_level = 0)
  const int num_children = cello::num_children();
  const int max_level_ = cello::hierarchy()->max_level();

  block_volume_.resize(max_level_ + 1);
  block_volume_[max_level_] = 1;

  for (int i = max_level_ - 1; i >= 0; i--) {
    block_volume_[i] = block_volume_[i+1] * num_children;
  }

  // Compute volume of the domain relative to finest-level blocks
  int nb3[3] = {1,1,1};
  cello::hierarchy()->root_blocks(nb3,nb3+1,nb3+2);
  max_volume_ = nb3[0]*nb3[1]*nb3[2];
  max_volume_ *= block_volume_[0];  // the index should probably just be -min_level_ 
                                    

  *volume_(block) = 0;



  // for testing with particles, we initialize particles in leaf Blocks  
  // if (block->is_leaf()) {

  //   int nprtls = 64; // number of particles
  //   double prtls[nprtls][7]; // particle attributes

  //   // attributes: mass, x, y, z, ax, ay, az
  //   // double prtls[nprtls][7] = {{1.0, 0.25, 0.25, 0.25, 0, 0, 0},
  //   //                            {1.0, 0.25, 0.75, 0.25, 0, 0, 0},
  //   //                            {1.0, 0.75, 0.25, 0.25, 0, 0, 0},
  //   //                            {1.0, 0.75, 0.75, 0.25, 0, 0, 0},
  //   //                            {1.0, 0.25, 0.25, 0.75, 0, 0, 0},
  //   //                            {1.0, 0.25, 0.75, 0.75, 0, 0, 0},
  //   //                            {1.0, 0.75, 0.25, 0.75, 0, 0, 0},
  //   //                            {1.0, 0.75, 0.75, 0.75, 0, 0, 0}};
    
  //   for (int i = 0; i < 4; i++) {
  //     for (int j = 0; j < 4; j++) {
  //       for (int k = 0; k < 4; k++) {
  //         int ip = i + 4 * (j + 4 * k);
  //         prtls[ip][0] = 1.0; 
  //         prtls[ip][1] = 0.125 + i * 0.25;
  //         prtls[ip][2] = 0.125 + j * 0.25;
  //         prtls[ip][3] = 0.125 + k * 0.25;
  //         prtls[ip][4] = 0;
  //         prtls[ip][5] = 0;
  //         prtls[ip][6] = 0;
  //       }
  //     }
  //   }

  //   // for (int i = 0; i < 2; i++) {
  //   //   for (int j = 0; j < 2; j++) {
  //   //     for (int k = 0; k < 2; k++) {
  //   //       int ip = i + 2 * (j + 2 * k);
  //   //       prtls[ip][0] = 1.0; 
  //   //       prtls[ip][1] = 0.25 + i * 0.5;
  //   //       prtls[ip][2] = 0.25 + j * 0.5;
  //   //       prtls[ip][3] = 0.25 + k * 0.5;
  //   //       prtls[ip][4] = 0;
  //   //       prtls[ip][5] = 0;
  //   //       prtls[ip][6] = 0;
  //   //     }
  //   //   }
  //   // }
                               
  //   InitializeParticles(block, nprtls, prtls);
  // }

  // reset particle accelerations to zero
  Particle particle = block->data()->particle();
  const int num_prtl_types = particle.num_types();

  for (int it = 0; it < num_prtl_types; it++) {
    
    for (int ib = 0; ib < particle.num_batches(it); ib++) {

      const int np = particle.num_particles(it,ib);

      const int ia_ax  = particle.attribute_index(it,"ax");
      const int ia_ay  = particle.attribute_index(it,"ay");
      const int ia_az  = particle.attribute_index(it,"az");

      enzo_float * axa =  (enzo_float *)particle.attribute_array (it,ia_ax,ib);
      enzo_float * aya =  (enzo_float *)particle.attribute_array (it,ia_ay,ib);
      enzo_float * aza =  (enzo_float *)particle.attribute_array (it,ia_az,ib);

      const int dax =  particle.stride(it,ia_ax);
      const int day =  particle.stride(it,ia_ay);
      const int daz =  particle.stride(it,ia_az);

      for (int ip=0; ip < np; ip++) {
        // CkPrintf("particle accel in compute: %.9f, %.9f, %.9f\n", axa[ip*dax], aya[ip*day], aza[ip*daz]);
        axa[ip*dax] = 0;
        aya[ip*day] = 0;
        aza[ip*daz] = 0;

      }
    }
  }


  // is there an easier way to check for periodicity?
  int is_periodic_x; int is_periodic_y; int is_periodic_z;
  cello::hierarchy()->get_periodicity(&is_periodic_x, &is_periodic_y, &is_periodic_z);

  // CkPrintf("before ewald?\n");

  if ((ewald_ == nullptr) && (is_periodic_x || is_periodic_y || is_periodic_z)) {
    // Call the primary constructor of EnzoEwald that takes the
    // dimensions of the downsampled interpolation grid as input parameters.
    // Then copy (or move if the compiler is smart) the result into this->ewald_

    ewald_ = new EnzoEwald (interp_xpoints_, interp_ypoints_, interp_zpoints_,
                            ewald_order_);

  } 
  
  compute_ (block);
}

//----------------------------------------------------------------------

double EnzoMethodMultipole::timestep ( Block * block ) throw()
{
  return timestep_(block);
}

//----------------------------------------------------------------------

double EnzoMethodMultipole::timestep_ (Block * block) throw()
{
  Field field = block->data()->field();

  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

#ifdef NEW_TIMESTEP  
  enzo_float * a3[3] =
    { (enzo_float*) field.values ("acceleration_x"),
      (enzo_float*) field.values ("acceleration_y"),
      (enzo_float*) field.values ("acceleration_z") };
#else
  enzo_float * ax = (enzo_float*) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float*) field.values ("acceleration_y");
  enzo_float * az = (enzo_float*) field.values ("acceleration_z");
#endif  

  const int rank = cello::rank();
  
  enzo_float dt = std::numeric_limits<enzo_float>::max();

#ifdef NEW_TIMESTEP  
  double h3[3];
  block->cell_width(h3,h3+1,h3+2);
#else  
  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);
#endif  
  
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (cosmology) {
    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    double dt = block->dt();
    double time = block->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
#ifdef NEW_TIMESTEP  
    if (rank >= 1) h3[0]*=cosmo_a;
    if (rank >= 2) h3[1]*=cosmo_a;
    if (rank >= 3) h3[2]*=cosmo_a;
#else
    if (rank >= 1) hx*=cosmo_a;
    if (rank >= 2) hy*=cosmo_a;
    if (rank >= 3) hz*=cosmo_a;
#endif 
  }
  
  double mean_cell_width;

#ifdef NEW_TIMESTEP
  if (rank == 1) mean_cell_width = h3[0];
  if (rank == 2) mean_cell_width = sqrt(h3[0]*h3[1]);
  if (rank == 3) mean_cell_width = cbrt(h3[0]*h3[1]*h3[2]);
#else
  if (rank == 1) mean_cell_width = hx;
  if (rank == 2) mean_cell_width = sqrt(hx*hy);
  if (rank == 3) mean_cell_width = cbrt(hx*hy*hz);
#endif
  
  // Timestep is sqrt(mean_cell_width / (a_mag_max + epsilon)),
  // where a_mag_max is the maximum acceleration magnitude
  // across all cells in the block, and epsilon defined as
  // mean_cell_width / dt_max_^2. This means that when acceleration
  // is zero everywhere, the timestep is equal to dt_max_
  
  const double epsilon = mean_cell_width / (dt_max_ * dt_max_);

  // Find the maximum of the square of the magnitude of acceleration
  // across all active cells, then get the square root of this value
  
  double a_mag_2_max = 0.0;
  double a_mag_2;

  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
	int i=ix + mx*(iy + iz*my);
#ifdef NEW_TIMESTEP
	if (rank == 1) a_mag_2 = a3[0][i] * a3[0][i];
	if (rank == 2) a_mag_2 = a3[0][i] * a3[0][i] + a3[1][i] * a3[1][i];
	if (rank == 3) a_mag_2 = a3[0][i] * a3[0][i] + a3[1][i] * a3[1][i]
			       + a3[2][i] * a3[2][i];
#else
	if (rank == 1) a_mag_2 = ax[i] * ax[i];
	if (rank == 2) a_mag_2 = ax[i] * ax[i] + ay[i] * ay[i];
	if (rank == 3) a_mag_2 = ax[i] * ax[i] + ay[i] * ay[i]
			       + az[i] * az[i];
#endif
	a_mag_2_max = std::max(a_mag_2_max,a_mag_2);
      }
    }
  }

  const double a_mag_max = sqrt(a_mag_2_max);
  dt = sqrt(mean_cell_width / (a_mag_max + epsilon)) ;

  return 0.5*dt;
}


//======================================================================


void EnzoMethodMultipole::compute_ (Block * block) throw()
{
  EnzoBlock * enzo_block = enzo::block(block);
  int level = enzo_block->level();
  
  if (block->is_leaf()) {
    
    compute_multipoles_ (block);
    begin_up_cycle_ (enzo_block);
    
  } 
  
  else {
    if (0 <= level && level < cello::hierarchy()->max_level() ) {
      restrict_recv(enzo_block, nullptr);
    }
  }
    
}


void EnzoMethodMultipole::begin_up_cycle_(EnzoBlock * enzo_block) throw()
{
  const int level = enzo_block->level();
  // CkPrintf("level = %d\n", level);
  // CkPrintf("min_level = %d\n", cello::hierarchy()->min_level());

  /*   for testing purposes     */
  // double * mass = pmass(enzo_block);
  // double * com = pcom(enzo_block);
  // double * quadrupole = pquadrupole(enzo_block);

  // if (level == 0) {
    
  //   CkPrintf("up cycle level: 0\n");
  //   CkPrintf("total mass: %f\n", *mass); 
  //   CkPrintf("COM: (%f, %f, %f)\n", com[0], com[1], com[2]);
  //   CkPrintf("quadrupole: ");
  //   for (int i = 0; i < 6; i++) {
  //     CkPrintf("%f  ", quadrupole[i]);
  //   }
  //   CkPrintf("\n");

  // }


  // should actually be if (level > min_level_)
  if (level > 0) {
    restrict_send (enzo_block);
  }


  // The volume is initialized before the barrier: after it, other Blocks may
  // already be sending update_volume() to this one, before this Block has
  // itself received the barrier's broadcast.
  *volume_(enzo_block) = block_volume_[level];

  CkCallback callback(CkIndex_EnzoBlock::r_method_multipole_dualwalk_barrier(nullptr),
		      enzo::block_array());
  enzo::block(enzo_block)->contribute(callback);

}


void EnzoBlock::r_method_multipole_dualwalk_barrier(CkReductionMsg* msg)
{
  EnzoMethodMultipole * method =
    static_cast<EnzoMethodMultipole*> (this->method());

  method->dual_walk_(this);
  
}

void EnzoMethodMultipole::dual_walk_(EnzoBlock * enzo_block) throw()
{

  const Index& index = enzo_block->index();
  int level = enzo_block->level();

  // if level == min_level_
  if (level == 0) {

    // call interact on all pairs of root blocks

    int nx, ny, nz;
    cello::hierarchy()->root_blocks(&nx, &ny, &nz);
    // const int rank = cello::rank();
    // if (rank < 3) nz = 1;
    // if (rank < 2) ny = 1;

    if (nx*ny*nz != 1) {

      // can I do iay = iaz, iaz < ny, etc?
      for (int iaz = 0; iaz < nz; iaz++) {
        for (int iay = 0; iay < ny; iay++) {
          for (int iax = 0; iax < nx; iax++) {

            const Index other_ind(iax, iay, iaz);

            if (index < other_ind) {
              continue;
            }
            else {
              // start the dual tree walk on all pairs of root blocks, including self-interactions
              // whether the other root block is a leaf is not known here
              // (-1): traverse() asks it
              enzo::block_array()[index].p_method_multipole_traverse(other_ind, -1);
            }
          }
        }
      }
    }

    else
    { 
      // start the dual tree walk with two copies of the root (if the root isn't a leaf)

      if (enzo_block->is_leaf()) {
	      begin_down_cycle_(enzo_block);
      }
      else {
	      enzo::block_array()[index].p_method_multipole_traverse(index, 0);
      }
    }

  }

  if (! enzo_block->is_leaf()) {

    // traverse terminates when all leaf blocks done. Since we use a
    // barrier on all blocks, we call the barrier now on all non-leaf
    // blocks
    CkCallback callback(CkIndex_EnzoBlock::r_method_multipole_traverse_complete(nullptr),
			enzo::block_array());
    enzo::block(enzo_block)->contribute(callback);
  }
  else {

    int nx, ny, nz;
    cello::hierarchy()->root_blocks(&nx, &ny, &nz);

    // A leaf is complete only after this callback has run (other Blocks
    // may have credited it all of its volume before): credit the extra
    // unit of volume that update_volume() waits for. (A lone root Block
    // that is a leaf does not wait for update_volume().)
    if (level > 0 || nx*ny*nz != 1) {
      update_volume (enzo_block, index, 1);
    }
  }
}

void EnzoBlock::r_method_multipole_traverse_complete(CkReductionMsg * msg)
{
  // delete msg; is this necessary?

  EnzoMethodMultipole * method =
    static_cast<EnzoMethodMultipole*> (this->method());

  // min_level
  if (level() == 0){ //only the block at the top invokes begin_down_cycle
      method->begin_down_cycle_(this);
  }
}

void EnzoMethodMultipole::begin_down_cycle_(EnzoBlock * enzo_block) throw() 
{

  // Field field = enzo_block->data()->field();
  // enzo_float * accel = (enzo_float *) field.values("acceleration_x"); 
  // Particle particle = enzo_block->data()->particle();
  // const int num_prtl_types = particle.num_types();


  // int mx, my, mz;
  // int gx, gy, gz;
  // double hx, hy, hz;
  // field.dimensions  (0, &mx, &my, &mz);
  // field.ghost_depth (0, &gx, &gy, &gz);
  
  if (enzo_block->is_leaf()) {

    evaluate_force_(enzo_block);
    // CkPrintf("accel_x at 10, 10, 10, pre-refresh: %.9f\n", accel[10 + 32 * (10 + 10 * 32)]);

    refresh_acceleration_(enzo_block);
  }

  else {
    prolong_send (enzo_block);
    // non-leaf Blocks take no further part in this Method
    enzo_block->compute_done();
  }
  
}

void EnzoMethodMultipole::refresh_acceleration_ (EnzoBlock * enzo_block) throw()
{
  cello::refresh(ir_exit_)->set_active(enzo_block->is_leaf());
  enzo_block->refresh_start
    (ir_exit_, CkIndex_EnzoBlock::p_method_multipole_end());
}

void EnzoBlock::p_method_multipole_end()
{

  // Particle particle = this->data()->particle();
  // const int num_prtl_types = particle.num_types();

  // for (int it = 0; it < num_prtl_types; it++) {
    
  //   for (int ib = 0; ib < particle.num_batches(it); ib++) {

  //     const int np = particle.num_particles(it,ib);

  //     const int ia_ax  = particle.attribute_index(it,"ax");
  //     const int ia_ay  = particle.attribute_index(it,"ay");
  //     const int ia_az  = particle.attribute_index(it,"az");

  //     enzo_float * axa =  (enzo_float *)particle.attribute_array (it,ia_ax,ib);
  //     enzo_float * aya =  (enzo_float *)particle.attribute_array (it,ia_ay,ib);
  //     enzo_float * aza =  (enzo_float *)particle.attribute_array (it,ia_az,ib);

  //     const int dax =  particle.stride(it,ia_ax);
  //     const int day =  particle.stride(it,ia_ay);
  //     const int daz =  particle.stride(it,ia_az);

  //     for (int ip=0; ip < np; ip++) {
  //       CkPrintf("particle accel: %.9f, %.9f, %.9f\n", axa[ip*dax], aya[ip*day], aza[ip*daz]);
  //     }
  //   }
  // }
  // Field field = this->data()->field();
  // enzo_float * accelx = (enzo_float *) field.values("acceleration_x"); 
  // enzo_float * accely = (enzo_float *) field.values("acceleration_y");
  // // enzo_float * vel = (enzo_float *) field.values("velocity_x"); 
  // CkPrintf("accel_x at 1, 1, 1, post-refresh: %.12f\n", accelx[1]); 
  // CkPrintf("accel_y at 1, 1, 1, post-refresh: %.12f\n", accely[1]);
  // CkPrintf("accel_x at 4, 4, 4, post-refresh: %.12f\n", accelx[2]); 
  // CkPrintf("accel_y at 4, 4, 4, post-refresh: %.12f\n", accely[2]);

  // CkPrintf("vel_x, post-refresh: %f\n", vel[0]);

  compute_done();
}



void EnzoMethodMultipole::restrict_send(EnzoBlock * enzo_block) throw()
{
  MultipoleMsg * msg = pack_multipole_(enzo_block);

  Index index_parent = enzo_block->index().index_parent(0);

  enzo::block_array()[index_parent].p_method_multipole_restrict_recv(msg);
}


void EnzoBlock::p_method_multipole_restrict_recv(MultipoleMsg * msg)
{
  EnzoMethodMultipole * method = 
    static_cast<EnzoMethodMultipole*> (this->method());

  method->restrict_recv(this,msg);
}


void EnzoMethodMultipole::restrict_recv
(EnzoBlock * enzo_block, MultipoleMsg * msg) throw()
{
  
  // Save multipole message from child
  if (msg != nullptr)
    {
      *pmsg_restrict(enzo_block, msg->child_index()) = msg;
    }
  
  // Continue if all expected messages received
  if ( psync_restrict(enzo_block)->next() ) {

    // Restore saved messages then clear
    for (int i = 0; i < cello::num_children(); i++) {
      msg = *pmsg_restrict(enzo_block, i);
      *pmsg_restrict(enzo_block, i) = nullptr;

      // Unpack multipoles from message then delete message
      unpack_multipole_(enzo_block, msg);
    }
  
    begin_up_cycle_ (enzo_block);
  }
  
}

void EnzoMethodMultipole::prolong_send(EnzoBlock * enzo_block) throw()
{

  ItChild it_child(cello::rank());
  int ic3[3]; 

  while (it_child.next(ic3)) {

    MultipoleMsg * msg = pack_coeffs_(enzo_block);

    Index index_child = enzo_block->index().index_child(ic3, cello::hierarchy()->min_level());

    // std::string string = enzo_block->name(index_child);
    // CkPrintf("sending to block name = %s\n", string.c_str());
    // fflush(stdout);

    enzo::block_array()[index_child].p_method_multipole_prolong_recv(msg);  

  }
}


void EnzoBlock::p_method_multipole_prolong_recv(MultipoleMsg * msg)
{
  EnzoMethodMultipole * method = 
    static_cast<EnzoMethodMultipole*> (this->method());

  method->prolong_recv(this, msg);
}


void EnzoMethodMultipole::prolong_recv
(EnzoBlock * enzo_block, MultipoleMsg * msg) throw()
{

  // Save multipole message from parent
  if (msg != nullptr)
    {

      MultipoleMsg ** temp = pmsg_prolong(enzo_block);

      if (temp == nullptr) {
        ERROR("early exit", "early exit");
      }
      *temp = msg;
    }

  
  // Return if not ready yet 
  //if (! psync_prolong(enzo_block)->next() ) return;

  // Restore saved message then clear
  msg = *pmsg_prolong(enzo_block);
  *pmsg_prolong(enzo_block) = nullptr;

  // Unpack coefficients from message then delete message
  unpack_coeffs_(enzo_block, msg);
    
  
  begin_down_cycle_ (enzo_block);
  
}


MultipoleMsg * EnzoMethodMultipole::pack_multipole_(EnzoBlock * enzo_block) throw()
{
 
  // Create a MultipoleMsg for sending data to parent
 
  MultipoleMsg * msg  = new MultipoleMsg;

  double * mass = pmass(enzo_block);
  double * com = pcom(enzo_block);
  double * quadrupole = pquadrupole(enzo_block);
 
  msg->mass = *mass;

  msg->com[0] = com[0];
  msg->com[1] = com[1];
  msg->com[2] = com[2];

  msg->quadrupole[0] = quadrupole[0];
  msg->quadrupole[1] = quadrupole[1];
  msg->quadrupole[2] = quadrupole[2];
  msg->quadrupole[3] = quadrupole[3];
  msg->quadrupole[4] = quadrupole[4];
  msg->quadrupole[5] = quadrupole[5];

  const int level = enzo_block->level();
  if (level > 0) {
    enzo_block->index().child(level,msg->ic3,msg->ic3+1,msg->ic3+2);
  }

  return msg;

}

MultipoleMsg * EnzoMethodMultipole::pack_coeffs_(EnzoBlock * enzo_block) throw()
{
 
  // Create a MultipoleMsg for sending data to children
 
  MultipoleMsg * msg  = new MultipoleMsg;

  double * com = pcom(enzo_block);
  double * c1 = pc1(enzo_block);
  double * c2 = pc2(enzo_block);
  double * c3 = pc3(enzo_block);

  msg->com[0] = com[0];
  msg->com[1] = com[1];
  msg->com[2] = com[2];

  msg->c1[0] = c1[0];
  msg->c1[1] = c1[1];
  msg->c1[2] = c1[2];

  msg->c2[0] = c2[0];
  msg->c2[1] = c2[1];
  msg->c2[2] = c2[2];
  msg->c2[3] = c2[3];
  msg->c2[4] = c2[4];
  msg->c2[5] = c2[5];

  msg->c3[0] = c3[0];
  msg->c3[1] = c3[1];
  msg->c3[2] = c3[2];
  msg->c3[3] = c3[3];
  msg->c3[4] = c3[4];
  msg->c3[5] = c3[5];
  msg->c3[6] = c3[6];
  msg->c3[7] = c3[7];
  msg->c3[8] = c3[8];
  msg->c3[9] = c3[9];

  return msg;

}


void EnzoMethodMultipole::unpack_multipole_
(EnzoBlock * enzo_block, MultipoleMsg * msg) throw()
{

  double * this_mass = pmass(enzo_block);
  double * this_com = pcom(enzo_block);
  double * this_quadrupole = pquadrupole(enzo_block);

  // copy data from msg to this EnzoBlock
  double child_mass = msg->mass;
  double * child_com = msg->com;
  double * child_quadrupole = msg->quadrupole;

  if (child_mass != 0) {

    double new_com[3];
    double tot_mass = *this_mass + child_mass;
    new_com[0] = (*this_mass * this_com[0] + child_mass * child_com[0]) / tot_mass;
    new_com[1] = (*this_mass * this_com[1] + child_mass * child_com[1]) / tot_mass;
    new_com[2] = (*this_mass * this_com[2] + child_mass * child_com[2]) / tot_mass;
    
    std::array<double, 6> shifted_this_quadrupole = shift_quadrupole_(this_quadrupole, *this_mass, this_com, new_com);
    std::array<double, 6> shifted_child_quadrupole = shift_quadrupole_(child_quadrupole, child_mass, child_com, new_com);

    this_quadrupole[0] = shifted_this_quadrupole[0] + shifted_child_quadrupole[0];
    this_quadrupole[1] = shifted_this_quadrupole[1] + shifted_child_quadrupole[1];
    this_quadrupole[2] = shifted_this_quadrupole[2] + shifted_child_quadrupole[2];
    this_quadrupole[3] = shifted_this_quadrupole[3] + shifted_child_quadrupole[3];
    this_quadrupole[4] = shifted_this_quadrupole[4] + shifted_child_quadrupole[4];
    this_quadrupole[5] = shifted_this_quadrupole[5] + shifted_child_quadrupole[5];

    this_com[0] = new_com[0];
    this_com[1] = new_com[1];
    this_com[2] = new_com[2];

    *this_mass += child_mass;

  }

  delete msg;
}

void EnzoMethodMultipole::unpack_coeffs_
(EnzoBlock * enzo_block, MultipoleMsg * msg) throw()
{
  // coefficient data from the current block
  double * pthis_com = pcom(enzo_block);
  double * pthis_c1 = pc1(enzo_block);
  double * pthis_c2 = pc2(enzo_block);
  double * pthis_c3 = pc3(enzo_block);

  std::array<double, 3> this_com;
  std::array<double, 3> this_c1;
  std::array<double, 6> this_c2;
  std::array<double, 10> this_c3;

  std::copy(pthis_com, pthis_com+3, this_com.begin());
  std::copy(pthis_c1, pthis_c1+3, this_c1.begin());
  std::copy(pthis_c2, pthis_c2+6, this_c2.begin());
  std::copy(pthis_c3, pthis_c3+10, this_c3.begin());
  
  // copy data from msg to this EnzoBlock
  double * pparent_com = msg->com;
  double * pparent_c1 = msg->c1;
  double * pparent_c2 = msg->c2;
  double * pparent_c3 = msg->c3;

  std::array<double, 3> parent_com;
  std::array<double, 3> parent_c1;
  std::array<double, 6> parent_c2;
  std::array<double, 10> parent_c3;

  std::copy(pparent_com, pparent_com+3, parent_com.begin());
  std::copy(pparent_c1, pparent_c1+3, parent_c1.begin());
  std::copy(pparent_c2, pparent_c2+6, parent_c2.begin());
  std::copy(pparent_c3, pparent_c3+10, parent_c3.begin());
  
  std::array<double, 3> com_shift = subtract_(parent_com, this_com);

  std::array<double, 3> shifted_parent_c1_secondterm = dot_12_(com_shift, parent_c2);
  std::array<double, 3> shifted_parent_c1_thirdterm = dot_23_(outer_11_(com_shift, com_shift), dot_scalar_3_(0.5, parent_c3));
  std::array<double, 3> shifted_parent_c1 = add_11_(add_11_(parent_c1, shifted_parent_c1_secondterm), shifted_parent_c1_thirdterm);
  
  std::array<double, 6> shifted_parent_c2_secondterm = dot_13_(com_shift, parent_c3);
  std::array<double, 6> shifted_parent_c2 = add_22_(parent_c2, shifted_parent_c2_secondterm);
   
  std::array<double, 3> new_c1 = add_11_(this_c1, shifted_parent_c1);   
  std::array<double, 6> new_c2 = add_22_(this_c2, shifted_parent_c2);
  std::array<double, 10> new_c3 = add_33_(this_c3, parent_c3);

  pc1(enzo_block)[0] = new_c1[0];
  pc1(enzo_block)[1] = new_c1[1];
  pc1(enzo_block)[2] = new_c1[2];
  
  pc2(enzo_block)[0] = new_c2[0];
  pc2(enzo_block)[1] = new_c2[1];
  pc2(enzo_block)[2] = new_c2[2];
  pc2(enzo_block)[3] = new_c2[3];
  pc2(enzo_block)[4] = new_c2[4];
  pc2(enzo_block)[5] = new_c2[5];

  pc3(enzo_block)[0] = new_c3[0];
  pc3(enzo_block)[1] = new_c3[1];
  pc3(enzo_block)[2] = new_c3[2];
  pc3(enzo_block)[3] = new_c3[3];
  pc3(enzo_block)[4] = new_c3[4];
  pc3(enzo_block)[5] = new_c3[5];
  pc3(enzo_block)[6] = new_c3[6];
  pc3(enzo_block)[7] = new_c3[7];
  pc3(enzo_block)[8] = new_c3[8];
  pc3(enzo_block)[9] = new_c3[9];
  

  delete msg;
}


/************************************************************************/
/// FMM functions


void EnzoMethodMultipole::compute_multipoles_ (Block * block) throw()
{  
  Field field = block->data()->field();

  enzo_float * density = (enzo_float *) field.values("density"); 
  enzo_float * accel = (enzo_float *) field.values("acceleration_x"); 
  enzo_float * vel = (enzo_float *) field.values("velocity_x"); 

  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;

  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);

  block->cell_width(&hx, &hy, &hz);
  double cell_vol = hx * hy * hz;

  double lo[3];
  block->lower(lo, lo+1, lo+2);

  Particle particle = block->data()->particle();
  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();
  const int num_is_grav = particle_groups->size("is_gravitating");

  // distinguish between mass as attribute vs. mass as constant
  enzo_float * prtmass = NULL;
  int dm;

  std::array<double, 3> com_sum = {};
  std::array<double, 6> quad_sum = {};

  double * mass = pmass(block);
  double * com = pcom(block);
  double * quadrupole = pquadrupole(block);

  // CkPrintf("mass 1: %f\n", *mass);
  // CkPrintf("dens: %f, %f, %f\n", density[10], density[11], density[12]);

  // loop over cells
  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {
        
        double dens = density[ix + mx * (iy + iz * my)];

        // if (iz == gz) {
        //   CkPrintf("(%d, %d): dens = %f\n", ix-gx, iy-gy, dens);
        // }

        // if (std::isnan(dens)) {
        //   CkPrintf("bad density: %d, %d, %d\n", ix, iy, iz);
        //   CkPrintf("bad acceleration: %f\n", accel[ix + mx * (iy + iz * my)]);
        //   CkPrintf("bad velocity: %f\n", vel[ix + mx * (iy + iz * my)]);
        // }

        double cell_mass = dens * cell_vol;
        double pos[3] = {lo[0] + (ix-gx + 0.5)*hx, lo[1] + (iy-gy + 0.5)*hy, lo[2] + (iz-gz + 0.5)*hz};

        *mass += cell_mass;

        com_sum[0] += cell_mass * pos[0];
        com_sum[1] += cell_mass * pos[1];
        com_sum[2] += cell_mass * pos[2];

        quad_sum[0] += cell_mass * pos[0] * pos[0];
        quad_sum[1] += cell_mass * pos[0] * pos[1];
        quad_sum[2] += cell_mass * pos[0] * pos[2];
        quad_sum[3] += cell_mass * pos[1] * pos[1];
        quad_sum[4] += cell_mass * pos[1] * pos[2];
        quad_sum[5] += cell_mass * pos[2] * pos[2];

      }
    }
  }

  // CkPrintf("mass 2: %f\n", *mass);

  // loop over particles
  for (int ipt = 0; ipt < num_is_grav; ipt++) {
    const int it = particle.type_index(particle_groups->item("is_gravitating",ipt));

    int imass = 0;

    // check correct precision for position
    int ia = particle.attribute_index(it,"x");
    int ba = particle.attribute_bytes(it,ia); // "bytes (actual)"
    int be = sizeof(enzo_float);                // "bytes (expected)"

      ASSERT4 ("EnzoMethodPmUpdate::compute()",
	       "Particle type %s attribute %s defined as %s but expecting %s",
	       particle.type_name(it).c_str(),
	       particle.attribute_name(it,ia).c_str(),
	       ((ba == 4) ? "single" : ((ba == 8) ? "double" : "quadruple")),
	       ((be == 4) ? "single" : ((be == 8) ? "double" : "quadruple")),
	       (ba == be));

    for (int ib = 0; ib < particle.num_batches(it); ib++) {

      const int np = particle.num_particles(it,ib);

      if (particle.has_attribute(it,"mass")) {

        // Particle type has an attribute called "mass".
        // In this case we set prtmass to point to the mass attribute array
        // Also set dm to be the stride for the "mass" attribute
        imass = particle.attribute_index(it,"mass");
        prtmass = (enzo_float *) particle.attribute_array( it, imass, ib);
        dm = particle.stride(it,imass);

	    } 
      else {

        // Particle type has a constant called "mass".
        // In this case we set prtmass to point to the value
        // of the mass constant.
        // dm is set to 0, which will mean that we can loop through an
        // "array" of length 1.
        imass = particle.constant_index(it,"mass");
        prtmass = (enzo_float*)particle.constant_value(it,imass);
        dm = 0;
      }

      const int ia_x  = particle.attribute_index(it,"x");
      const int ia_y  = particle.attribute_index(it,"y");
      const int ia_z  = particle.attribute_index(it,"z");

      enzo_float * xa =  (enzo_float *)particle.attribute_array (it,ia_x,ib);
      enzo_float * ya =  (enzo_float *)particle.attribute_array (it,ia_y,ib);
      enzo_float * za =  (enzo_float *)particle.attribute_array (it,ia_z,ib);

      const int dx =  particle.stride(it,ia_x);
      const int dy =  particle.stride(it,ia_y);
      const int dz =  particle.stride(it,ia_z);

      for (int ip=0; ip < np; ip++) {
        
        *mass += prtmass[ip*dm];

        double pos[3] = {xa[ip*dx], ya[ip*dy], za[ip*dz]};
        
        com_sum[0] += prtmass[ip*dm] * pos[0];
        com_sum[1] += prtmass[ip*dm] * pos[1];
        com_sum[2] += prtmass[ip*dm] * pos[2];
        
        quad_sum[0] += prtmass[ip*dm] * pos[0] * pos[0];
        quad_sum[1] += prtmass[ip*dm] * pos[0] * pos[1];
        quad_sum[2] += prtmass[ip*dm] * pos[0] * pos[2];
        quad_sum[3] += prtmass[ip*dm] * pos[1] * pos[1];
        quad_sum[4] += prtmass[ip*dm] * pos[1] * pos[2];
        quad_sum[5] += prtmass[ip*dm] * pos[2] * pos[2];

      }
    }
  }

  if (*mass != 0) {

    com[0] = com_sum[0] / *mass;
    com[1] = com_sum[1] / *mass;
    com[2] = com_sum[2] / *mass;

    quadrupole[0] = quad_sum[0] - *mass * com[0] * com[0];
    quadrupole[1] = quad_sum[1] - *mass * com[0] * com[1];
    quadrupole[2] = quad_sum[2] - *mass * com[0] * com[2];
    quadrupole[3] = quad_sum[3] - *mass * com[1] * com[1];
    quadrupole[4] = quad_sum[4] - *mass * com[1] * com[2];
    quadrupole[5] = quad_sum[5] - *mass * com[2] * com[2];

  }
}
 
void EnzoMethodMultipole::evaluate_force_(Block * block) throw()
{
  Data * data = block->data();
  Field field = data->field();
  
  enzo_float * accel_x = (enzo_float*) field.values ("acceleration_x");
  enzo_float * accel_y = (enzo_float*) field.values ("acceleration_y");
  enzo_float * accel_z = (enzo_float*) field.values ("acceleration_z");

  // CkPrintf("accel_x before evaluate_force: %f\n", accel_x[0]);

  enzo_float * dens = (enzo_float*) field.values ("density");

  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;

  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);

  // hz = 1 by default
  block->cell_width(&hx, &hy, &hz);
  double cell_vol = hx * hy * hz;
  // CkPrintf("hx, hy, hz: %f, %f, %f\n", hx, hy, hz);

  double lo[3];
  block->lower(lo, lo+1, lo+2);

  Particle particle = block->data()->particle();
  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();
  const int num_is_grav = particle_groups->size("is_gravitating");
  const int num_prtl_types = particle_descr->num_types();

  // distinguish between mass as attribute vs. mass as constant
  enzo_float * prtmass2 = NULL;
  int dm2;

  double mass = *pmass(block);
  double * pthis_com = pcom(block);
  double * pthis_quadrupole = pquadrupole(block);
  double * pthis_c1 = pc1(block);
  double * pthis_c2 = pc2(block);
  double * pthis_c3 = pc3(block);

  std::array<double, 3> com;
  std::array<double, 6> quadrupole;
  std::array<double, 3> c1;
  std::array<double, 6> c2;
  std::array<double, 10> c3;

  std::copy(pthis_com, pthis_com+3, com.begin());
  std::copy(pthis_quadrupole, pthis_quadrupole+6, quadrupole.begin());
  std::copy(pthis_c1, pthis_c1+3, c1.begin());
  std::copy(pthis_c2, pthis_c2+6, c2.begin());
  std::copy(pthis_c3, pthis_c3+10, c3.begin());

  // compute long-range contribution from periodic images of this Block 
  /* 
  if (ewald_ != nullptr) {

    // CkPrintf("ewald in evaluate_force\n");

    CelloView<double, 1> d2_gadget = ewald_->d2_gadget(0, 0, 0);

    std::array<double, 3> d1_ewald = ewald_->interp_d1(0, 0, 0); // d1 contribution from periodicity
    std::array<double, 6> d2_ewald = ewald_->interp_d2(0, 0, 0); // d2 contribution from periodicity
    std::array<double, 10> d3_ewald = ewald_->interp_d3(0, 0, 0); // d3 contribution from periodicity

    CkPrintf("d1_ewald in evaluate_force: %f, %f, %f\n", d1_ewald[0], d1_ewald[1], d1_ewald[2]);
    CkPrintf("d2_ewald in evaluate_force: %f, %f, %f, %f, %f, %f\n", 
          d2_ewald[0], d2_ewald[1], d2_ewald[2], d2_ewald[3], d2_ewald[4], d2_ewald[5]);
    CkPrintf("d3_ewald in evaluate_force: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
          d3_ewald[0], d3_ewald[1], d3_ewald[2], d3_ewald[3], d3_ewald[4], d3_ewald[5],
          d3_ewald[6], d3_ewald[7], d3_ewald[8], d3_ewald[9]);

    // compute the coefficients of the Taylor expansion of acceleration due to the particles in Block b
    std::array<double, 3> delta_c1 = add_11_(dot_scalar_1_(mass, d1_ewald), dot_23_(quadrupole, dot_scalar_3_(0.5, d3_ewald)));
    std::array<double, 6> delta_c2 = dot_scalar_2_(mass, d2_ewald);
    std::array<double, 10> delta_c3 = dot_scalar_3_(mass, d3_ewald);

    // CkPrintf("delta_c1 in evaluate_force: %f, %f, %f\n", delta_c1[0], delta_c1[1], delta_c1[2]);
    // CkPrintf("delta_c2 in evaluate_force: %f, %f, %f, %f, %f, %f\n", delta_c2[0], delta_c2[1], delta_c2[2], delta_c2[3], delta_c2[4], delta_c2[5]);
    // CkPrintf("delta_c3 in evaluate_force: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", delta_c3[0], delta_c3[1], delta_c3[2], delta_c3[3], delta_c3[4], delta_c3[5], delta_c3[6], delta_c3[7], delta_c3[8], delta_c3[9]);

    
    // add the coefficients for the new interaction to the coefficients already associated with this Block
    std::array<double, 3> new_c1 = add_11_(c1, delta_c1);  
    std::array<double, 6> new_c2 = add_22_(c2, delta_c2);
    std::array<double, 10> new_c3 = add_33_(c3, delta_c3);

    c1[0] = new_c1[0];
    c1[1] = new_c1[1];
    c1[2] = new_c1[2];
    
    c2[0] = new_c2[0];
    c2[1] = new_c2[1];
    c2[2] = new_c2[2];
    c2[3] = new_c2[3];
    c2[4] = new_c2[4];
    c2[5] = new_c2[5];

    c3[0] = new_c3[0];
    c3[1] = new_c3[1];
    c3[2] = new_c3[2];
    c3[3] = new_c3[3];
    c3[4] = new_c3[4];
    c3[5] = new_c3[5];
    c3[6] = new_c3[6];
    c3[7] = new_c3[7];
    c3[8] = new_c3[8];
    c3[9] = new_c3[9];
  }
  */

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  double grav_const; 
  if (cosmology) {
    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    double dt = block->dt();
    double time = block->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    grav_const = 1./(4 * cello::pi * cosmo_a * cosmo_a);
  }
  else {
    grav_const = enzo::grav_constant_codeU();
  }

  // CkPrintf("grav_constant: %f\n", grav_const);
  
  // loop over all cells
  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {

        int i = ix + mx*(iy + my*iz);

        // compute force sourced from other Blocks
        std::array<double, 3> a;
        a[0] = (lo[0] + (ix-gx + 0.5)*hx) - com[0]; 
        a[1] = (lo[1] + (iy-gy + 0.5)*hy) - com[1];
        a[2] = (lo[2] + (iz-gz + 0.5)*hz) - com[2];
        
        std::array<double, 3> second_term = dot_12_(a, c2);
        std::array<double, 3> third_term = dot_23_(outer_11_(a, a), dot_scalar_3_(0.5, c3));

        std::array<double, 3> block_force = add_11_(subtract_(c1, second_term), third_term);

        // CkPrintf("first term: %f, %f, %f\n", c1[0], c1[1], c1[2]);
        // CkPrintf("c2: %f, %f, %f, %f, %f, %f\n", c2[0], c2[1], c2[2], c2[3], c2[4], c2[5]);
	      // CkPrintf("a: %f, %f, %f\n", a[0], a[1], a[2]);
        // CkPrintf("second term: %f, %f, %f\n", second_term[0], second_term[1], second_term[2]);
        // CkPrintf("third term: %f, %f, %f\n", third_term[0], third_term[1], third_term[2]);


        // subtracting rather than adding since block_force is multiplied by -G 
        accel_x[i] -= grav_const * block_force[0];
        accel_y[i] -= grav_const * block_force[1];
        accel_z[i] -= grav_const * block_force[2];

        // CkPrintf("position: %f, %f, %f\n", (lo[0] + (ix-gx + 0.5)*hx), (lo[1] + (iy-gy + 0.5)*hy), (lo[2] + (iz-gz + 0.5)*hz));
        // CkPrintf("Block force: %f, %f, %f\n", -1*grav_const*block_force[0], -1*grav_const*block_force[1], -1*grav_const*block_force[2]);

        std::array<double, 3> tot_cell_force = {}; // for debugging purposes

        // compute force sourced from other cells
        for (int iz2 = gz; iz2 < mz-gz; iz2++) {
          for (int iy2 = gy; iy2 < my-gy; iy2++) {
            for (int ix2 = gx; ix2 < mx-gx; ix2++) {

              int i2 = ix2 + mx*(iy2 + my*iz2);

              if (i != i2) {

                std::array<double, 3> disp;

                double b_mass = dens[i2]*cell_vol;
                // CkPrintf("b_mass: %f\n", b_mass);
                
                if (ewald_ != nullptr) {

                  double a_loc[3]{lo[0] + (ix-gx + 0.5)*hx, lo[1] + (iy-gy + 0.5)*hy, lo[2] + (iz-gz + 0.5)*hz};
                  double b_loc[3]{lo[0] + (ix2-gx + 0.5)*hx, lo[1] + (iy2-gy + 0.5)*hy, lo[2] + (iz2-gz + 0.5)*hz};
                  double b_image[3];
                  cello::hierarchy()->get_nearest_periodic_image(b_loc, a_loc, b_image);

                  // CkPrintf("a_loc: %f, %f, %f\n", a_loc[0], a_loc[1], a_loc[2]);
                  // CkPrintf("b_loc: %f, %f, %f\n", b_loc[0], b_loc[1], b_loc[2]);
                  // CkPrintf("b_image: %f, %f, %f\n", b_image[0], b_image[1], b_image[2]);
                
                  // displacement vector pointing from current cell to interacting cell
                  const int rank = cello::rank();
                  disp[0] = b_image[0] - a_loc[0];
                  disp[1] = (rank < 2) ? 0 : b_image[1] - a_loc[1];
                  disp[2] = (rank < 3) ? 0 : b_image[2] - a_loc[2];
                  
                  // acceleration contribution from periodic images of particle b
                  std::array<double, 3> d1_ewald = ewald_->interp_d1(disp[0], disp[1], disp[2]); 
                  
                  // should i be subtracting or adding?
                  // minus sign from negative gradient * minus sign from chain rule
                  accel_x[i] += grav_const * b_mass * d1_ewald[0];
                  accel_y[i] += grav_const * b_mass * d1_ewald[1];
                  accel_z[i] += grav_const * b_mass * d1_ewald[2];

                  tot_cell_force[0] += grav_const * b_mass * d1_ewald[0];
                  tot_cell_force[1] += grav_const * b_mass * d1_ewald[1];
                  tot_cell_force[2] += grav_const * b_mass * d1_ewald[2];

                  // if ((iz == gz) && (iz2==gz) && (ix-gx==3) && (iy-gy==20)) {
                  //    CkPrintf("(%d, %d, %f) periodic force: (%f, %f); accel: (%.9f, %.9f)\n",
                  //    ix2-gx, iy2-gy, dens[i2],
                  //    grav_const * b_mass * d1_ewald[0], 
                  //    grav_const * b_mass * d1_ewald[1], 
                  //    accel_x[i], accel_y[i]);
                  // }
                }
                else {
                  disp[0] = (ix2 - ix) * hx;
                  disp[1] = (iy2 - iy) * hy;
                  disp[2] = (iz2 - iz) * hz;
                }
                
                std::array<double, 3> cell_force = newton_force_(b_mass, disp); 

                // CkPrintf("nonperiodic force: %f, %f, %f\n", 
                //     grav_const*cell_force[0], grav_const*cell_force[1], grav_const*cell_force[2]);

                accel_x[i] += grav_const * cell_force[0];
                accel_y[i] += grav_const * cell_force[1];
                accel_z[i] += grav_const * cell_force[2];

                tot_cell_force[0] += grav_const * cell_force[0];
                tot_cell_force[1] += grav_const * cell_force[1];
                tot_cell_force[2] += grav_const * cell_force[2];

                // if ((iz == gz) && (iz2==gz) && (ix-gx==3) && (iy-gy==20)) {
                //      CkPrintf("(%d, %d, %f) nonperiodic force: (%f, %f); accel: (%.9f, %.9f)\n",
                //      ix2-gx, iy2-gy, dens[i2], 
                //      grav_const * cell_force[0], 
                //      grav_const * cell_force[1], 
                //      accel_x[i], accel_y[i]);
                //   }

              }
            }
          }
        }
	      // CkPrintf("Cell loc: %f, %f, %f\n", lo[0] + (ix-gx + 0.5)*hx, lo[1] + (iy-gy + 0.5)*hy, lo[2] + (iz-gz + 0.5)*hz);   
        // CkPrintf("Cell force: %f, %f, %f\n", tot_cell_force[0], tot_cell_force[1], tot_cell_force[2]);
        // CkPrintf("Cell accel: %f, %f, %f\n", accel_x[i], accel_y[i], accel_z[i]);

        std::array<double, 3> tot_prt_force = {}; // for debugging purposes

        // compute cell force sourced from particles
        for (int ipt = 0; ipt < num_is_grav; ipt++) {
          const int it = particle.type_index(particle_groups->item("is_gravitating",ipt));

          int imass2 = 0;

          for (int ib = 0; ib < particle.num_batches(it); ib++) {

            const int np = particle.num_particles(it,ib);

            if (particle.has_attribute(it,"mass")) {
              imass2 = particle.attribute_index(it,"mass");
              prtmass2 = (enzo_float *) particle.attribute_array( it, imass2, ib);
              dm2 = particle.stride(it,imass2);
            } 
            else {
              imass2 = particle.constant_index(it,"mass");
              prtmass2 = (enzo_float*)particle.constant_value(it,imass2);
              dm2 = 0;
            }

            const int ia_x2  = particle.attribute_index(it,"x");
            const int ia_y2  = particle.attribute_index(it,"y");
            const int ia_z2  = particle.attribute_index(it,"z");

            enzo_float * xa2 =  (enzo_float *)particle.attribute_array (it,ia_x2,ib);
            enzo_float * ya2 =  (enzo_float *)particle.attribute_array (it,ia_y2,ib);
            enzo_float * za2 =  (enzo_float *)particle.attribute_array (it,ia_z2,ib);

            const int dx2 =  particle.stride(it,ia_x2);
            const int dy2 =  particle.stride(it,ia_y2);
            const int dz2 =  particle.stride(it,ia_z2);

            for (int ip=0; ip < np; ip++) {
              
              // disp points from current cell to particle
              std::array<double, 3> disp;
              
              if (ewald_ != nullptr) {
                
                double a_loc[3]{lo[0] + (ix-gx + 0.5)*hx, lo[1] + (iy-gy + 0.5)*hy, lo[2] + (iz-gz + 0.5)*hz};
                double b_loc[3]{xa2[ip*dx2], ya2[ip*dy2], za2[ip*dz2]};
                double b_image[3];
                cello::hierarchy()->get_nearest_periodic_image(b_loc, a_loc, b_image);

                const int rank = cello::rank();
                disp[0] = b_image[0] - a_loc[0];
                disp[1] = (rank < 2) ? 0 : b_image[1] - a_loc[1];
                disp[2] = (rank < 3) ? 0 : b_image[2] - a_loc[2];
                  
                // acceleration contribution from periodic images of Particle b
                std::array<double, 3> d1_ewald = ewald_->interp_d1(disp[0], disp[1], disp[2]); 
                  
                accel_x[i] += grav_const * prtmass2[ip*dm2] * d1_ewald[0];
                accel_y[i] += grav_const * prtmass2[ip*dm2] * d1_ewald[1];
                accel_z[i] += grav_const * prtmass2[ip*dm2] * d1_ewald[2];

                // tot_prt_force[0] += grav_const * prtmass2[ip*dm2] * d1_ewald[0];
                // tot_prt_force[1] += grav_const * prtmass2[ip*dm2] * d1_ewald[1];
                // tot_prt_force[2] += grav_const * prtmass2[ip*dm2] * d1_ewald[2];
              }
              else {
                disp[0] = xa2[ip*dx2] - (lo[0] + (ix-gx + 0.5)*hx); 
                disp[1] = ya2[ip*dy2] - (lo[1] + (iy-gy + 0.5)*hy);
                disp[2] = za2[ip*dz2] - (lo[2] + (iz-gz + 0.5)*hz);
              }
                
              std::array<double, 3> prtcell_force = newton_force_(prtmass2[ip*dm2], disp); 

              accel_x[i] += grav_const * prtcell_force[0];
              accel_y[i] += grav_const * prtcell_force[1];
              accel_z[i] += grav_const * prtcell_force[2];

              // tot_prt_force[0] += grav_const * prtcell_force[0];
              // tot_prt_force[1] += grav_const * prtcell_force[1];
              // tot_prt_force[2] += grav_const * prtcell_force[2];

              // ALSO UPDATE GRAVITATING PARTICLE ACCELS HERE?
            }
          }
        }

        // CkPrintf("Particle force: %f, %f, %f\n", tot_prt_force[0], tot_prt_force[1], tot_prt_force[2]);
        // CkPrintf("Accel: %f, %f, %f\n\n", accel_x[0], accel_y[0], accel_z[0]);
      }
    }
  }


  // loop over all particles
  // change this to loop over all particles, not just the ones that are gravitating (done?)

  for (int it = 0; it < num_prtl_types; it++) {
    
    for (int ib = 0; ib < particle.num_batches(it); ib++) {

      const int np = particle.num_particles(it,ib);

      const int ia_x  = particle.attribute_index(it,"x");
      const int ia_y  = particle.attribute_index(it,"y");
      const int ia_z  = particle.attribute_index(it,"z");
      const int ia_ax  = particle.attribute_index(it,"ax");
      const int ia_ay  = particle.attribute_index(it,"ay");
      const int ia_az  = particle.attribute_index(it,"az");

      enzo_float * xa =  (enzo_float *)particle.attribute_array (it,ia_x,ib);
      enzo_float * ya =  (enzo_float *)particle.attribute_array (it,ia_y,ib);
      enzo_float * za =  (enzo_float *)particle.attribute_array (it,ia_z,ib);
      enzo_float * axa =  (enzo_float *)particle.attribute_array (it,ia_ax,ib);
      enzo_float * aya =  (enzo_float *)particle.attribute_array (it,ia_ay,ib);
      enzo_float * aza =  (enzo_float *)particle.attribute_array (it,ia_az,ib);

      const int dx =  particle.stride(it,ia_x);
      const int dy =  particle.stride(it,ia_y);
      const int dz =  particle.stride(it,ia_z);
      const int dax =  particle.stride(it,ia_ax);
      const int day =  particle.stride(it,ia_ay);
      const int daz =  particle.stride(it,ia_az);

      for (int ip=0; ip < np; ip++) {

        // compute force sourced by other blocks
        std::array<double, 3> a;
        a[0] =  xa[ip*dx] - com[0]; 
        a[1] =  ya[ip*dy] - com[1];
        a[2] =  za[ip*dz] - com[2];

        std::array<double, 3> second_term = dot_12_(a, c2);
        std::array<double, 3> third_term = dot_23_(outer_11_(a, a), dot_scalar_3_(0.5, c3));

        // how does the code treat G?
        std::array<double, 3> block_force = add_11_(subtract_(c1, second_term), third_term);

        // subtracting rather than adding since block_force is multiplied by -G
        axa[ip*dax] -= grav_const * block_force[0];
        aya[ip*day] -= grav_const * block_force[1];
        aza[ip*daz] -= grav_const * block_force[2];

        //CkPrintf("position (prt): %f, %f, %f\n", xa[ip*dx], ya[ip*dy], za[ip*dz]);
        //CkPrintf("Block force (prt): %f, %f, %f\n", -1.0*block_force[0], -1.0*block_force[1], -1.0*block_force[2]);

        std::array<double, 3> tot_cell_force = {}; // for debugging purposes
        
        // compute force sourced by cells
        for (int iz2 = gz; iz2 < mz-gz; iz2++) {
          for (int iy2 = gy; iy2 < my-gy; iy2++) {
            for (int ix2 = gx; ix2 < mx-gx; ix2++) {

              int i2 = ix2 + mx*(iy2 + my*iz2);

              // disp points from current particle to interacting cell
              std::array<double, 3> disp;

              double b_mass = dens[i2]*cell_vol;

              // force from periodic images of Cell b
              if (ewald_ != nullptr) {

                double a_loc[3]{xa[ip*dx], ya[ip*dy], za[ip*dz]};
                double b_loc[3]{lo[0] + (ix2-gx + 0.5)*hx, lo[1] + (iy2-gy + 0.5)*hy, lo[2] + (iz2-gz + 0.5)*hz};
                double b_image[3];
                cello::hierarchy()->get_nearest_periodic_image(b_loc, a_loc, b_image);

                const int rank = cello::rank();
                disp[0] = b_image[0] - a_loc[0];
                disp[1] = (rank < 2) ? 0 : b_image[1] - a_loc[1];
                disp[2] = (rank < 3) ? 0 : b_image[2] - a_loc[2];
                  
                // acceleration contribution from periodic images of Cell b
                std::array<double, 3> d1_ewald = ewald_->interp_d1(disp[0], disp[1], disp[2]); 
                  
                axa[ip*dax] += grav_const * b_mass * d1_ewald[0];
                aya[ip*day] += grav_const * b_mass * d1_ewald[1];
                aza[ip*daz] += grav_const * b_mass * d1_ewald[2];

                // CkPrintf("cell->particle force (periodic): %f, %f, %f\n",
                //     grav_const * b_mass * d1_ewald[0], grav_const * b_mass * d1_ewald[1], grav_const * b_mass * d1_ewald[2]);

                // tot_cell_force[0] += grav_const * b_mass * d1_ewald[0];
                // tot_cell_force[1] += grav_const * b_mass * d1_ewald[1];
                // tot_cell_force[2] += grav_const * b_mass * d1_ewald[2];
              }
              else {
                disp[0] = (lo[0] + (ix2-gx + 0.5)*hx) - xa[ip*dx]; 
                disp[1] = (lo[1] + (iy2-gy + 0.5)*hy) - ya[ip*dy];
                disp[2] = (lo[2] + (iz2-gz + 0.5)*hz) - za[ip*dz];
              }
                
              std::array<double, 3> cellprt_force = newton_force_(b_mass, disp); 

              axa[ip*dax] += grav_const * cellprt_force[0];
              aya[ip*day] += grav_const * cellprt_force[1];
              aza[ip*daz] += grav_const * cellprt_force[2];

              // CkPrintf("cell->particle force (nonperiodic): %f, %f, %f\n", 
              //     grav_const * cellprt_force[0], grav_const * cellprt_force[1], grav_const * cellprt_force[2]);

              // tot_cell_force[0] += grav_const * cellprt_force[0];
              // tot_cell_force[1] += grav_const * cellprt_force[1];
              // tot_cell_force[2] += grav_const * cellprt_force[2];

            }
          }
        }

        //CkPrintf("Cell force (prt): %f, %f, %f\n", tot_cell_force[0], tot_cell_force[1], tot_cell_force[2]);

        std::array<double, 3> tot_prt_force = {}; // for debugging purposes

        // compute force sourced by other particles
        for (int ipt2 = 0; ipt2 < num_is_grav; ipt2++) {
          const int it2 = particle.type_index(particle_groups->item("is_gravitating",ipt2));

          int imass2 = 0;

          for (int ib2 = 0; ib2 < particle.num_batches(it2); ib2++) {

            const int np2 = particle.num_particles(it2,ib2);

            if (particle.has_attribute(it2,"mass")) {
              imass2 = particle.attribute_index(it2,"mass");
              prtmass2 = (enzo_float *) particle.attribute_array( it2, imass2, ib2);
              dm2 = particle.stride(it2,imass2);
            } 
            else {
              imass2 = particle.constant_index(it2,"mass");
              prtmass2 = (enzo_float*)particle.constant_value(it2,imass2);
              dm2 = 0;
            }

            const int ia_x2  = particle.attribute_index(it2,"x");
            const int ia_y2  = particle.attribute_index(it2,"y");
            const int ia_z2  = particle.attribute_index(it2,"z");

            enzo_float * xa2 =  (enzo_float *)particle.attribute_array (it2,ia_x2,ib2);
            enzo_float * ya2 =  (enzo_float *)particle.attribute_array (it2,ia_y2,ib2);
            enzo_float * za2 =  (enzo_float *)particle.attribute_array (it2,ia_z2,ib2);

            const int dx2 =  particle.stride(it2,ia_x2);
            const int dy2 =  particle.stride(it2,ia_y2);
            const int dz2 =  particle.stride(it2,ia_z2);

            for (int ip2=0; ip2 < np2; ip2++) {

              if (!(it == it2 && ib == ib2 && ip == ip2)) {
              
                // disp points from current particle to interacting particle
                std::array<double, 3> disp;

                // force from periodic images of Particle b
                if (ewald_ != nullptr) {

                  double a_loc[3]{xa[ip*dx], ya[ip*dy], za[ip*dz]};
                  double b_loc[3]{xa2[ip2*dx2], ya2[ip2*dy2], za2[ip2*dz2]};
                  double b_image[3];
                  cello::hierarchy()->get_nearest_periodic_image(b_loc, a_loc, b_image);

                  const int rank = cello::rank();
                  disp[0] = b_image[0] - a_loc[0];
                  disp[1] = (rank < 2) ? 0 : b_image[1] - a_loc[1];
                  disp[2] = (rank < 3) ? 0 : b_image[2] - a_loc[2];
                    
                  // acceleration contribution from periodic images of Particle b
                  std::array<double, 3> d1_ewald = ewald_->interp_d1(disp[0], disp[1], disp[2]); 
                    
                  axa[ip*dax] += grav_const * prtmass2[ip2*dm2] * d1_ewald[0];
                  aya[ip*day] += grav_const * prtmass2[ip2*dm2] * d1_ewald[1];
                  aza[ip*daz] += grav_const * prtmass2[ip2*dm2] * d1_ewald[2];

                  // CkPrintf("particle->particle force (periodic): %f, %f, %f\n", 
                  //     grav_const * prtmass2[ip2*dm2] * d1_ewald[0], grav_const * prtmass2[ip2*dm2] * d1_ewald[1], grav_const * prtmass2[ip2*dm2] * d1_ewald[2]);

                  tot_prt_force[0] += grav_const * prtmass2[ip2*dm2] * d1_ewald[0];
                  tot_prt_force[1] += grav_const * prtmass2[ip2*dm2] * d1_ewald[1];
                  tot_prt_force[2] += grav_const * prtmass2[ip2*dm2] * d1_ewald[2];
                }
                else {
                  disp[0] = xa2[ip2*dx2] - xa[ip*dx]; 
                  disp[1] = ya2[ip2*dy2] - ya[ip*dy];
                  disp[2] = za2[ip2*dz2] - za[ip*dz];
                }
                  
                std::array<double, 3> prtprt_force = newton_force_(prtmass2[ip2*dm2], disp); 

                axa[ip*dax] += grav_const * prtprt_force[0];
                aya[ip*day] += grav_const * prtprt_force[1];
                aza[ip*daz] += grav_const * prtprt_force[2];

                // CkPrintf("particle->particle force (nonperiodic): %f, %f, %f\n", 
                //       grav_const * prtprt_force[0], grav_const * prtprt_force[1], grav_const * prtprt_force[2]);

                tot_prt_force[0] += grav_const * prtprt_force[0];
                tot_prt_force[1] += grav_const * prtprt_force[1];
                tot_prt_force[2] += grav_const * prtprt_force[2];
              }

            }
          }
        }

        // CkPrintf("TOTAL PARTICLE FORCE: %.12f, %.12f, %.12f\n", tot_prt_force[0], tot_prt_force[1], tot_prt_force[2]);
        //CkPrintf("Accel (prt): %f, %f, %f\n\n", axa[ip*dax], aya[ip*day], aza[ip*daz]);

      }
    }
  }

} 

/************************************************************************/

/// Dual Tree Walk helper functions

void EnzoBlock::p_method_multipole_traverse(Index index, int type)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->traverse(this, index, type);
}


void EnzoMethodMultipole::traverse
(EnzoBlock * enzo_block, Index index_b, int type_b)
{

  Index index_a = enzo_block->index();
  const bool known_leaf_b = (type_b != -1);
  const bool is_leaf_a = enzo_block->is_leaf();

  if (! known_leaf_b) {
    // If unknown whether B is leaf or not, flip arguments and send to B
    enzo::block_array()[index_b].p_method_multipole_traverse
      (index_a, is_leaf_a ? 1 : 0);
    return;
  }

  // Determine if blocks are "far", and save radii for later
  double ra,rb;
  bool mac = is_far_(enzo_block,index_b,&ra,&rb);

  // is_leaf values: 0 no, 1 yes -1 unknown
  const bool is_leaf_b = (type_b == 1);

  ItChild it_child_a(cello::rank());
  ItChild it_child_b(cello::rank());
  int ica3[3],icb3[3];

  int level_a = index_a.level();
  int level_b = index_b.level();
  int volume_a = block_volume_[level_a];
  int volume_b = block_volume_[level_b];


  // this stalls if max_level > 0 but there's no refinement ==> fixed?
  if (index_a == index_b) {
    // CkPrintf("root double loop\n");

    // A leaf has no children to pair up. Messages sent to the (not yet
    // existing) children of a leaf are held by Charm++ until they are
    // created by a later refinement, and then delivered to them.
    if (is_leaf_a) return;

    // loop 1 through A children
    while (it_child_a.next(ica3)) {
      int ica = ica3[0] + 2*(ica3[1] + 2*ica3[2]);
      Index index_a_child = index_a.index_child(ica3, 0);

      // loop 2 through A children
      while (it_child_b.next(icb3)) {
        int icb = icb3[0] + 2*(icb3[1] + 2*icb3[2]);
        Index index_b_child = index_b.index_child(icb3, 0);
        // avoid calling both traverse (A,B) and traverse (B,A)

        // CkPrintf("%f, %f\n", ica, icb);

        if (ica <= icb) {  
          // call traverse on child 1 and child 2
          enzo::block_array()[index_a_child].p_method_multipole_traverse
            (index_b_child,-1);
        }
      }
    }
  }

  
  else if (mac) {

    // mac satisfied -- interact approximately
    traverse_approx_pair (enzo_block,index_a,volume_a,index_b,volume_b);

  } 


  // should this go before or after the MAC check?
  else if (is_leaf_a && is_leaf_b) {

    // two leafs -- interact directly
    traverse_direct_pair (enzo_block,index_a,volume_a,index_b,volume_b);

  }
  
  
  else if (is_leaf_a) {

    while (it_child_b.next(icb3)) {

      Index index_b_child = index_b.index_child(icb3, 0);
      traverse(enzo_block,index_b_child,-1);
    }

  } 
  
  else if (is_leaf_b) {

    while (it_child_a.next(ica3)) {

      Index index_a_child = index_a.index_child(ica3, 0);
      enzo::block_array()[index_a_child].p_method_multipole_traverse
        (index_b,type_b);
    }
  } 

  else if (ra < rb) { // open the larger block (b)

    while (it_child_b.next(icb3)) {
      Index index_b_child = index_b.index_child(icb3, 0);
      traverse(enzo_block,index_b_child,-1);
    }
  } 
  
  else {  // open the larger block (a)

    while (it_child_a.next(ica3)) {
      Index index_a_child = index_a.index_child(ica3, 0);
      enzo::block_array()[index_a_child].p_method_multipole_traverse
        (index_b,type_b);   
    }
  }
}


void EnzoMethodMultipole::traverse_approx_pair
(EnzoBlock * enzo_block,
 Index index_a, int volume_a,
 Index index_b, int volume_b)
{

  // update volumes for each block to detect when to terminate traverse.
  // The volume of each Block is carried by the interaction message to the
  // other one, which credits it after applying the interaction (see
  // p_method_multipole_interact_approx()), so that a Block cannot be counted as
  // complete before the interaction has been applied. (Messages between
  // two Blocks are not guaranteed to arrive in the order they were sent.)

  // compute the a->b interaction
  MultipoleMsg * msg_a = pack_multipole_(enzo_block);
  // (note blocks may be equal, so don't double-count)
  msg_a->volume = (index_a != index_b) ? volume_a : 0;
  enzo::block_array()[index_b].p_method_multipole_interact_approx(msg_a);

  // compute the b->a interaction; can we exploit symmetry here?
  enzo::block_array()[index_b].p_method_multipole_interact_approx_send (index_a);

  if (index_a == index_b) {
    update_volume (enzo_block,index_b,volume_b);
  }
}

void EnzoBlock::p_method_multipole_interact_approx (MultipoleMsg * msg)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  const int volume = msg->volume;
  method->interact_approx_ (this, msg);
  if (volume != 0) {
    method->update_volume (this, this->index(), volume);
  }
}

void EnzoBlock::p_method_multipole_interact_approx_send (Index receiver)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->interact_approx_send (this, receiver);
}

void EnzoMethodMultipole::interact_approx_send(EnzoBlock * enzo_block, Index receiver) throw()
{
  MultipoleMsg * msg = pack_multipole_(enzo_block);
  // credited by the receiver after it has applied the multipole
  msg->volume = block_volume_[enzo_block->level()];
  enzo::block_array()[receiver].p_method_multipole_interact_approx(msg);
}


void EnzoMethodMultipole::interact_approx_(Block * block, MultipoleMsg * msg_b) throw()
{
  // retrieve data from this block
  double * pcom_a = pcom(block);
  double * pc1_a = pc1(block);
  double * pc2_a = pc2(block);
  double * pc3_a = pc3(block);
  
  std::array<double, 3> com_a;
  std::array<double, 3> c1_a;
  std::array<double, 6> c2_a;
  std::array<double, 10> c3_a;

  std::copy(pcom_a, pcom_a+3, com_a.begin());
  std::copy(pc1_a, pc1_a+3, c1_a.begin());
  std::copy(pc2_a, pc2_a+6, c2_a.begin());
  std::copy(pc3_a, pc3_a+10, c3_a.begin());

  // retrieve data from msg_b
  double mass_b = msg_b->mass;
  std::array<double, 3> com_b; 
  std::array<double, 6> quadrupole_b;

  std::copy(msg_b->com, msg_b->com+3, com_b.begin());
  std::copy(msg_b->quadrupole, msg_b->quadrupole+6, quadrupole_b.begin());
  

  std::array<double, 3> rvec;

  if (ewald_ != nullptr) {
    double b_image[3];
    
    // should i be taking the nearest periodic image here? yes for ewald contribution, but what about 
    // for the non-periodic d2 and d3? does it matter, since ewald will capture all other periodic images?
    cello::hierarchy()->get_nearest_periodic_image(com_b.data(), com_a.data(), b_image);
    
    const int rank = cello::rank();
    rvec[0] = b_image[0] - com_a[0];
    rvec[1] = (rank < 2) ? 0 : b_image[1] - com_a[1];
    rvec[2] = (rank < 3) ? 0 : b_image[2] - com_a[2];

    //CkPrintf("rvec in interact_approx: %f, %f, %f\n", rvec[0], rvec[1], rvec[2]);
  }
  else {
    rvec = subtract_(com_b, com_a);        // displacement vector between com_b and com_a
  }
  
  double r2 = dot_11_(rvec, rvec);                          // magnitude of displacement vector
  double r3 = r2 * sqrt(r2);
  double r5 = r2 * r3;
  double r7 = r2 * r5;
  
  std::array<double, 3> d1 = dot_scalar_1_(-1.0/r3, rvec);  // derivative tensor d1
  std::array<double, 6> d2;                                 // derivative tensor d2
  std::array<double, 10> d3;                                // derivative tensor d3
          
  // compute the components of d2
  d2[0] = 3.0/r5 * rvec[0] * rvec[0] - 1.0/r3;
  d2[1] = 3.0/r5 * rvec[0] * rvec[1];
  d2[2] = 3.0/r5 * rvec[0] * rvec[2];
  d2[3] = 3.0/r5 * rvec[1] * rvec[1] - 1.0/r3;
  d2[4] = 3.0/r5 * rvec[1] * rvec[2];
  d2[5] = 3.0/r5 * rvec[2] * rvec[2] - 1.0/r3;

  // compute the components of d3
  // precompute coefficients like -15/r7
  d3[0] = -15.0/r7 * rvec[0] * rvec[0] * rvec[0] + 9.0/r5 * rvec[0];
  d3[1] = -15.0/r7 * rvec[0] * rvec[0] * rvec[1] + 3.0/r5 * rvec[1];
  d3[2] = -15.0/r7 * rvec[0] * rvec[0] * rvec[2] + 3.0/r5 * rvec[2];
  d3[3] = -15.0/r7 * rvec[0] * rvec[1] * rvec[1] + 3.0/r5 * rvec[0];
  d3[4] = -15.0/r7 * rvec[0] * rvec[1] * rvec[2];
  d3[5] = -15.0/r7 * rvec[0] * rvec[2] * rvec[2] + 3.0/r5 * rvec[0];
  d3[6] = -15.0/r7 * rvec[1] * rvec[1] * rvec[1] + 9.0/r5 * rvec[1];
  d3[7] = -15.0/r7 * rvec[1] * rvec[1] * rvec[2] + 3.0/r5 * rvec[2];
  d3[8] = -15.0/r7 * rvec[1] * rvec[2] * rvec[2] + 3.0/r5 * rvec[1];
  d3[9] = -15.0/r7 * rvec[2] * rvec[2] * rvec[2] + 9.0/r5 * rvec[2];
      

  if (ewald_ != nullptr) {

    // CkPrintf("ewald in interact_approx\n");

    std::array<double, 3> d1_ewald = ewald_->interp_d1(rvec[0], rvec[1], rvec[2]); // d1 contribution from periodicity
    std::array<double, 6> d2_ewald = ewald_->interp_d2(rvec[0], rvec[1], rvec[2]); // d2 contribution from periodicity
    std::array<double, 10> d3_ewald = ewald_->interp_d3(rvec[0], rvec[1], rvec[2]); // d3 contribution from periodicity
      
    //CkPrintf("d1_ewald in interact_approx: %f, %f, %f\n", d1_ewald[0], d1_ewald[1], d1_ewald[2]);
    
    // Combine the derivative tensors from the periodic and non-periodic contributions    
    // The tables hold the derivatives of n(W) = 1/|W| - V_per(W), so those
    // of the periodic Green's function V_per are the Newtonian ones minus
    // the tabulated ones (the direct sums add +grad n to the acceleration
    // for the same reason).
    d1 = add_11_(d1, dot_scalar_1_(-1.0, d1_ewald));
    d2 = add_22_(d2, dot_scalar_2_(-1.0, d2_ewald));
    d3 = add_33_(d3, dot_scalar_3_(-1.0, d3_ewald));
  }
  
  // compute the coefficients of the Taylor expansion of acceleration due to the particles in Block b
  std::array<double, 3> delta_c1 = add_11_(dot_scalar_1_(mass_b, d1), dot_23_(quadrupole_b, dot_scalar_3_(0.5, d3)));
  std::array<double, 6> delta_c2 = dot_scalar_2_(mass_b, d2);
  std::array<double, 10> delta_c3 = dot_scalar_3_(mass_b, d3);
   
  // add the coefficients for the new interaction to the coefficients already associated with this Block
  std::array<double, 3> new_c1 = add_11_(c1_a, delta_c1);  
  std::array<double, 6> new_c2 = add_22_(c2_a, delta_c2);
  std::array<double, 10> new_c3 = add_33_(c3_a, delta_c3);

  // unravel?
  for (int i = 0; i < 3; i++) {
    pc1(block)[i] = new_c1[i];
  }

  for (int i = 0; i < 6; i++) {
    pc2(block)[i] = new_c2[i];
  }

  for (int i = 0; i < 10; i++) {
    pc3(block)[i] = new_c3[i];
  }

}



void EnzoMethodMultipole::traverse_direct_pair
(EnzoBlock * enzo_block,
 Index index_a, int volume_a,
 Index index_b, int volume_b)
{
  // update volumes for each block to detect when to terminate traverse.
  // The volume is carried by the interaction message and credited by the
  // receiver after applying it (see traverse_approx_pair()).

  // compute the a->b interaction
  // (note blocks may be equal, so don't double-count)
  pack_dens_(enzo_block, index_b, (index_a != index_b) ? volume_a : 0);

  // compute the b->a interaction; any way to exploit symmetry here?
  enzo::block_array()[index_b].p_method_multipole_interact_direct_send (index_a);

  if (index_a == index_b) {
    update_volume (enzo_block,index_b,volume_b);
  }

}

void EnzoBlock::p_method_multipole_interact_direct (int volume, int fldsize, int prtsize, char * fldbuffer, char * prtbuffer)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->interact_direct_ (this, fldbuffer, prtbuffer);
  if (volume != 0) {
    method->update_volume (this, this->index(), volume);
  }
}

void EnzoBlock::p_method_multipole_interact_direct_send (Index receiver)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->interact_direct_send (this, receiver);
}

void EnzoMethodMultipole::interact_direct_send(EnzoBlock * enzo_block, Index receiver) throw()
{
  // the receiver credits this volume after applying the data
  pack_dens_ (enzo_block, receiver, block_volume_[enzo_block->level()]);
}



void EnzoMethodMultipole::pack_dens_(EnzoBlock * enzo_block, Index index_b, int volume) throw()
{

  Data * data = enzo_block->data();
  Field field = data->field();
  Particle particle = data->particle();

  enzo_float * dens = (enzo_float*) field.values ("density");

  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;
  double lo[3];

  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);

  enzo_block->cell_width(&hx, &hy, &hz);

  enzo_block->lower(lo, lo+1, lo+2);

  double * mass = pmass(enzo_block);
  double * com = pcom(enzo_block);
  double * quadrupole = pquadrupole(enzo_block);

  int fldsize = 0;

  SIZE_SCALAR_TYPE(fldsize, int, mx);
  SIZE_SCALAR_TYPE(fldsize, int, my);
  SIZE_SCALAR_TYPE(fldsize, int, mz);
  SIZE_SCALAR_TYPE(fldsize, double, hx);
  SIZE_SCALAR_TYPE(fldsize, double, hy);
  SIZE_SCALAR_TYPE(fldsize, double, hz);
  SIZE_ARRAY_TYPE(fldsize, double, lo, 3);
  SIZE_ARRAY_TYPE(fldsize, enzo_float, dens, mx*my*mz);

  //CkPrintf("is ewald a nullptr in pack_dens? %d\n", ewald_ != nullptr);
  if (ewald_ != nullptr) {
    SIZE_SCALAR_TYPE(fldsize, double, mass);
    SIZE_ARRAY_TYPE(fldsize, double, com, 3);
    SIZE_ARRAY_TYPE(fldsize, double, quadrupole, 6);
  }

  char * fldbuffer = new char[fldsize];

  char * pc; 
  pc = fldbuffer;

  SAVE_SCALAR_TYPE(pc, int, mx);
  SAVE_SCALAR_TYPE(pc, int, my);
  SAVE_SCALAR_TYPE(pc, int, mz);
  SAVE_SCALAR_TYPE(pc, double, hx);
  SAVE_SCALAR_TYPE(pc, double, hy);
  SAVE_SCALAR_TYPE(pc, double, hz);
  SAVE_ARRAY_TYPE(pc, double, lo, 3);
  SAVE_ARRAY_TYPE(pc, enzo_float, dens, mx*my*mz);

  if (ewald_ != nullptr) {
    SAVE_SCALAR_TYPE(pc, double, mass);
    SAVE_ARRAY_TYPE(pc, double, com, 3);
    SAVE_ARRAY_TYPE(pc, double, quadrupole, 6);
  }

  int prtsize = particle.data_size();
  char * prtbuffer = new char[prtsize];
  char * prtbuffer_next = particle.save_data(prtbuffer);

  enzo::block_array()[index_b].p_method_multipole_interact_direct (volume, fldsize, prtsize, fldbuffer, prtbuffer);

  // the entry method call copies the arrays into its message
  delete [] fldbuffer;
  delete [] prtbuffer;
}


void EnzoMethodMultipole::interact_direct_(Block * block, char * fldbuffer_b, char * prtbuffer_b) throw()
{

  Data * data = block->data();
  Field field = data->field();

  enzo_float * accel_x = (enzo_float*) field.values ("acceleration_x");
  enzo_float * accel_y = (enzo_float*) field.values ("acceleration_y");
  enzo_float * accel_z = (enzo_float*) field.values ("acceleration_z");

  int mx, my, mz;
  int gx, gy, gz;
  double hx, hy, hz;
  
  field.dimensions  (0, &mx, &my, &mz);
  field.ghost_depth (0, &gx, &gy, &gz);
  

  block->cell_width(&hx, &hy, &hz);
  // hz = 1.0; // this is to make math easier -- delete later
  double cell_vol = hx * hy * hz;

  double lo[3];
  block->lower(lo, lo+1, lo+2);

  // unpacking data from fldbuffer_b
  int mx2, my2, mz2;
  double hx2, hy2, hz2;
  double lo2[3];

  char * pc;
  pc = fldbuffer_b;

  LOAD_SCALAR_TYPE(pc, int, mx2);
  LOAD_SCALAR_TYPE(pc, int, my2);
  LOAD_SCALAR_TYPE(pc, int, mz2);
  LOAD_SCALAR_TYPE(pc, double, hx2);
  LOAD_SCALAR_TYPE(pc, double, hy2);
  LOAD_SCALAR_TYPE(pc, double, hz2);
  LOAD_ARRAY_TYPE(pc, double, lo2, 3);

  double dens[mx2*my2*mz2];
  LOAD_ARRAY_TYPE(pc, double, dens, mx2*my2*mz2);


  Particle particle = data->particle();
  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();
  const int num_is_grav = particle_groups->size("is_gravitating");
  const int num_prtl_types = particle_descr->num_types();

  enzo_float * prtmass2 = NULL;
  int dm2;

  // unpacking particles from prtbuffer_b
  ParticleData new_p_data;
  Particle particle2 (particle_descr, &new_p_data);
  char * prtbuffer_next = particle2.load_data(prtbuffer_b);

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  double grav_const; 
  if (cosmology) {
    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    double dt = block->dt();
    double time = block->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    grav_const = 1./(4 * cello::pi * cosmo_a * cosmo_a);
  }
  else {
    grav_const = enzo::grav_constant_codeU();
  }


  // loop over all cells in this Block
  for (int iz = gz; iz < mz-gz; iz++) {
    for (int iy = gy; iy < my-gy; iy++) {
      for (int ix = gx; ix < mx-gx; ix++) {

        int i = ix + mx*(iy + my*iz);

	// CkPrintf("position: %f, %f, %f\n", (lo[0] + (ix-gx + 0.5)*hx), (lo[1] + (iy-gy + 0.5)*hy), (lo[2] + (iz-gz + 0.5)*hz));

        std::array<double, 3> tot_cell_force = {}; // for debugging

        // compute force sourced from cells in Block b
        for (int iz2 = gz; iz2 < mz2-gz; iz2++) {
          for (int iy2 = gy; iy2 < my2-gy; iy2++) {
            for (int ix2 = gx; ix2 < mx2-gx; ix2++) {

              int i2 = ix2 + mx2*(iy2 + my2*iz2);

              // CkPrintf("position2: %f, %f, %f\n", (lo2[0] + (ix2-gx + 0.5)*hx2), (lo2[1] + (iy2-gy + 0.5)*hy2), (lo2[2] + (iz2-gz + 0.5)*hz2));

              // disp points from cell in current Block to cell in Block b
              std::array<double, 3> disp;

              double b_mass = dens[i2]*hx2*hy2*hz2;
                
              // force from periodic images of Cell b
              if (ewald_ != nullptr) {

                double a_loc[3]{lo[0] + (ix-gx + 0.5)*hx, lo[1] + (iy-gy + 0.5)*hy, lo[2] + (iz-gz + 0.5)*hz};
                double b_loc[3]{lo2[0] + (ix2-gx + 0.5)*hx2, lo2[1] + (iy2-gy + 0.5)*hy2, lo2[2] + (iz2-gz + 0.5)*hz2};
                double b_image[3];
                cello::hierarchy()->get_nearest_periodic_image(b_loc, a_loc, b_image);

                // CkPrintf("a_loc: %f, %f, %f\n", a_loc[0], a_loc[1], a_loc[2]);
                // CkPrintf("b_loc: %f, %f, %f\n", b_loc[0], b_loc[1], b_loc[2]);
                // CkPrintf("b_image: %f, %f, %f\n", b_image[0], b_image[1], b_image[2]);
                
                // displacement vector pointing from current cell to interacting cell
                const int rank = cello::rank();
                disp[0] = b_image[0] - a_loc[0];
                disp[1] = (rank < 2) ? 0 : b_image[1] - a_loc[1];
                disp[2] = (rank < 3) ? 0 : b_image[2] - a_loc[2];
                  
                // acceleration contribution from periodic images of Cell b
                std::array<double, 3> d1_ewald = ewald_->interp_d1(disp[0], disp[1], disp[2]); 
                  
                // should i be subtracting or adding?
                // minus sign from negative gradient * minus sign from chain rule
                accel_x[i] += grav_const * b_mass * d1_ewald[0];
                accel_y[i] += grav_const * b_mass * d1_ewald[1];
                accel_z[i] += grav_const * b_mass * d1_ewald[2];

                tot_cell_force[0] += grav_const * b_mass * d1_ewald[0];
                tot_cell_force[1] += grav_const * b_mass * d1_ewald[1];
                tot_cell_force[2] += grav_const * b_mass * d1_ewald[2];
              }
              else {
                disp[0] = lo2[0] - lo[0] + (ix2-gx + 0.5)*hx2 - (ix-gx + 0.5)*hx; 
                disp[1] = lo2[1] - lo[1] + (iy2-gy + 0.5)*hy2 - (iy-gy + 0.5)*hy;
                disp[2] = lo2[2] - lo[2] + (iz2-gz + 0.5)*hz2 - (iz-gz + 0.5)*hz;
              }
              

              // CkPrintf("disp: %f, %f, %f\n", disp[0], disp[1], disp[2]);
                
              std::array<double, 3> cell_force = newton_force_(b_mass, disp); 

              accel_x[i] += grav_const * cell_force[0];
              accel_y[i] += grav_const * cell_force[1];
              accel_z[i] += grav_const * cell_force[2];

              tot_cell_force[0] += grav_const * cell_force[0];
              tot_cell_force[1] += grav_const * cell_force[1];
              tot_cell_force[2] += grav_const * cell_force[2];
            }
          }
        }
	      // CkPrintf("Cell loc: %f, %f, %f\n", lo[0] + (ix-gx + 0.5)*hx, lo[1] + (iy-gy + 0.5)*hy, lo[2] + (iz-gz + 0.5)*hz);
        // CkPrintf("Leaf-leaf cell force: %f, %f, %f\n", tot_cell_force[0], tot_cell_force[1], tot_cell_force[2]);
        //CkPrintf("Acceleration: %f, %f, %f\n\n", accel_x[i], accel_y[i], accel_z[i]);

        // compute cell force sourced from particles in Block b
        for (int ipt = 0; ipt < num_is_grav; ipt++) {

          const int it = particle2.type_index(particle_groups->item("is_gravitating",ipt));

          int imass2 = 0;

          for (int ib = 0; ib < particle2.num_batches(it); ib++) {

            const int np = particle2.num_particles(it,ib);

            if (particle2.has_attribute(it,"mass")) {
              imass2 = particle2.attribute_index(it,"mass");
              prtmass2 = (enzo_float *) particle2.attribute_array( it, imass2, ib);
              dm2 = particle2.stride(it,imass2);
            } 
            else {
              imass2 = particle2.constant_index(it,"mass");
              prtmass2 = (enzo_float*)particle2.constant_value(it,imass2);
              dm2 = 0;
            }

            const int ia_x2  = particle2.attribute_index(it,"x");
            const int ia_y2  = particle2.attribute_index(it,"y");
            const int ia_z2  = particle2.attribute_index(it,"z");

            enzo_float * xa2 =  (enzo_float *)particle2.attribute_array (it,ia_x2,ib);
            enzo_float * ya2 =  (enzo_float *)particle2.attribute_array (it,ia_y2,ib);
            enzo_float * za2 =  (enzo_float *)particle2.attribute_array (it,ia_z2,ib);

            const int dx2 =  particle2.stride(it,ia_x2);
            const int dy2 =  particle2.stride(it,ia_y2);
            const int dz2 =  particle2.stride(it,ia_z2);

            for (int ip=0; ip < np; ip++) {
              
              // disp points from cell in current Block to particle in Block b
              std::array<double, 3> disp;
              
              // force from periodic images of Particle b
              if (ewald_ != nullptr) {
                
                double a_loc[3]{lo[0] + (ix-gx + 0.5)*hx, lo[1] + (iy-gy + 0.5)*hy, lo[2] + (iz-gz + 0.5)*hz};
                double b_loc[3]{xa2[ip*dx2], ya2[ip*dy2], za2[ip*dz2]};
                double b_image[3];
                cello::hierarchy()->get_nearest_periodic_image(b_loc, a_loc, b_image);

                const int rank = cello::rank();
                disp[0] = b_image[0] - a_loc[0];
                disp[1] = (rank < 2) ? 0 : b_image[1] - a_loc[1];
                disp[2] = (rank < 3) ? 0 : b_image[2] - a_loc[2];
                  
                // acceleration contribution from periodic images of Particle b
                std::array<double, 3> d1_ewald = ewald_->interp_d1(disp[0], disp[1], disp[2]); 
                  
                accel_x[i] += grav_const * prtmass2[ip*dm2] * d1_ewald[0];
                accel_y[i] += grav_const * prtmass2[ip*dm2] * d1_ewald[1];
                accel_z[i] += grav_const * prtmass2[ip*dm2] * d1_ewald[2];
              }
              else {
                disp[0] = xa2[ip*dx2] - (lo[0] + (ix-gx + 0.5)*hx); 
                disp[1] = ya2[ip*dy2] - (lo[1] + (iy-gy + 0.5)*hy);
                disp[2] = za2[ip*dz2] - (lo[2] + (iz-gz + 0.5)*hz);
              }
                
              std::array<double, 3> prtcell_force = newton_force_(prtmass2[ip*dm2], disp); 

              accel_x[i] += grav_const * prtcell_force[0];
              accel_y[i] += grav_const * prtcell_force[1];
              accel_z[i] += grav_const * prtcell_force[2];

            }
          }
        }
      }
    }
  }


  // loop over all particles in this Block
  for (int it = 0; it < num_prtl_types; it++) {

    for (int ib = 0; ib < particle.num_batches(it); ib++) {

      const int np = particle.num_particles(it,ib);

      const int ia_x  = particle.attribute_index(it,"x");
      const int ia_y  = particle.attribute_index(it,"y");
      const int ia_z  = particle.attribute_index(it,"z");
      const int ia_ax  = particle.attribute_index(it,"ax");
      const int ia_ay  = particle.attribute_index(it,"ay");
      const int ia_az  = particle.attribute_index(it,"az");

      enzo_float * xa =  (enzo_float *)particle.attribute_array (it,ia_x,ib);
      enzo_float * ya =  (enzo_float *)particle.attribute_array (it,ia_y,ib);
      enzo_float * za =  (enzo_float *)particle.attribute_array (it,ia_z,ib);
      enzo_float * axa =  (enzo_float *)particle.attribute_array (it,ia_ax,ib);
      enzo_float * aya =  (enzo_float *)particle.attribute_array (it,ia_ay,ib);
      enzo_float * aza =  (enzo_float *)particle.attribute_array (it,ia_az,ib);

      const int dx =  particle.stride(it,ia_x);
      const int dy =  particle.stride(it,ia_y);
      const int dz =  particle.stride(it,ia_z);
      const int dax =  particle.stride(it,ia_ax);
      const int day =  particle.stride(it,ia_ay);
      const int daz =  particle.stride(it,ia_az);

      for (int ip=0; ip < np; ip++) {
        
        // compute force sourced by cells in Block b
        for (int iz2 = gz; iz2 < mz2-gz; iz2++) {
          for (int iy2 = gy; iy2 < my2-gy; iy2++) {
            for (int ix2 = gx; ix2 < mx2-gx; ix2++) {

              int i2 = ix2 + mx2*(iy2 + my2*iz2);

              // disp points from particle in current Block to cell in Block b
              std::array<double, 3> disp;

              double b_mass = dens[i2]*hx2*hy2*hz2;

              // force from periodic images of Cell b
              if (ewald_ != nullptr) {

                double a_loc[3]{xa[ip*dx], ya[ip*dy], za[ip*dz]};
                double b_loc[3]{lo2[0] + (ix2-gx + 0.5)*hx2, lo2[1] + (iy2-gy + 0.5)*hy2, lo2[2] + (iz2-gz + 0.5)*hz2};
                double b_image[3];
                cello::hierarchy()->get_nearest_periodic_image(b_loc, a_loc, b_image);

                const int rank = cello::rank();
                disp[0] = b_image[0] - a_loc[0];
                disp[1] = (rank < 2) ? 0 : b_image[1] - a_loc[1];
                disp[2] = (rank < 3) ? 0 : b_image[2] - a_loc[2];
                  
                // acceleration contribution from periodic images of Cell b
                std::array<double, 3> d1_ewald = ewald_->interp_d1(disp[0], disp[1], disp[2]); 
                  
                axa[ip*dax] += grav_const * b_mass * d1_ewald[0];
                aya[ip*day] += grav_const * b_mass * d1_ewald[1];
                aza[ip*daz] += grav_const * b_mass * d1_ewald[2];
              }
              else {
                disp[0] = (lo2[0] + (ix2-gx + 0.5)*hx2) - xa[ip*dx]; 
                disp[1] = (lo2[1] + (iy2-gy + 0.5)*hy2) - ya[ip*dy];
                disp[2] = (lo2[2] + (iz2-gz + 0.5)*hz2) - za[ip*dz];
              }

              std::array<double, 3> cellprt_force = newton_force_(b_mass, disp); 

              axa[ip*dax] += grav_const * cellprt_force[0];
              aya[ip*day] += grav_const * cellprt_force[1];
              aza[ip*daz] += grav_const * cellprt_force[2];

            }
          }
        }

        // compute force sourced by particles in Block b
        for (int ipt2 = 0; ipt2 < num_is_grav; ipt2++) {
          const int it2 = particle2.type_index(particle_groups->item("is_gravitating",ipt2));

          int imass2 = 0;

          for (int ib2 = 0; ib2 < particle2.num_batches(it2); ib2++) {

            const int np2 = particle2.num_particles(it2,ib2);

            if (particle2.has_attribute(it2,"mass")) {
              imass2 = particle2.attribute_index(it2,"mass");
              prtmass2 = (enzo_float *) particle2.attribute_array( it2, imass2, ib2);
              dm2 = particle2.stride(it2,imass2);
            } 
            else {
              imass2 = particle2.constant_index(it2,"mass");
              prtmass2 = (enzo_float*)particle2.constant_value(it2,imass2);
              dm2 = 0;
            }

            const int ia_x2  = particle2.attribute_index(it2,"x");
            const int ia_y2  = particle2.attribute_index(it2,"y");
            const int ia_z2  = particle2.attribute_index(it2,"z");

            enzo_float * xa2 =  (enzo_float *)particle2.attribute_array (it2,ia_x2,ib2);
            enzo_float * ya2 =  (enzo_float *)particle2.attribute_array (it2,ia_y2,ib2);
            enzo_float * za2 =  (enzo_float *)particle2.attribute_array (it2,ia_z2,ib2);

            const int dx2 =  particle2.stride(it2,ia_x2);
            const int dy2 =  particle2.stride(it2,ia_y2);
            const int dz2 =  particle2.stride(it2,ia_z2);

            for (int ip2=0; ip2 < np2; ip2++) {
              
              // disp points from particle in current Block to particle in Block b
              std::array<double, 3> disp;

              // force from periodic images of Particle b
              if (ewald_ != nullptr) {

                double a_loc[3]{xa[ip*dx], ya[ip*dy], za[ip*dz]};
                double b_loc[3]{xa2[ip2*dx2], ya2[ip2*dy2], za2[ip2*dz2]};
                double b_image[3];
                cello::hierarchy()->get_nearest_periodic_image(b_loc, a_loc, b_image);

                const int rank = cello::rank();
                disp[0] = b_image[0] - a_loc[0];
                disp[1] = (rank < 2) ? 0 : b_image[1] - a_loc[1];
                disp[2] = (rank < 3) ? 0 : b_image[2] - a_loc[2];
                  
                // acceleration contribution from periodic images of Particle b
                std::array<double, 3> d1_ewald = ewald_->interp_d1(disp[0], disp[1], disp[2]); 
                  
                axa[ip*dax] += grav_const * prtmass2[ip2*dm2] * d1_ewald[0];
                aya[ip*day] += grav_const * prtmass2[ip2*dm2] * d1_ewald[1];
                aza[ip*daz] += grav_const * prtmass2[ip2*dm2] * d1_ewald[2];
              }
              else {
                disp[0] = xa2[ip2*dx2] - xa[ip*dx]; 
                disp[1] = ya2[ip2*dy2] - ya[ip*dy];
                disp[2] = za2[ip2*dz2] - za[ip*dz];
              }

              std::array<double, 3> prtprt_force = newton_force_(prtmass2[ip2*dm2], disp); 

              axa[ip*dax] += grav_const * prtprt_force[0];
              aya[ip*day] += grav_const * prtprt_force[1];
              aza[ip*daz] += grav_const * prtprt_force[2];

            }
          }
        }

      }
    }
  }

  // compute long-range contribution from periodic images of Block b
  /* 
  if (ewald_ != nullptr) {

    double * com_a = pcom(block);
    double * pc1_a = pc1(block);
    double * pc2_a = pc2(block);
    double * pc3_a = pc3(block);
    
    std::array<double, 3> c1_a;
    std::array<double, 6> c2_a;
    std::array<double, 10> c3_a;

    std::copy(pc1_a, pc1_a+3, c1_a.begin());
    std::copy(pc2_a, pc2_a+6, c2_a.begin());
    std::copy(pc3_a, pc3_a+10, c3_a.begin());

    double mass_b;
    double com_b[3];
    std::array<double, 6> quadrupole_b;
    double quadrupole_array[6];

    // CkPrintf("is the segfault before load scalar?\n");
    LOAD_SCALAR_TYPE(pc, double, mass_b);
    LOAD_ARRAY_TYPE(pc, double, com_b, 3);
    LOAD_ARRAY_TYPE(pc, double, quadrupole_array, 6);
    // CkPrintf("is the segfault after load scalar?\n");

    // is this necessary?
    for (int i=0; i<6; i++) {
      quadrupole_b[i] = quadrupole_array[i];
    }

    // replace rvec with scalars?
    double b_image[3]; 
    double rvec[3];
    cello::hierarchy()->get_nearest_periodic_image(com_b, com_a, b_image);
    CkPrintf("com_a: %f %f %f\n", com_a[0], com_a[1], com_a[2]);
    CkPrintf("com_b: %f %f %f\n", com_b[0], com_b[1], com_b[2]);
    CkPrintf("b_image: %f %f %f\n", b_image[0], b_image[1], b_image[2]); 
    // b_image[2] = 0.5; 
    // get_nearest_periodic_image has a check for dimensionality which i don't have, making the z-component garbage
    
    const int rank = cello::rank();
    rvec[0] = b_image[0] - com_a[0];
    rvec[1] = (rank < 2) ? 0 : b_image[1] - com_a[1];
    rvec[2] = (rank < 3) ? 0 : b_image[2] - com_a[2];
   
    // CkPrintf("rvec in interact_direct: %f, %f, %f\n", rvec[0], rvec[1], rvec[2]);

    std::array<double, 3> d1_ewald = ewald_->interp_d1(rvec[0], rvec[1], rvec[2]); // d1 contribution from periodicity
    std::array<double, 6> d2_ewald = ewald_->interp_d2(rvec[0], rvec[1], rvec[2]); // d2 contribution from periodicity
    std::array<double, 10> d3_ewald = ewald_->interp_d3(rvec[0], rvec[1], rvec[2]); // d3 contribution from periodicity

    // CkPrintf("d1_ewald in interact_direct: %f, %f, %f\n", d1_ewald[0], d1_ewald[1], d1_ewald[2]);
    // CkPrintf("d2_ewald in interact_direct: %f, %f, %f, %f, %f, %f\n", d2_ewald[0], d2_ewald[1], d2_ewald[2], d2_ewald[3], d2_ewald[4], d2_ewald[5]);
    
    // compute the coefficients of the Taylor expansion of acceleration due to the particles in Block b
    std::array<double, 3> delta_c1 = add_11_(dot_scalar_1_(mass_b, d1_ewald), dot_23_(quadrupole_b, dot_scalar_3_(0.5, d3_ewald)));
    std::array<double, 6> delta_c2 = dot_scalar_2_(mass_b, d2_ewald);
    std::array<double, 10> delta_c3 = dot_scalar_3_(mass_b, d3_ewald);

    //CkPrintf("c1 in interact_direct: %f, %f, %f\n", delta_c1[0], delta_c1[1], delta_c1[2]);
    //CkPrintf("c2 in interact_direct: %f, %f, %f, %f, %f, %f\n", delta_c2[0], delta_c2[1], delta_c2[2], delta_c2[3], delta_c2[4], delta_c2[5]);
    //CkPrintf("c3 in interact_direct: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", delta_c3[0], delta_c3[1], delta_c3[2], delta_c3[3], delta_c3[4], delta_c3[5], delta_c3[6], delta_c3[7], delta_c3[8], delta_c3[9]);
    
    // add the coefficients for the new interaction to the coefficients already associated with this Block
    std::array<double, 3> new_c1 = add_11_(c1_a, delta_c1);  
    std::array<double, 6> new_c2 = add_22_(c2_a, delta_c2);
    std::array<double, 10> new_c3 = add_33_(c3_a, delta_c3);

    for (int i = 0; i < 3; i++) {
      pc1(block)[i] = new_c1[i];
    }

    for (int i = 0; i < 6; i++) {
      pc2(block)[i] = new_c2[i];
    }

    for (int i = 0; i < 10; i++) {
      pc3(block)[i] = new_c3[i];
    }
  }
  */

}


void EnzoBlock::p_method_multipole_update_volume (Index index, int volume)
{
  EnzoMethodMultipole * method = static_cast<EnzoMethodMultipole*> (this->method());
  method->update_volume (this, index, volume);
}

// i don't think we need the index?
void EnzoMethodMultipole::update_volume
(EnzoBlock * enzo_block, Index index, int volume)
{

  std::string string = enzo_block->name();
   

  if (enzo_block->is_leaf()) {

    *volume_(enzo_block) += volume;

    // max_volume_ plus the unit credited by this Block's own dual_walk_()
    if (*volume_(enzo_block) == max_volume_ + 1) {

      CkCallback callback(CkIndex_EnzoBlock::r_method_multipole_traverse_complete(nullptr),
			        enzo::block_array());
      enzo_block->contribute(callback);

    }

  } 
  
  else {

    const int min_level = 0; //cello::hierarchy()->min_level();
    ItChild it_child(cello::rank());
    int ic3[3];

    while (it_child.next(ic3)) {

      Index index_child = enzo_block->index().index_child(ic3,min_level);
      enzo::block_array()[index_child].p_method_multipole_update_volume (index,volume);
    }
  }
}

bool EnzoMethodMultipole::is_far_ (EnzoBlock * enzo_block,
                             Index index_b, double *ra, double *rb) const
{
  double iam3[3],iap3[3],ibm3[3],ibp3[3],na3[3],nb3[3];
  enzo_block->lower(iam3,iam3+1,iam3+2);
  enzo_block->upper(iap3,iap3+1,iap3+2);
  enzo_block->lower(ibm3,ibm3+1,ibm3+2,&index_b);  // overloaded lower/upper? (see Cello/mesh_Block.hpp)
  enzo_block->upper(ibp3,ibp3+1,ibp3+2,&index_b);
  *ra=0;
  *rb=0;
  double ra3[3]={0,0,0}, rb3[3]={0,0,0};
  double ca[3]={0,0,0}, cb[3]={0,0,0};
  double d3[3]={0,0,0};
  double d = 0;

  for (int i=0; i<cello::rank(); i++) {
    ra3[i] = (iap3[i] - iam3[i]);   // length of side i of Block a
    rb3[i] = (ibp3[i] - ibm3[i]);   // length of side i of Block b
    *ra += ra3[i]*ra3[i];           // square of the length of side i of Block a
    *rb += rb3[i]*rb3[i];           // square of the length of side i of Block b

    ca[i] = (iam3[i]+0.5*ra3[i]);  // i coordinate of center of Block a
    cb[i] = (ibm3[i]+0.5*rb3[i]);  // i coordinate of center of Block b

    // double ca = (iam3[i]+0.5*ra3[i]);
    // double cb = (ibm3[i]+0.5*rb3[i]);
    // d3[i] = (ca - cb); // check if this is the nearest periodic distance
    // d += d3[i]*d3[i];
  }

  // compute the displacement between the centers of the two Blocks
  if (ewald_ != nullptr) {
    double ca_image[3];

    // be sure to check for 1D/2D/3D and set the junk components accordingly
    cello::hierarchy()->get_nearest_periodic_image(ca, cb, ca_image);

    const int rank = cello::rank();
    d3[0] = ca_image[0] - cb[0];
    d3[1] = (rank < 2) ? 0 : ca_image[1] - cb[1];
    d3[2] = (rank < 3) ? 0 : ca_image[2] - cb[2];

  }
  else {
    d3[0] = ca[0] - cb[0];
    d3[1] = ca[1] - cb[1];
    d3[2] = ca[2] - cb[2]; 
  }

  d = sqrt(d3[0]*d3[0] + d3[1]*d3[1] + d3[2]*d3[2]);  // distance between centers of Blocks
  *ra = 0.5*sqrt(*ra);  // half the length of the diagonal of Block a
  *rb = 0.5*sqrt(*rb);  // half the length of the diagonal of Block b 

  if (d * theta_ > (*ra + *rb)) {
    CkPrintf("angle: %f\n", (*ra + *rb)/d);
  }

  return  (d * theta_ > (*ra + *rb));
}

/************************************************************************/

/**** for particle debugging  ****/

// see EnzoInitialIsolatedGalaxy for more detailed initialization

 void EnzoMethodMultipole::InitializeParticles(Block * block, int nprtls, double prtls[][7]){

  Particle particle = block->data()->particle();
  ParticleDescr * particle_descr = cello::particle_descr();

  double lo[3];
  double hi[3];
  block->lower(lo, lo+1, lo+2);
  block->upper(hi, hi+1, hi+2);

  // CkPrintf("lo: %f, %f, %f\n", lo[0], lo[1], lo[2]);
  // CkPrintf("hi: %f, %f, %f\n", hi[0], hi[1], hi[2]);

  //
  // Loop through all particle types and initialize their positions and
  // accelerations.
  //

  int ntypes = 1;
  int * particleIcTypes = new int[ntypes];

  particleIcTypes[0] = particle_descr->type_index("star");
  // particleIcTypes[1] = particle_descr->type_index("trace");

  // particleIcFileNames.push_back("halo.dat");
  // nparticles_ = std::max(nparticles_, nlines("halo.dat"));
  // ipt++;

  // Loop over all particle types and initialize
  for(int ipt = 0; ipt < ntypes; ipt++){

    int it   = particleIcTypes[ipt];

    // obtain particle attribute indexes for this type
    int ia_x = particle.attribute_index (it, "x");
    int ia_y = particle.attribute_index (it, "y");
    int ia_z = particle.attribute_index (it, "z");
    int ia_ax = particle.attribute_index (it, "ax");
    int ia_ay = particle.attribute_index (it, "ay");
    int ia_az = particle.attribute_index (it, "az");

    int ib  = 0; // batch counter
    int ipp = 0; // particle counter

    // this will point to the particular value in the
    // particle attribute array

    enzo_float * px   = 0;
    enzo_float * py   = 0;
    enzo_float * pz   = 0;
    enzo_float * pax  = 0;
    enzo_float * pay  = 0;
    enzo_float * paz  = 0;

    // change to ntypes-1 if you want to test tracer particles
    if (ipt != ntypes) {

      int ia_m = particle.attribute_index (it, "mass");
      enzo_float * prtmass = 0;

      for (int i = 0; i < nprtls; i++){

        // CkPrintf("pos: %f, %f, %f\n", prtls[i][1], prtls[i][2], prtls[i][3]);

        if (prtls[i][1] >= lo[0] && prtls[i][1] < hi[0] &&
            prtls[i][2] >= lo[1] && prtls[i][2] < hi[1] &&
            prtls[i][3] >= lo[2] && prtls[i][3] < hi[2]) {

          // CkPrintf("added a particle!\n");

          int new_particle = particle.insert_particles(it, 1);
          particle.index(new_particle,&ib,&ipp);

          // get pointers to each of the associated arrays
          prtmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
          px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
          py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
          pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
          pax   = (enzo_float *) particle.attribute_array(it, ia_ax, ib);
          pay   = (enzo_float *) particle.attribute_array(it, ia_ay, ib);
          paz   = (enzo_float *) particle.attribute_array(it, ia_az, ib);


          // set the particle values
          prtmass[ipp] = prtls[i][0];
          px[ipp]    = prtls[i][1];
          py[ipp]    = prtls[i][2];
          pz[ipp]    = prtls[i][3];
          pax[ipp]   = prtls[i][4];
          pay[ipp]   = prtls[i][5];
          paz[ipp]   = prtls[i][6];
        }
      }
    } 

    else {

      int i = nprtls - 1;

      if (prtls[i][1] >= lo[0] && prtls[i][1] < hi[0] &&
            prtls[i][2] >= lo[1] && prtls[i][2] < hi[1] &&
            prtls[i][3] >= lo[2] && prtls[i][3] < hi[2]) {

          // CkPrintf("added a particle!\n");

          int new_particle = particle.insert_particles(it, 1);
          particle.index(new_particle,&ib,&ipp);

          // get pointers to each of the associated arrays
          px    = (enzo_float *) particle.attribute_array(it, ia_x, ib);
          py    = (enzo_float *) particle.attribute_array(it, ia_y, ib);
          pz    = (enzo_float *) particle.attribute_array(it, ia_z, ib);
          pax   = (enzo_float *) particle.attribute_array(it, ia_ax, ib);
          pay   = (enzo_float *) particle.attribute_array(it, ia_ay, ib);
          paz   = (enzo_float *) particle.attribute_array(it, ia_az, ib);


          // set the particle values
          px[ipp]    = prtls[i][1];
          py[ipp]    = prtls[i][2];
          pz[ipp]    = prtls[i][3];
          pax[ipp]   = prtls[i][4];
          pay[ipp]   = prtls[i][5];
          paz[ipp]   = prtls[i][6];
      }
    }

  } // end loop over particle types


  return;
}

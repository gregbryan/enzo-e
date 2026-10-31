// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMethodMultipole.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu) 
/// @date     Tuesday July 19, 2022
/// @brief    [\ref Enzo] Declaration of EnzoMethodMultipole
/// Compute accelerations on gas and particles using FMM

#ifndef ENZO_METHOD_MULTIPOLE_HPP
#define ENZO_METHOD_MULTIPOLE_HPP

#include "Enzo/gravity/EnzoEwald.hpp"

class EnzoMethodMultipole : public Method {

  /// @class    EnzoMethodMultipole
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Compute multipoles and pass multipoles up octree

public: // interface -- which methods should be public and which protected?

  /// Create a new EnzoMethodMultipole object
  EnzoMethodMultipole(ParameterGroup p);

  EnzoMethodMultipole() = delete;


  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodMultipole);

  /// Charm++ PUP::able migration constructor
  EnzoMethodMultipole (CkMigrateMessage *m)
    : Method (m),
      i_sync_restrict_(-1),
      //i_sync_prolong_(-1),
      i_msg_restrict_(),
      i_msg_prolong_(-1),
      theta_(0),
      eps0_(0),
      // r0_(0),
      is_volume_(-1),
      block_volume_(),
      max_volume_(0),
      interp_xpoints_(32),
      interp_ypoints_(32),
      interp_zpoints_(32),
      ewald_order_(3),
      ewald_(nullptr),
      dt_max_(1.0e10),
      leaf_chunk_(0),
      leaf_theta_(0.5),
      ir_exit_(-1)
  { for (int i = 0; i < cello::num_children(); i++) i_msg_restrict_[i] = -1; }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  ~EnzoMethodMultipole()
  {
    if (ewald_ != nullptr) {
      delete ewald_;
   }
  }

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "multipole"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) throw();
  virtual double timestep_( Block * block) throw();

  void begin_up_cycle_ (EnzoBlock * enzo_block) throw();
  void begin_down_cycle_ (EnzoBlock * enzo_block) throw();

  void restrict_send (EnzoBlock * enzo_block) throw();
  void prolong_send (EnzoBlock * enzo_block) throw();

  void restrict_recv (EnzoBlock * enzo_block, MultipoleMsg * msg) throw();
  void prolong_recv (EnzoBlock * enzo_block, MultipoleMsg * msg) throw();

  MultipoleMsg * pack_multipole_ (EnzoBlock * enzo_block) throw();
  MultipoleMsg * pack_coeffs_ (EnzoBlock * enzo_block) throw();

  void unpack_multipole_ (EnzoBlock * enzo_block, MultipoleMsg * msg) throw();
  void unpack_coeffs_ (EnzoBlock * enzo_block, MultipoleMsg * msg) throw();

  /*****************************************************/

  /// Access the restrict Sync Scalar value for the Block
  Sync * psync_restrict(Block * block)
  {
    ScalarData<Sync> * scalar_data = block->data()->scalar_data_sync();
    ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
    return scalar_data->value(scalar_descr, i_sync_restrict_);
  }

  /// Access the prolong Sync Scalar value for the Block
  //  Sync * psync_prolong(Block * block)
  // {
  //  ScalarData<Sync> * scalar_data = block->data()->scalar_data_sync();
  //  ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
  //  return scalar_data->value(scalar_descr, i_sync_prolong_);
  //}

  /// Access the Field message for buffering restriction data
  MultipoleMsg ** pmsg_restrict(Block * block, int ic)
  {
    ScalarData<void *> * scalar_data = block->data()->scalar_data_void();
    ScalarDescr *        scalar_descr = cello::scalar_descr_void();
    return (MultipoleMsg **)scalar_data->value(scalar_descr, i_msg_restrict_[ic]);
  }

  /// Access the Field message for buffering prolongation data
  MultipoleMsg ** pmsg_prolong(Block * block)
  {
    ScalarData<void *> * scalar_data = block->data()->scalar_data_void();
    ScalarDescr *        scalar_descr = cello::scalar_descr_void();
    
    void ** data = scalar_data->value(scalar_descr, i_msg_prolong_);

    if (data == nullptr){
      ERROR("data is null", "data is null");
    }

    return (MultipoleMsg **)data;
  }

  /// returns a pointer to this Block's mass
  double * pmass(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_mass_);
  }

  /// returns an array of three elements corresponding to the components of the center-of-mass
  double * pcom(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_com_);
  }

  /// returns a 6-element array corresponding to the components of the quadrupole tensor
  double * pquadrupole(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_quadrupole_);
  }

  /// returns a 3-element array corresponding to the first Taylor coefficient of the long-range force
  double * pc1(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_c1_);
  }

  /// returns a 6-element array corresponding to the second Taylor coefficient of the long-range force
  double * pc2(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_c2_);
  }

  /// returns a 10-element array corresponding to the third Taylor coefficient of the long-range force
  double * pc3(Block * block)
  {
    ScalarData<double> * scalar_data = block->data()->scalar_data_double();
    ScalarDescr *      scalar_descr = cello::scalar_descr_double();
    return scalar_data->value(scalar_descr, i_c3_);
  }

  /*****************************************************/


protected: // methods

  void compute_ (Block * block) throw();


  /**********    FMM specific functions       **********/

  void compute_multipoles_ (Block * block) throw();

  void evaluate_force_ (Block * block) throw();  // compute force from taylor coeffs

  void refresh_acceleration_ (EnzoBlock * enzo_block) throw();

  /********** functions for dual tree walk **********/

public:

  void interact_approx_ (Block * block, MultipoleMsg * msg_b) throw(); // compute Taylor coeffs for two interacting blocks
  void interact_direct_ (Block * block, char * fldbuffer_b, char * prtbuffer_b) throw(); // compute Newtonian force directly for two interacting leaves

  void interact_approx_send(EnzoBlock * enzo_block, Index receiver) throw();

  void interact_direct_send(EnzoBlock * enzo_block, Index receiver) throw();

  void dual_walk_ (EnzoBlock * enzo_block) throw();

  void traverse (EnzoBlock * enzo_block, Index index_b, int type_b);

  /// Identified pair of blocks considered to be "distant"
  void traverse_approx_pair (EnzoBlock * enzo_block,
                            Index index_a, int volume_a,
                            Index index_b, int volume_b);
                              
  /// Identified pair of blocks considered to be "close". (Block
  /// A may be the same as block B but won't be called twice.)
  void traverse_direct_pair (EnzoBlock * enzo_block,
                            Index index_a, int volume_a,
                            Index index_b, int volume_b);

  void pack_dens_ (EnzoBlock * enzo_block, Index index_b, int volume) throw();

  void update_volume (EnzoBlock * enzo_block, Index index, int volume);

protected:

  bool is_far_ (EnzoBlock * enzo_block, Index index_b,
              double * ra, double * rb) const;

  long long * volume_(Block * block) {
    return block->data()->scalar_long_long().value(is_volume_);
  }

  /**********    functions for debugging particles          ************/
  void InitializeParticles (Block * block, int nprtls, double prtls[][7]);


public:
  /**********   convenience functions for tensor arithmetic ***********/

  // shift a quadrupole tensor (old_quadrupole) from center-of-mass old_com to the new center-of-mass new_com
  std::array<double, 6> shift_quadrupole_(double * old_quadrupole, double tot_mass, double * old_com, double * new_com) throw()
  {
    double xdisp = new_com[0] - old_com[0];
    double ydisp = new_com[1] - old_com[1];
    double zdisp = new_com[2] - old_com[2];

    std::array<double, 6> new_quadrupole;

    new_quadrupole[0] = old_quadrupole[0] + tot_mass * xdisp * xdisp;
    new_quadrupole[1] = old_quadrupole[1] + tot_mass * xdisp * ydisp;
    new_quadrupole[2] = old_quadrupole[2] + tot_mass * xdisp * zdisp;
    new_quadrupole[3] = old_quadrupole[3] + tot_mass * ydisp * ydisp;
    new_quadrupole[4] = old_quadrupole[4] + tot_mass * ydisp * zdisp;
    new_quadrupole[5] = old_quadrupole[5] + tot_mass * zdisp * zdisp;
    
    return new_quadrupole;
  }

  // subtract two 3-vectors
  static std::array<double, 3> subtract_(std::array<double, 3> a, std::array<double, 3> b) throw()
  {
    std::array<double, 3> diff;
    
    diff[0] = a[0] - b[0];
    diff[1] = a[1] - b[1];
    diff[2] = a[2] - b[2];

    return diff;
  }

  // add two 3-vectors
  static std::array<double, 3> add_11_(std::array<double, 3> a, std::array<double, 3> b) throw()
  {
    std::array<double, 3> sum;

    sum[0] = a[0] + b[0];
    sum[1] = a[1] + b[1];
    sum[2] = a[2] + b[2];

    return sum;
  }

  // add two rank-2 tensors
  static std::array<double, 6> add_22_(std::array<double, 6> a, std::array<double, 6> b) throw()
  {
    std::array<double, 6> sum;

    sum[0] = a[0] + b[0];
    sum[1] = a[1] + b[1];
    sum[2] = a[2] + b[2];
    sum[3] = a[3] + b[3];
    sum[4] = a[4] + b[4];
    sum[5] = a[5] + b[5];

    return sum;
  }

  // add two rank-3 tensors
  static std::array<double, 10> add_33_(std::array<double, 10> a, std::array<double, 10> b) throw()
  {
    std::array<double, 10> sum;

    sum[0] = a[0] + b[0];
    sum[1] = a[1] + b[1];
    sum[2] = a[2] + b[2];
    sum[3] = a[3] + b[3];
    sum[4] = a[4] + b[4];
    sum[5] = a[5] + b[5];
    sum[6] = a[6] + b[6];
    sum[7] = a[7] + b[7];
    sum[8] = a[8] + b[8];
    sum[9] = a[9] + b[9];

    return sum;
  }

  // compute the outer product of two 3-vectors 
  static std::array<double, 6> outer_11_(std::array<double, 3> a, std::array<double, 3> b) throw()
  {
    std::array<double, 6> prod;

    prod[0] = a[0] * b[0];
    prod[1] = a[0] * b[1];
    prod[2] = a[0] * b[2];
    prod[3] = a[1] * b[1];
    prod[4] = a[1] * b[2];
    prod[5] = a[2] * b[2];

    return prod;
  }

  // compute the outer product of a 3-vector 'a' with a rank-2 tensor 'b'
  static std::array<double, 10> outer_12_(std::array<double, 3> a, std::array<double, 6> b) throw()
  {
    std::array<double, 10> prod;

    prod[0] = a[0] * b[0];  // 000
    prod[1] = a[0] * b[1];  // 001
    prod[2] = a[0] * b[2];  // 002
    prod[3] = a[0] * b[3];  // 011
    prod[4] = a[0] * b[4];  // 012
    prod[5] = a[0] * b[5];  // 022
    prod[6] = a[1] * b[3];  // 111
    prod[7] = a[1] * b[4];  // 112
    prod[8] = a[1] * b[5];  // 122
    prod[9] = a[2] * b[5];  // 222

    return prod;
  }

  // compute the outer product of a 3-vector 'a' with a rank-3 tensor 'b'
  static std::array<double, 15> outer_13_(std::array<double, 3> a, std::array<double, 10> b) throw()
  {
    std::array<double, 15> prod;

    prod[0] = a[0] * b[0];  // 0000
    prod[1] = a[0] * b[1];  // 0001
    prod[2] = a[0] * b[2];  // 0002
    prod[3] = a[0] * b[3];  // 0011
    prod[4] = a[0] * b[4];  // 0012
    prod[5] = a[0] * b[5];  // 0022
    prod[6] = a[0] * b[6];  // 0111
    prod[7] = a[0] * b[7];  // 0112
    prod[8] = a[0] * b[8];  // 0122
    prod[9] = a[0] * b[9];  // 0222
    prod[10] = a[1] * b[6];  // 1111
    prod[11] = a[1] * b[7];  // 1112
    prod[12] = a[1] * b[8];  // 1122
    prod[13] = a[1] * b[9];  // 1222
    prod[14] = a[2] * b[9];  // 2222


    return prod;
  }

  // multiply a 3-vector by a scalar
  static std::array<double, 3> dot_scalar_1_(double a, std::array<double, 3> b) throw()
  {
    std::array<double, 3> prod;

    prod[0] = a * b[0];
    prod[1] = a * b[1];
    prod[2] = a * b[2];

    return prod;
  }

  // multiply a rank-2 tensor by a scalar
  static std::array<double, 6> dot_scalar_2_(double a, std::array<double, 6> b) throw()
  {
    std::array<double, 6> prod;

    prod[0] = a * b[0];
    prod[1] = a * b[1];
    prod[2] = a * b[2];
    prod[3] = a * b[3];
    prod[4] = a * b[4];
    prod[5] = a * b[5];

    return prod;
  }

  // multiply a rank-3 tensor by a scalar
  static std::array<double, 10> dot_scalar_3_(double a, std::array<double, 10> b) throw()
  {
    std::array<double, 10> prod;

    prod[0] = a * b[0];
    prod[1] = a * b[1];
    prod[2] = a * b[2];
    prod[3] = a * b[3];
    prod[4] = a * b[4];
    prod[5] = a * b[5];
    prod[6] = a * b[6];
    prod[7] = a * b[7];
    prod[8] = a * b[8];
    prod[9] = a * b[9];

    return prod;
  }

  // compute the dot product of two 3-vectors
  static double dot_11_(std::array<double, 3> a, std::array<double, 3> b) throw()
  {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }

  // compute the dot product of a 3-vector 'a' with a rank-2 tensor 'b'
  static std::array<double, 3> dot_12_(std::array<double, 3> a, std::array<double, 6> b) throw()
  {
    std::array<double, 3> prod;

    prod[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];   // 00 + 01 + 02
    prod[1] = a[0]*b[1] + a[1]*b[3] + a[2]*b[4];   // 01 + 11 + 12
    prod[2] = a[0]*b[2] + a[1]*b[4] + a[2]*b[5];   // 02 + 12 + 22
    
    return prod;
  }

  // compute the dot product of a 3-vector 'a' with a rank-3 tensor 'b'
  static std::array<double, 6> dot_13_(std::array<double, 3> a, std::array<double, 10> b) throw()
  {
    std::array<double, 6> prod;

    prod[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    prod[1] = a[0]*b[1] + a[1]*b[3] + a[2]*b[4];
    prod[2] = a[0]*b[2] + a[1]*b[4] + a[2]*b[5];
    prod[3] = a[0]*b[3] + a[1]*b[6] + a[2]*b[7];
    prod[4] = a[0]*b[4] + a[1]*b[7] + a[2]*b[8];
    prod[5] = a[0]*b[5] + a[1]*b[8] + a[2]*b[9];

    return prod;
  }

  // compute the dot product of a 3-vector 'a' with a rank-4 tensor 'b'
  static std::array<double, 10> dot_14_(std::array<double, 3> a, std::array<double, 15> b) throw()
  {
    std::array<double, 10> prod;

    prod[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    prod[1] = a[0]*b[1] + a[1]*b[3] + a[2]*b[4];
    prod[2] = a[0]*b[2] + a[1]*b[4] + a[2]*b[5];
    prod[3] = a[0]*b[3] + a[1]*b[6] + a[2]*b[7];
    prod[4] = a[0]*b[4] + a[1]*b[7] + a[2]*b[8];
    prod[5] = a[0]*b[5] + a[1]*b[8] + a[2]*b[9];
    prod[6] = a[0]*b[6] + a[1]*b[10] + a[2]*b[11];
    prod[7] = a[0]*b[7] + a[1]*b[11] + a[2]*b[12];
    prod[8] = a[0]*b[8] + a[1]*b[12] + a[2]*b[13];
    prod[9] = a[0]*b[9] + a[1]*b[13] + a[2]*b[14];

    return prod;
  }

  // compute the dot product of a rank-2 tensor 'a' with a rank-2 tensor 'b'
  // only need this if we want the potential (i.e., only useful for d0)
  static double dot_22_(std::array<double, 6> a, std::array<double, 6> b) throw()
  {
    return (a[0]*b[0] + 2*a[1]*b[1] + 2*a[2]*b[2] 
           + a[3]*b[3] + 2*a[4]*b[4] 
           + a[5]*b[5]);
  }

  // compute the dot product of a rank-2 tensor 'a' with a rank-3 tensor 'b'
  static std::array<double, 3> dot_23_(std::array<double, 6> a, std::array<double, 10> b) throw()
  {
    std::array<double, 3> prod;

    prod[0] = a[0]*b[0] + 2*a[1]*b[1] + 2*a[2]*b[2] + a[3]*b[3] + 2*a[4]*b[4] + a[5]*b[5];
    prod[1] = a[0]*b[1] + 2*a[1]*b[3] + 2*a[2]*b[4] + a[3]*b[6] + 2*a[4]*b[7] + a[5]*b[8];
    prod[2] = a[0]*b[2] + 2*a[1]*b[4] + 2*a[2]*b[5] + a[3]*b[7] + 2*a[4]*b[8] + a[5]*b[9];
 
    return prod;
  }

  // compute the dot product of a rank-2 tensor 'a' with a rank-4 tensor 'b'
  static std::array<double, 6> dot_24_(std::array<double, 6> a, std::array<double, 15> b) throw()
  {
    std::array<double, 6> prod;

    prod[0] = a[0]*b[0] + 2*a[1]*b[1] + 2*a[2]*b[2] + a[3]*b[3] + 2*a[4]*b[4] + a[5]*b[5];
    prod[1] = a[0]*b[1] + 2*a[1]*b[3] + 2*a[2]*b[4] + a[3]*b[6] + 2*a[4]*b[7] + a[5]*b[8];
    prod[2] = a[0]*b[2] + 2*a[1]*b[4] + 2*a[2]*b[5] + a[3]*b[7] + 2*a[4]*b[8] + a[5]*b[9];
    prod[3] = a[0]*b[3] + 2*a[1]*b[6] + 2*a[2]*b[7] + a[3]*b[10] + 2*a[4]*b[11] + a[5]*b[12];
    prod[4] = a[0]*b[4] + 2*a[1]*b[7] + 2*a[2]*b[8] + a[3]*b[11] + 2*a[4]*b[12] + a[5]*b[13];
    prod[5] = a[0]*b[5] + 2*a[1]*b[8] + 2*a[2]*b[9] + a[3]*b[12] + 2*a[4]*b[13] + a[5]*b[14];

    return prod;
  }

  // compute the dot product of a rank-2 tensor 'a' with a rank-5 tensor 'b'
  static std::array<double, 10> dot_25_(std::array<double, 6> a, std::array<double, 21> b) throw()
  {
    std::array<double, 10> prod;

    prod[0] = a[0]*b[0] + 2*a[1]*b[1] + 2*a[2]*b[2] + a[3]*b[3] + 2*a[4]*b[4] + a[5]*b[5];
    prod[1] = a[0]*b[1] + 2*a[1]*b[3] + 2*a[2]*b[4] + a[3]*b[6] + 2*a[4]*b[7] + a[5]*b[8];
    prod[2] = a[0]*b[2] + 2*a[1]*b[4] + 2*a[2]*b[5] + a[3]*b[7] + 2*a[4]*b[8] + a[5]*b[9];
    prod[3] = a[0]*b[3] + 2*a[1]*b[6] + 2*a[2]*b[7] + a[3]*b[10] + 2*a[4]*b[11] + a[5]*b[12];
    prod[4] = a[0]*b[4] + 2*a[1]*b[7] + 2*a[2]*b[8] + a[3]*b[11] + 2*a[4]*b[12] + a[5]*b[13];
    prod[5] = a[0]*b[5] + 2*a[1]*b[8] + 2*a[2]*b[9] + a[3]*b[12] + 2*a[4]*b[13] + a[5]*b[14];
    prod[6] = a[0]*b[6] + 2*a[1]*b[10] + 2*a[2]*b[11] + a[3]*b[15] + 2*a[4]*b[16] + a[5]*b[17];
    prod[7] = a[0]*b[7] + 2*a[1]*b[11] + 2*a[2]*b[12] + a[3]*b[16] + 2*a[4]*b[17] + a[5]*b[18];
    prod[8] = a[0]*b[8] + 2*a[1]*b[12] + 2*a[2]*b[13] + a[3]*b[17] + 2*a[4]*b[18] + a[5]*b[19];
    prod[9] = a[0]*b[9] + 2*a[1]*b[13] + 2*a[2]*b[14] + a[3]*b[18] + 2*a[4]*b[19] + a[5]*b[20];

    return prod;
  }

  // compute the dot product of a rank-3 tensor 'a' with a rank-3 tensor 'b'
  // only need this if we want the potential (i.e., only useful for d0)
  static double dot_33_(std::array<double, 10> a, std::array<double, 10> b) throw()
  {
    return (a[0]*b[0] + 3*a[1]*b[1] + 3*a[2]*b[2]
           + 3*a[3]*b[3] + 6*a[4]*b[4] + 3*a[5]*b[5]
           + a[6]*b[6] + 3*a[7]*b[7] + 3*a[8]*b[8] + a[9]*b[9]);
  }

  // compute the dot product of a rank-3 tensor 'a' with a rank-4 tensor 'b'
  static std::array<double, 3> dot_34_(std::array<double, 10> a, std::array<double, 15> b) throw()
  {
    std::array<double, 3> prod;

    prod[0] = a[0]*b[0] + 3*a[1]*b[1] + 3*a[2]*b[2] + 3*a[3]*b[3] + 6*a[4]*b[4] + 3*a[5]*b[5] + a[6]*b[6] + 3*a[7]*b[7] + 3*a[8]*b[8] + a[9]*b[9];
    prod[1] = a[0]*b[1] + 3*a[1]*b[3] + 3*a[2]*b[4] + 3*a[3]*b[6] + 6*a[4]*b[7] + 3*a[5]*b[8] + a[6]*b[10] + 3*a[7]*b[11] + 3*a[8]*b[12] + a[9]*b[13];
    prod[2] = a[0]*b[2] + 3*a[1]*b[4] + 3*a[2]*b[5] + 3*a[3]*b[7] + 6*a[4]*b[8] + 3*a[5]*b[9] + a[6]*b[11] + 3*a[7]*b[12] + 3*a[8]*b[13] + a[9]*b[14];

    return prod;
  }

  // compute the dot product of a rank-3 tensor 'a' with a rank-5 tensor 'b'
  static std::array<double, 6> dot_35_(std::array<double, 10> a, std::array<double, 21> b) throw()
  {
    std::array<double, 6> prod;

    prod[0] = a[0]*b[0] + 3*a[1]*b[1] + 3*a[2]*b[2] + 3*a[3]*b[3] + 6*a[4]*b[4] + 3*a[5]*b[5] + a[6]*b[6] + 3*a[7]*b[7] + 3*a[8]*b[8] + a[9]*b[9];
    prod[1] = a[0]*b[1] + 3*a[1]*b[3] + 3*a[2]*b[4] + 3*a[3]*b[6] + 6*a[4]*b[7] + 3*a[5]*b[8] + a[6]*b[10] + 3*a[7]*b[11] + 3*a[8]*b[12] + a[9]*b[13];
    prod[2] = a[0]*b[2] + 3*a[1]*b[4] + 3*a[2]*b[5] + 3*a[3]*b[7] + 6*a[4]*b[8] + 3*a[5]*b[9] + a[6]*b[11] + 3*a[7]*b[12] + 3*a[8]*b[13] + a[9]*b[14];
    prod[3] = a[0]*b[3] + 3*a[1]*b[6] + 3*a[2]*b[7] + 3*a[3]*b[10] + 6*a[4]*b[11] + 3*a[5]*b[12] + a[6]*b[15] + 3*a[7]*b[16] + 3*a[8]*b[17] + a[9]*b[18];
    prod[4] = a[0]*b[4] + 3*a[1]*b[7] + 3*a[2]*b[8] + 3*a[3]*b[11] + 6*a[4]*b[12] + 3*a[5]*b[13] + a[6]*b[16] + 3*a[7]*b[17] + 3*a[8]*b[18] + a[9]*b[19];
    prod[5] = a[0]*b[5] + 3*a[1]*b[8] + 3*a[2]*b[9] + 3*a[3]*b[12] + 6*a[4]*b[13] + 3*a[5]*b[14] + a[6]*b[17] + 3*a[7]*b[18] + 3*a[8]*b[19] + a[9]*b[20];

    return prod;
  }

  // compute the dot product of a rank-3 tensor 'a' with a rank-6 tensor 'b'
  static std::array<double, 10> dot_36_(std::array<double, 10> a, std::array<double, 28> b) throw()
  {
    std::array<double, 10> prod;

    prod[0] = a[0]*b[0] + 3*a[1]*b[1] + 3*a[2]*b[2] + 3*a[3]*b[3] + 6*a[4]*b[4] + 3*a[5]*b[5] + a[6]*b[6] + 3*a[7]*b[7] + 3*a[8]*b[8] + a[9]*b[9];
    prod[1] = a[0]*b[1] + 3*a[1]*b[3] + 3*a[2]*b[4] + 3*a[3]*b[6] + 6*a[4]*b[7] + 3*a[5]*b[8] + a[6]*b[10] + 3*a[7]*b[11] + 3*a[8]*b[12] + a[9]*b[13];
    prod[2] = a[0]*b[2] + 3*a[1]*b[4] + 3*a[2]*b[5] + 3*a[3]*b[7] + 6*a[4]*b[8] + 3*a[5]*b[9] + a[6]*b[11] + 3*a[7]*b[12] + 3*a[8]*b[13] + a[9]*b[14];
    prod[3] = a[0]*b[3] + 3*a[1]*b[6] + 3*a[2]*b[7] + 3*a[3]*b[10] + 6*a[4]*b[11] + 3*a[5]*b[12] + a[6]*b[15] + 3*a[7]*b[16] + 3*a[8]*b[17] + a[9]*b[18];
    prod[4] = a[0]*b[4] + 3*a[1]*b[7] + 3*a[2]*b[8] + 3*a[3]*b[11] + 6*a[4]*b[12] + 3*a[5]*b[13] + a[6]*b[16] + 3*a[7]*b[17] + 3*a[8]*b[18] + a[9]*b[19];
    prod[5] = a[0]*b[5] + 3*a[1]*b[8] + 3*a[2]*b[9] + 3*a[3]*b[12] + 6*a[4]*b[13] + 3*a[5]*b[14] + a[6]*b[17] + 3*a[7]*b[18] + 3*a[8]*b[19] + a[9]*b[20];
    prod[6] = a[0]*b[6] + 3*a[1]*b[10] + 3*a[2]*b[11] + 3*a[3]*b[15] + 6*a[4]*b[16] + 3*a[5]*b[17] + a[6]*b[21] + 3*a[7]*b[22] + 3*a[8]*b[23] + a[9]*b[24];
    prod[7] = a[0]*b[7] + 3*a[1]*b[11] + 3*a[2]*b[12] + 3*a[3]*b[16] + 6*a[4]*b[17] + 3*a[5]*b[18] + a[6]*b[22] + 3*a[7]*b[23] + 3*a[8]*b[24] + a[9]*b[25];
    prod[8] = a[0]*b[8] + 3*a[1]*b[12] + 3*a[2]*b[13] + 3*a[3]*b[17] + 6*a[4]*b[18] + 3*a[5]*b[19] + a[6]*b[23] + 3*a[7]*b[24] + 3*a[8]*b[25] + a[9]*b[26];
    prod[9] = a[0]*b[9] + 3*a[1]*b[13] + 3*a[2]*b[14] + 3*a[3]*b[18] + 6*a[4]*b[19] + 3*a[5]*b[20] + a[6]*b[24] + 3*a[7]*b[25] + 3*a[8]*b[26] + a[9]*b[27];

    return prod;
  }


protected: // attributes
  
  // int i_sync_prolong_;
  int i_sync_restrict_;
  int i_msg_prolong_;
  int i_msg_restrict_[8];

  /// Refresh ids
  int ir_exit_;

  /// indices to multipoles for this block
  int i_mass_;
  int i_com_;
  int i_quadrupole_;

  /// indices to Taylor coefficients for this block
  int i_c1_;
  int i_c2_;
  int i_c3_;

  /// Parameter controlling the multipole acceptance criterion, defined
  /// as theta * distance(A,B) > (radius(A) + radius(B))
  double theta_;

  /// Gravitational softening length
  double eps0_;

  /// Cutoff distance for gravitational softening -- unnecessary, since r0 = 2.8*eps0_
  // double r0_;

  // /// Minimum/maximum mesh refinement level (saved for efficiency)  !!! include this!
  // int min_level_;
  // int max_level_;

  /// ScalarData<long long> index for summing block volumes to
  /// test for termination
  int is_volume_;

  /// volume of block in each level, with finest level normalized to 1. Used
  /// to test for termination
  std::vector<int> block_volume_;

  /// Volume of domain in terms of finest blocks
  long long max_volume_;

  /// dimensions of interpolation grid for Ewald summation (periodic BCs)
  int interp_xpoints_;
  int interp_ypoints_;
  int interp_zpoints_;

  /// order (3 or 4) of the Taylor expansion about the nearest grid point
  /// used to interpolate the first derivative of the Ewald correction
  int ewald_order_;

  /// object for storing interpolation grid needed for Ewald summation (periodic BCs)
  EnzoEwald * ewald_;

  /// Maximum timestep
  double dt_max_;

  /// Tree-based leaf-leaf interaction (see EnzoMethodMultipoleSubtree.hpp).
  /// Each leaf Block is divided into chunks of leaf_chunk_^3 cells that form a
  /// local octree, and a dual tree walk with opening angle leaf_theta_
  /// replaces the brute-force sums. leaf_chunk_ = 0 (default) turns it off.
  int leaf_chunk_;
  double leaf_theta_;

};

#endif /* ENZO_METHOD_MULTIPOLE_HPP */

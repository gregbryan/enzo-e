// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MultipoleMsg.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu)
/// @date     2022-7-20
/// @brief    [\ref Charm] Declaration of the MultipoleMsg Charm++ Message

#ifndef CHARM_MULTIPOLE_MSG_HPP
#define CHARM_MULTIPOLE_MSG_HPP

class MultipoleMsg : public CMessage_MultipoleMsg {

// public:
//   MultipoleMsg();
//   ~MultipoleMsg();

//   MultipoleMsg(const MultipoleMsg&) {
//     ERROR ("MultipoleMsg::MultipoleMsg",
// 	     "copy constructor not supported");
//   }
  
//   MultipoleMsg & operator= (const MultipoleMsg & data_msg) throw() {
//      ERROR ("MultipoleMsg::MultipoleMsg",
// 	     "copy assignment not supported");
//   } 

// public: // static methods

//   /// Pack data to serialize
//   static void * pack (MultipoleMsg*);

//   /// Unpack data to de-serialize
//   static MultipoleMsg * unpack(void *);

//   void set_data_msg  (DataMsg * data_msg);

//   void update (Data * data);
  

public: // attributes


  int child_index() const { return ic3[0] + 2*(ic3[1] + 2*(ic3[2])); }
  
  /// total mass
  double mass;

  /// center of mass
  double com[3];

  /// quadrupole moment
  double quadrupole[6];

  /// first Taylor coefficient
  double c1[3];

  /// second Taylor coefficient
  double c2[6];

  /// third Taylor coefficient
  double c3[10];

  /// Child indices
  int ic3[3];

  /// Volume to be credited to the receiver once the message was applied
  /// (interact_approx only)
  int volume;


  /// lower left corner of Block
  // double lo[3];

  // double up[3];

  // int size[3];

  /// density field
  // DataMsg * data_msg_;

  // void * buffer_;
  
};

#endif /* CHARM_MULTIPOLE_MSG_HPP */

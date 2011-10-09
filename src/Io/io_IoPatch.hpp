// See LICENSE_CELLO file for license and copyright information

#ifndef IO_IO_PATCH_HPP
#define IO_IO_PATCH_HPP

/// @file     io_IoPatch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoPatch class
///

class Patch;

class IoPatch : public Io {

  /// @class    IoPatch
  /// @ingroup  Io
  /// @brief    [\ref Io] 

public: // interface

  /// Constructor
  IoPatch(const Patch * patch) throw();

  /// Return the ith metadata item associated with the Patch object
  void meta_value 
  (int index, 
   void ** buffer, const char ** name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  /// Return the ith data item associated with the Patch object
  void data_value 
  (int index, 
   void ** buffer, const char ** name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  
private: // functions

  const Patch * patch_;

private: // attributes


};

#endif /* IO_IO_PATCH_HPP */

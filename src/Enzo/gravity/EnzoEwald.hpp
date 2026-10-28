// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoEwald.hpp
/// @author   Ryan Golant (ryan.golant@columbia.edu) 
/// @date     Thursday September 28, 2023
/// @brief    [\ref Enzo] Declaration of EnzoEwald
///           Compute Ewald sums for periodic boundary conditions

#ifndef ENZO_EWALD_HPP
#define ENZO_EWALD_HPP

// #include "view.hpp"

class EnzoEwald {

  /// @class    EnzoEwald
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Compute Ewald sums for periodic boundary conditions

public: 

  /// Create a new EnzoEwald object. "order" is the order of the Taylor
  /// expansion about the nearest grid point used by interp_d1(): 3 or 4
  /// (the fourth-order term needs the 5th derivatives, which are tabulated)
  EnzoEwald (int interp_xpoints, int interp_ypoints, int interp_zpoints,
             int order = 3);

  EnzoEwald()
    : interp_xpoints_(32),  // number of interpolation points in the x-direction
      interp_ypoints_(32),  // number of interpolation points in the y-direction
      interp_zpoints_(32),  // number of interpolation points in the z-direction
      order_(3)             // order of the Taylor expansion in interp_d1()
  { }

  void pup(PUP::er &p) {
    p | interp_xpoints_;
    p | interp_ypoints_;
    p | interp_zpoints_;
    p | order_;
  }

  void init_interpolate_() throw();

  void find_nearest_interp_point(double x, double y, double z, 
                                 double * interp_x, double * interp_y, double * interp_z, int * i) throw()
  {
    double lox, loy, loz; 
    double hix, hiy, hiz;
    cello::hierarchy()->lower(&lox, &loy, &loz);
    cello::hierarchy()->upper(&hix, &hiy, &hiz);

    double Lx = hix - lox;
    double Ly = hiy - loy;
    double Lz = hiz - loz;

    double dx = Lx / (interp_xpoints_ - 1);
    double dy = Ly / (interp_ypoints_ - 1);
    double dz = Lz / (interp_zpoints_ - 1);

    int ix = round((x+Lx/2.)/dx);
    int iy = round((y+Ly/2.)/dy);
    int iz = round((z+Lz/2.)/dz);

    *i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);

    *interp_x = -Lx/2. + ix*dx;
    *interp_y = -Ly/2. + iy*dy;
    *interp_z = -Lz/2. + iz*dz;

    // CkPrintf("%f, %f, %f, %d, %f, %f, %f\n", x, y, z, *i, *interp_x, *interp_y, *interp_z);

    return;
  }

  // double interp_d0(double x, double y, double z) throw();  // this function is not necessary for computing acceleration
  std::array<double, 3> interp_d1(double x, double y, double z) throw();
  std::array<double, 6> interp_d2(double x, double y, double z) throw();
  std::array<double, 10> interp_d3(double x, double y, double z) throw();

  // double d0(double x, double y, double z) throw();  // this function is not necessary for computing acceleration 
  // std::array<double, 3> d1(double x, double y, double z) throw();
  // std::array<double, 6> d2(double x, double y, double z) throw();
  // std::array<double, 10> d3(double x, double y, double z) throw();
  // std::array<double, 15> d4(double x, double y, double z) throw();
  // std::array<double, 21> d5(double x, double y, double z) throw();
  // std::array<double, 28> d6(double x, double y, double z) throw();

  CelloView<double, 1> d1(double x, double y, double z) throw();
  CelloView<double, 1> d2(double x, double y, double z) throw();
  CelloView<double, 1> d2_gadget(double x, double y, double z) throw();
  CelloView<double, 1> d3(double x, double y, double z) throw();
  CelloView<double, 1> d4(double x, double y, double z) throw();
  CelloView<double, 1> d5(double x, double y, double z) throw();
  CelloView<double, 1> d6(double x, double y, double z) throw();

protected: //attributes

  // dimensions of interpolation grid (Nx x Ny x Nz)
  int interp_xpoints_;
  int interp_ypoints_;
  int interp_zpoints_;

  // order of the Taylor expansion about the nearest grid point in interp_d1()
  int order_;

  // int tot_points_ = interp_xpoints_ * interp_ypoints_ * interp_zpoints_;
  //const int tot_points_ = 10

  /// derivative interpolation arrays
  // std::array<double, tot_points_> d0_array_; 
  // std::array<std::array<double, 3>, tot_points_> d1_array_; 
  // std::array<std::array<double, 6>, tot_points_> d2_array_; 
  // std::array<std::array<double, 10>, tot_points_> d3_array_; 
  // std::array<std::array<double, 15>, tot_points_> d4_array_; 
  // std::array<std::array<double, 21>, tot_points_> d5_array_; 
  // std::array<std::array<double, 28>, tot_points_> d6_array_; 

  CelloView<double, 2> d1_array_; 
  CelloView<double, 2> d2_array_; 
  CelloView<double, 2> d3_array_; 
  CelloView<double, 2> d4_array_; 
  CelloView<double, 2> d5_array_; 
  CelloView<double, 2> d6_array_; 


  
  

};

#endif /* ENZO_EWALD_HPP */

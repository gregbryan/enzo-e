// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoEwald.cpp
/// @author   Ryan Golant (ryan.golant@columbia.edu)
/// @date     September 28, 2023
/// @brief    Compute Ewald sums for periodic boundary conditions


#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/gravity/gravity.hpp"
// #include "view.hpp"


EnzoEwald::EnzoEwald (int interp_xpoints, int interp_ypoints, int interp_zpoints, int order)
  : interp_xpoints_(interp_xpoints),  // number of interpolation points in the x-direction
    interp_ypoints_(interp_ypoints),  // number of interpolation points in the y-direction
    interp_zpoints_(interp_zpoints),  // number of interpolation points in the z-direction
    order_(order)                     // order of the Taylor expansion in interp_d1()
    
{ 
  // EnzoEwald constructor is called in *compute* of EnzoMethodMultipole

  int tot_points_ = interp_xpoints_ * interp_ypoints_ * interp_zpoints_;

  // d0_array_ = {};    // Nx x Ny x Nz x 1 (on down-sampled grid of dimension Nx x Ny x Nz)
  d1_array_ = CelloView<double, 2> (tot_points_, 3);  // Nx x Ny x Nz x 3
  d2_array_ = CelloView<double, 2> (tot_points_, 6);  // Nx x Ny x Nz x 6
  d3_array_ = CelloView<double, 2> (tot_points_, 10); // Nx x Ny x Nz x 10
  d4_array_ = CelloView<double, 2> (tot_points_, 15); // Nx x Ny x Nz x 15
  d5_array_ = CelloView<double, 2> (tot_points_, 21); // Nx x Ny x Nz x 21
  d6_array_ = CelloView<double, 2> (tot_points_, 28); // Nx x Ny x Nz x 28
  
  CkPrintf("before init interpolate\n");

  init_interpolate_();

  CkPrintf("after init interpolate\n");

}

// set up interpolation arrays on the primary domain
void EnzoEwald::init_interpolate_() throw()
{

  Hierarchy * hierarchy = enzo::simulation()->hierarchy();

  double lox, loy, loz; 
  double hix, hiy, hiz;
  
  hierarchy->lower(&lox, &loy, &loz);
  hierarchy->upper(&hix, &hiy, &hiz);

  double Lx = hix - lox;
  double Ly = hiy - loy;
  double Lz = hiz - loz;
  
  double dx = (interp_xpoints_ > 1) ? Lx / (interp_xpoints_ - 1) : 0;
  double dy = (interp_ypoints_ > 1) ? Ly / (interp_ypoints_ - 1) : 0;
  double dz = (interp_zpoints_ > 1) ? Lz / (interp_zpoints_ - 1) : 0;

  int midx = floor((interp_xpoints_ - 1)/2);
  int midy = floor((interp_ypoints_ - 1)/2);
  int midz = floor((interp_zpoints_ - 1)/2);

  
  for (int iz = 0; iz <= midz; iz++) {
    for (int iy = 0; iy <= midy; iy++) {
      for (int ix = 0; ix <= midx; ix++) {
        
        int i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);
        double x = -Lx/2. + ix*dx;  
        double y = -Ly/2. + iy*dy;
        double z = -Lz/2. + iz*dz;

        CelloView<double, 1> d1_slice_ = d1(x, y, z);
        CelloView<double, 1> d2_slice_ = d2(x, y, z);
        CelloView<double, 1> d3_slice_ = d3(x, y, z);
        CelloView<double, 1> d4_slice_ = d4(x, y, z);
        CelloView<double, 1> d5_slice_ = d5(x, y, z);
        CelloView<double, 1> d6_slice_ = d6(x, y, z);

        for (int j = 0; j < 3; j++) {
          d1_array_(i, j) = d1_slice_(j);
        }

        for (int j = 0; j < 6; j++) {
          d2_array_(i, j) = d2_slice_(j);
        }

        for (int j = 0; j < 10; j++) {
          d3_array_(i, j) = d3_slice_(j);
        }

        for (int j = 0; j < 15; j++) {
          d4_array_(i, j) = d4_slice_(j);
        }

        for (int j = 0; j < 21; j++) {
          d5_array_(i, j) = d5_slice_(j);
        }

        for (int j = 0; j < 28; j++) {
          d6_array_(i, j) = d6_slice_(j);
        }


      }
    }
  }

  for (int iz = 0; iz <= midz; iz++) {
    for (int iy = 0; iy <= midy; iy++) {
      for (int ix = midx+1; ix < interp_xpoints_; ix++) {
        
        int i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);
        int ix_flip = interp_xpoints_ - 1 - ix;
        int i_flip = ix_flip + interp_xpoints_ * (iy + iz * interp_ypoints_);
        
        d1_array_(i, 0) = -d1_array_(i_flip, 0); // x
        d1_array_(i, 1) = d1_array_(i_flip, 1);  // y
        d1_array_(i, 2) = d1_array_(i_flip, 2);  // z
        
        d2_array_(i, 0) = d2_array_(i_flip, 0);  // xx
        d2_array_(i, 1) = -d2_array_(i_flip, 1); // xy
        d2_array_(i, 2) = -d2_array_(i_flip, 2); // xz
        d2_array_(i, 3) = d2_array_(i_flip, 3); // yy
        d2_array_(i, 4) = d2_array_(i_flip, 4); // yz
        d2_array_(i, 5) = d2_array_(i_flip, 5); // zz
        
        d3_array_(i, 0) = -d3_array_(i_flip, 0); // xxx
        d3_array_(i, 1) = d3_array_(i_flip, 1); // xxy
        d3_array_(i, 2) = d3_array_(i_flip, 2); // xxz
        d3_array_(i, 3) = -d3_array_(i_flip, 3); // xyy
        d3_array_(i, 4) = -d3_array_(i_flip, 4); // xyz
        d3_array_(i, 5) = -d3_array_(i_flip, 5); // xzz
        d3_array_(i, 6) = d3_array_(i_flip, 6); // yyy
        d3_array_(i, 7) = d3_array_(i_flip, 7); // yyz
        d3_array_(i, 8) = d3_array_(i_flip, 8); // yzz
        d3_array_(i, 9) = d3_array_(i_flip, 9); // zzz
        
        d4_array_(i, 0) = d4_array_(i_flip, 0); // xxxx
        d4_array_(i, 1) = -d4_array_(i_flip, 1); // xxxy
        d4_array_(i, 2) = -d4_array_(i_flip, 2); // xxxz
        d4_array_(i, 3) = d4_array_(i_flip, 3); // xxyy
        d4_array_(i, 4) = d4_array_(i_flip, 4); // xxyz
        d4_array_(i, 5) = d4_array_(i_flip, 5); // xxzz
        d4_array_(i, 6) = -d4_array_(i_flip, 6); // xyyy
        d4_array_(i, 7) = -d4_array_(i_flip, 7); // xyyz
        d4_array_(i, 8) = -d4_array_(i_flip, 8); // xyzz
        d4_array_(i, 9) = -d4_array_(i_flip, 9); // xzzz
        d4_array_(i, 10) = d4_array_(i_flip, 10); // yyyy
        d4_array_(i, 11) = d4_array_(i_flip, 11); // yyyz
        d4_array_(i, 12) = d4_array_(i_flip, 12); // yyzz
        d4_array_(i, 13) = d4_array_(i_flip, 13); // yzzz
        d4_array_(i, 14) = d4_array_(i_flip, 14); // zzzz
        
        d5_array_(i, 0) = -d5_array_(i_flip, 0); // xxxxx
        d5_array_(i, 1) = d5_array_(i_flip, 1); // xxxxy
        d5_array_(i, 2) = d5_array_(i_flip, 2); // xxxxz
        d5_array_(i, 3) = -d5_array_(i_flip, 3); // xxxyy  
        d5_array_(i, 4) = -d5_array_(i_flip, 4); // xxxyz
        d5_array_(i, 5) = -d5_array_(i_flip, 5); // xxxzz
        d5_array_(i, 6) = d5_array_(i_flip, 6); // xxyyy
        d5_array_(i, 7) = d5_array_(i_flip, 7); // xxyyz
        d5_array_(i, 8) = d5_array_(i_flip, 8); // xxyzz
        d5_array_(i, 9) = d5_array_(i_flip, 9); // xxzzz
        d5_array_(i, 10) = -d5_array_(i_flip, 10); // xyyyy
        d5_array_(i, 11) = -d5_array_(i_flip, 11); // xyyyz
        d5_array_(i, 12) = -d5_array_(i_flip, 12); // xyyzz
        d5_array_(i, 13) = -d5_array_(i_flip, 13); // xyzzz
        d5_array_(i, 14) = -d5_array_(i_flip, 14); // xzzzz
        d5_array_(i, 15) = d5_array_(i_flip, 15); // yyyyy
        d5_array_(i, 16) = d5_array_(i_flip, 16); // yyyyz
        d5_array_(i, 17) = d5_array_(i_flip, 17); // yyyzz
        d5_array_(i, 18) = d5_array_(i_flip, 18); // yyzzz
        d5_array_(i, 19) = d5_array_(i_flip, 19); // yzzzz
        d5_array_(i, 20) = d5_array_(i_flip, 20); // zzzzz
        
        d6_array_(i, 0) = d6_array_(i_flip, 0); // xxxxxx
        d6_array_(i, 1) = -d6_array_(i_flip, 1); // xxxxxy
        d6_array_(i, 2) = -d6_array_(i_flip, 2); // xxxxxz
        d6_array_(i, 3) = d6_array_(i_flip, 3); // xxxxyy
        d6_array_(i, 4) = d6_array_(i_flip, 4); // xxxxyz
        d6_array_(i, 5) = d6_array_(i_flip, 5); // xxxxzz
        d6_array_(i, 6) = -d6_array_(i_flip, 6); // xxxyyy
        d6_array_(i, 7) = -d6_array_(i_flip, 7); // xxxyyz
        d6_array_(i, 8) = -d6_array_(i_flip, 8); // xxxyzz
        d6_array_(i, 9) = -d6_array_(i_flip, 9); // xxxzzz
        d6_array_(i, 10) = d6_array_(i_flip, 10); // xxyyyy
        d6_array_(i, 11) = d6_array_(i_flip, 11); // xxyyyz
        d6_array_(i, 12) = d6_array_(i_flip, 12); // xxyyzz
        d6_array_(i, 13) = d6_array_(i_flip, 13); // xxyzzz
        d6_array_(i, 14) = d6_array_(i_flip, 14); // xxzzzz
        d6_array_(i, 15) = -d6_array_(i_flip, 15); // xyyyyy
        d6_array_(i, 16) = -d6_array_(i_flip, 16); // xyyyyz
        d6_array_(i, 17) = -d6_array_(i_flip, 17); // xyyyzz
        d6_array_(i, 18) = -d6_array_(i_flip, 18); // xyyzzz
        d6_array_(i, 19) = -d6_array_(i_flip, 19); // xyzzzz
        d6_array_(i, 20) = -d6_array_(i_flip, 20); // xzzzzz
        d6_array_(i, 21) = d6_array_(i_flip, 21); // yyyyyy
        d6_array_(i, 22) = d6_array_(i_flip, 22); // yyyyyz
        d6_array_(i, 23) = d6_array_(i_flip, 23); // yyyyzz
        d6_array_(i, 24) = d6_array_(i_flip, 24); // yyyzzz
        d6_array_(i, 25) = d6_array_(i_flip, 25); // yyzzzz
        d6_array_(i, 26) = d6_array_(i_flip, 26); // yzzzzz
        d6_array_(i, 27) = d6_array_(i_flip, 27); // zzzzzz

      }
    }
  }

  for (int iz = 0; iz <= midz; iz++) {
    for (int iy = midy+1; iy < interp_ypoints_; iy++) {
      for (int ix = 0; ix <= midx; ix++) {
        
        int i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);
        int iy_flip = interp_ypoints_ - 1 - iy;
        int i_flip = ix + interp_xpoints_ * (iy_flip + iz * interp_ypoints_);

        d1_array_(i, 0) = d1_array_(i_flip, 0); // x
        d1_array_(i, 1) = -d1_array_(i_flip, 1); // y
        d1_array_(i, 2) = d1_array_(i_flip, 2); // z
        
        d2_array_(i, 0) = d2_array_(i_flip, 0); // xx
        d2_array_(i, 1) = -d2_array_(i_flip, 1); // xy
        d2_array_(i, 2) = d2_array_(i_flip, 2); // xz
        d2_array_(i, 3) = d2_array_(i_flip, 3); // yy
        d2_array_(i, 4) = -d2_array_(i_flip, 4); // yz
        d2_array_(i, 5) = d2_array_(i_flip, 5); // zz
        
        d3_array_(i, 0) = d3_array_(i_flip, 0); // xxx
        d3_array_(i, 1) = -d3_array_(i_flip, 1); // xxy
        d3_array_(i, 2) = d3_array_(i_flip, 2); // xxz
        d3_array_(i, 3) = d3_array_(i_flip, 3); // xyy
        d3_array_(i, 4) = -d3_array_(i_flip, 4); // xyz
        d3_array_(i, 5) = d3_array_(i_flip, 5); // xzz
        d3_array_(i, 6) = -d3_array_(i_flip, 6); // yyy
        d3_array_(i, 7) = d3_array_(i_flip, 7); // yyz
        d3_array_(i, 8) = -d3_array_(i_flip, 8); // yzz
        d3_array_(i, 9) = d3_array_(i_flip, 9); // zzz
        
        d4_array_(i, 0) = d4_array_(i_flip, 0); // xxxx
        d4_array_(i, 1) = -d4_array_(i_flip, 1); // xxxy
        d4_array_(i, 2) = d4_array_(i_flip, 2); // xxxz
        d4_array_(i, 3) = d4_array_(i_flip, 3); // xxyy
        d4_array_(i, 4) = -d4_array_(i_flip, 4); // xxyz
        d4_array_(i, 5) = d4_array_(i_flip, 5); // xxzz
        d4_array_(i, 6) = -d4_array_(i_flip, 6); // xyyy
        d4_array_(i, 7) = d4_array_(i_flip, 7); // xyyz
        d4_array_(i, 8) = -d4_array_(i_flip, 8); // xyzz
        d4_array_(i, 9) = d4_array_(i_flip, 9); // xzzz
        d4_array_(i, 10) = d4_array_(i_flip, 10); // yyyy
        d4_array_(i, 11) = -d4_array_(i_flip, 11); // yyyz
        d4_array_(i, 12) = d4_array_(i_flip, 12); // yyzz
        d4_array_(i, 13) = -d4_array_(i_flip, 13); // yzzz
        d4_array_(i, 14) = d4_array_(i_flip, 14); // zzzz
        
        d5_array_(i, 0) = d5_array_(i_flip, 0); // xxxxx
        d5_array_(i, 1) = -d5_array_(i_flip, 1); // xxxxy
        d5_array_(i, 2) = d5_array_(i_flip, 2); // xxxxz
        d5_array_(i, 3) = d5_array_(i_flip, 3); // xxxyy
        d5_array_(i, 4) = -d5_array_(i_flip, 4); // xxxyz
        d5_array_(i, 5) = d5_array_(i_flip, 5); // xxxzz
        d5_array_(i, 6) = -d5_array_(i_flip, 6); // xxyyy
        d5_array_(i, 7) = d5_array_(i_flip, 7); // xxyyz
        d5_array_(i, 8) = -d5_array_(i_flip, 8); // xxyzz
        d5_array_(i, 9) = d5_array_(i_flip, 9); // xxzzz
        d5_array_(i, 10) = d5_array_(i_flip, 10); // xyyyy
        d5_array_(i, 11) = -d5_array_(i_flip, 11); // xyyyz
        d5_array_(i, 12) = d5_array_(i_flip, 12); // xyyzz
        d5_array_(i, 13) = -d5_array_(i_flip, 13); // xyzzz
        d5_array_(i, 14) = d5_array_(i_flip, 14); // xzzzz
        d5_array_(i, 15) = -d5_array_(i_flip, 15); // yyyyy
        d5_array_(i, 16) = d5_array_(i_flip, 16); // yyyyz
        d5_array_(i, 17) = -d5_array_(i_flip, 17); // yyyzz
        d5_array_(i, 18) = d5_array_(i_flip, 18); // yyzzz
        d5_array_(i, 19) = -d5_array_(i_flip, 19); // yzzzz
        d5_array_(i, 20) = d5_array_(i_flip, 20); // zzzzz
        
        d6_array_(i, 0) = d6_array_(i_flip, 0); // xxxxxx
        d6_array_(i, 1) = -d6_array_(i_flip, 1); // xxxxxy
        d6_array_(i, 2) = d6_array_(i_flip, 2); // xxxxxz
        d6_array_(i, 3) = d6_array_(i_flip, 3); // xxxxyy
        d6_array_(i, 4) = -d6_array_(i_flip, 4); // xxxxyz
        d6_array_(i, 5) = d6_array_(i_flip, 5); // xxxxzz
        d6_array_(i, 6) = -d6_array_(i_flip, 6); // xxxyyy
        d6_array_(i, 7) = d6_array_(i_flip, 7); // xxxyyz
        d6_array_(i, 8) = -d6_array_(i_flip, 8); // xxxyzz
        d6_array_(i, 9) = d6_array_(i_flip, 9); // xxxzzz
        d6_array_(i, 10) = d6_array_(i_flip, 10); // xxyyyy
        d6_array_(i, 11) = -d6_array_(i_flip, 11); // xxyyyz
        d6_array_(i, 12) = d6_array_(i_flip, 12); // xxyyzz
        d6_array_(i, 13) = -d6_array_(i_flip, 13); // xxyzzz
        d6_array_(i, 14) = d6_array_(i_flip, 14); // xxzzzz
        d6_array_(i, 15) = -d6_array_(i_flip, 15); // xyyyyy
        d6_array_(i, 16) = d6_array_(i_flip, 16); // xyyyyz
        d6_array_(i, 17) = -d6_array_(i_flip, 17); // xyyyzz
        d6_array_(i, 18) = d6_array_(i_flip, 18); // xyyzzz
        d6_array_(i, 19) = -d6_array_(i_flip, 19); // xyzzzz
        d6_array_(i, 20) = d6_array_(i_flip, 20); // xzzzzz
        d6_array_(i, 21) = d6_array_(i_flip, 21); // yyyyyy
        d6_array_(i, 22) = -d6_array_(i_flip, 22); // yyyyyz
        d6_array_(i, 23) = d6_array_(i_flip, 23); // yyyyzz
        d6_array_(i, 24) = -d6_array_(i_flip, 24); // yyyzzz
        d6_array_(i, 25) = d6_array_(i_flip, 25); // yyzzzz
        d6_array_(i, 26) = -d6_array_(i_flip, 26); // yzzzzz
        d6_array_(i, 27) = d6_array_(i_flip, 27); // zzzzzz
      }
    }
  }

  for (int iz = midz+1; iz < interp_zpoints_; iz++) {
    for (int iy = 0; iy <= midy; iy++) {
      for (int ix = 0; ix <= midx; ix++) {
          
          int i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);
          int iz_flip = interp_zpoints_ - 1 - iz;
          int i_flip = ix + interp_xpoints_ * (iy + iz_flip * interp_ypoints_);

          d1_array_(i, 0) = d1_array_(i_flip, 0); // x
          d1_array_(i, 1) = d1_array_(i_flip, 1); // y
          d1_array_(i, 2) = -d1_array_(i_flip, 2); // z
          
          d2_array_(i, 0) = d2_array_(i_flip, 0); // xx
          d2_array_(i, 1) = d2_array_(i_flip, 1); // xy
          d2_array_(i, 2) = -d2_array_(i_flip, 2); // xz
          d2_array_(i, 3) = d2_array_(i_flip, 3); // yy
          d2_array_(i, 4) = -d2_array_(i_flip, 4); // yz
          d2_array_(i, 5) = d2_array_(i_flip, 5); // zz
          
          d3_array_(i, 0) = d3_array_(i_flip, 0); // xxx
          d3_array_(i, 1) = d3_array_(i_flip, 1); // xxy
          d3_array_(i, 2) = -d3_array_(i_flip, 2); // xxz
          d3_array_(i, 3) = d3_array_(i_flip, 3); // xyy
          d3_array_(i, 4) = -d3_array_(i_flip, 4); // xyz
          d3_array_(i, 5) = d3_array_(i_flip, 5); // xzz
          d3_array_(i, 6) = d3_array_(i_flip, 6); // yyy
          d3_array_(i, 7) = -d3_array_(i_flip, 7); // yyz
          d3_array_(i, 8) = d3_array_(i_flip, 8); // yzz
          d3_array_(i, 9) = -d3_array_(i_flip, 9); // zzz
          
          d4_array_(i, 0) = d4_array_(i_flip, 0); // xxxx
          d4_array_(i, 1) = d4_array_(i_flip, 1); // xxxy
          d4_array_(i, 2) = -d4_array_(i_flip, 2); // xxxz
          d4_array_(i, 3) = d4_array_(i_flip, 3); // xxyy
          d4_array_(i, 4) = -d4_array_(i_flip, 4); // xxyz
          d4_array_(i, 5) = d4_array_(i_flip, 5); // xxzz
          d4_array_(i, 6) = d4_array_(i_flip, 6); // xyyy
          d4_array_(i, 7) = -d4_array_(i_flip, 7); // xyyz
          d4_array_(i, 8) = d4_array_(i_flip, 8); // xyzz
          d4_array_(i, 9) = -d4_array_(i_flip, 9); // xzzz
          d4_array_(i, 10) = d4_array_(i_flip, 10); // yyyy
          d4_array_(i, 11) = -d4_array_(i_flip, 11); // yyyz
          d4_array_(i, 12) = d4_array_(i_flip, 12); // yyzz
          d4_array_(i, 13) = -d4_array_(i_flip, 13); // yzzz
          d4_array_(i, 14) = d4_array_(i_flip, 14); // zzzz
          
          d5_array_(i, 0) = d5_array_(i_flip, 0); // xxxxx
          d5_array_(i, 1) = d5_array_(i_flip, 1); // xxxxy
          d5_array_(i, 2) = -d5_array_(i_flip, 2); // xxxxz
          d5_array_(i, 3) = d5_array_(i_flip, 3); // xxxyy
          d5_array_(i, 4) = -d5_array_(i_flip, 4); // xxxyz
          d5_array_(i, 5) = d5_array_(i_flip, 5); // xxxzz
          d5_array_(i, 6) = d5_array_(i_flip, 6); // xxyyy
          d5_array_(i, 7) = -d5_array_(i_flip, 7); // xxyyz
          d5_array_(i, 8) = d5_array_(i_flip, 8); // xxyzz
          d5_array_(i, 9) = -d5_array_(i_flip, 9); // xxzzz
          d5_array_(i, 10) = d5_array_(i_flip, 10); // xyyyy
          d5_array_(i, 11) = -d5_array_(i_flip, 11); // xyyyz
          d5_array_(i, 12) = d5_array_(i_flip, 12); // xyyzz
          d5_array_(i, 13) = -d5_array_(i_flip, 13); // xyzzz
          d5_array_(i, 14) = d5_array_(i_flip, 14); // xzzzz
          d5_array_(i, 15) = d5_array_(i_flip, 15); // yyyyy
          d5_array_(i, 16) = -d5_array_(i_flip, 16); // yyyyz
          d5_array_(i, 17) = d5_array_(i_flip, 17); // yyyzz
          d5_array_(i, 18) = -d5_array_(i_flip, 18); // yyzzz
          d5_array_(i, 19) = d5_array_(i_flip, 19); // yzzzz
          d5_array_(i, 20) = -d5_array_(i_flip, 20); // zzzzz
          
          d6_array_(i, 0) = d6_array_(i_flip, 0); // xxxxxx
          d6_array_(i, 1) = d6_array_(i_flip, 1); // xxxxxy
          d6_array_(i, 2) = -d6_array_(i_flip, 2); // xxxxxz
          d6_array_(i, 3) = d6_array_(i_flip, 3); // xxxxyy
          d6_array_(i, 4) = -d6_array_(i_flip, 4); // xxxxyz
          d6_array_(i, 5) = d6_array_(i_flip, 5); // xxxxzz
          d6_array_(i, 6) = d6_array_(i_flip, 6); // xxxyyy
          d6_array_(i, 7) = -d6_array_(i_flip, 7); // xxxyyz
          d6_array_(i, 8) = d6_array_(i_flip, 8); // xxxyzz
          d6_array_(i, 9) = -d6_array_(i_flip, 9); // xxxzzz
          d6_array_(i, 10) = d6_array_(i_flip, 10); // xxyyyy
          d6_array_(i, 11) = -d6_array_(i_flip, 11); // xxyyyz
          d6_array_(i, 12) = d6_array_(i_flip, 12); // xxyyzz
          d6_array_(i, 13) = -d6_array_(i_flip, 13); // xxyzzz
          d6_array_(i, 14) = d6_array_(i_flip, 14); // xxzzzz
          d6_array_(i, 15) = d6_array_(i_flip, 15); // xyyyyy
          d6_array_(i, 16) = -d6_array_(i_flip, 16); // xyyyyz
          d6_array_(i, 17) = d6_array_(i_flip, 17); // xyyyzz
          d6_array_(i, 18) = -d6_array_(i_flip, 18); // xyyzzz
          d6_array_(i, 19) = d6_array_(i_flip, 19); // xyzzzz
          d6_array_(i, 20) = -d6_array_(i_flip, 20); // xzzzzz
          d6_array_(i, 21) = d6_array_(i_flip, 21); // yyyyyy
          d6_array_(i, 22) = -d6_array_(i_flip, 22); // yyyyyz
          d6_array_(i, 23) = d6_array_(i_flip, 23); // yyyyzz
          d6_array_(i, 24) = -d6_array_(i_flip, 24); // yyyzzz
          d6_array_(i, 25) = d6_array_(i_flip, 25); // yyzzzz
          d6_array_(i, 26) = -d6_array_(i_flip, 26); // yzzzzz
          d6_array_(i, 27) = d6_array_(i_flip, 27); // zzzzzz
      }
    }
  }
            
  for (int iz = 0; iz <= midz; iz++) {
    for (int iy = midy+1; iy < interp_ypoints_; iy++) {
      for (int ix = midx+1; ix < interp_xpoints_; ix++) {
          
          int i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);
          int ix_flip = interp_xpoints_ - 1 - ix;
          int iy_flip = interp_ypoints_ - 1 - iy;
          int i_flip = ix_flip + interp_xpoints_ * (iy_flip + iz * interp_ypoints_);

          d1_array_(i, 0) = -d1_array_(i_flip, 0); // x
          d1_array_(i, 1) = -d1_array_(i_flip, 1); // y
          d1_array_(i, 2) = d1_array_(i_flip, 2); // z
          
          d2_array_(i, 0) = d2_array_(i_flip, 0); // xx
          d2_array_(i, 1) = d2_array_(i_flip, 1); // xy
          d2_array_(i, 2) = -d2_array_(i_flip, 2); // xz
          d2_array_(i, 3) = d2_array_(i_flip, 3); // yy
          d2_array_(i, 4) = -d2_array_(i_flip, 4); // yz
          d2_array_(i, 5) = d2_array_(i_flip, 5); // zz
          
          d3_array_(i, 0) = -d3_array_(i_flip, 0); // xxx
          d3_array_(i, 1) = -d3_array_(i_flip, 1); // xxy
          d3_array_(i, 2) = d3_array_(i_flip, 2); // xxz
          d3_array_(i, 3) = -d3_array_(i_flip, 3); // xyy
          d3_array_(i, 4) = d3_array_(i_flip, 4); // xyz
          d3_array_(i, 5) = -d3_array_(i_flip, 5); // xzz
          d3_array_(i, 6) = -d3_array_(i_flip, 6); // yyy
          d3_array_(i, 7) = d3_array_(i_flip, 7); // yyz
          d3_array_(i, 8) = -d3_array_(i_flip, 8); // yzz
          d3_array_(i, 9) = d3_array_(i_flip, 9); // zzz
          
          d4_array_(i, 0) = d4_array_(i_flip, 0); // xxxx
          d4_array_(i, 1) = d4_array_(i_flip, 1); // xxxy
          d4_array_(i, 2) = -d4_array_(i_flip, 2); // xxxz
          d4_array_(i, 3) = d4_array_(i_flip, 3); // xxyy
          d4_array_(i, 4) = -d4_array_(i_flip, 4); // xxyz
          d4_array_(i, 5) = d4_array_(i_flip, 5); // xxzz
          d4_array_(i, 6) = d4_array_(i_flip, 6); // xyyy
          d4_array_(i, 7) = -d4_array_(i_flip, 7); // xyyz
          d4_array_(i, 8) = d4_array_(i_flip, 8); // xyzz
          d4_array_(i, 9) = -d4_array_(i_flip, 9); // xzzz
          d4_array_(i, 10) = d4_array_(i_flip, 10); // yyyy
          d4_array_(i, 11) = -d4_array_(i_flip, 11); // yyyz
          d4_array_(i, 12) = d4_array_(i_flip, 12); // yyzz
          d4_array_(i, 13) = -d4_array_(i_flip, 13); // yzzz
          d4_array_(i, 14) = d4_array_(i_flip, 14); // zzzz
          
          d5_array_(i, 0) = -d5_array_(i_flip, 0); // xxxxx
          d5_array_(i, 1) = -d5_array_(i_flip, 1); // xxxxy
          d5_array_(i, 2) = d5_array_(i_flip, 2); // xxxxz
          d5_array_(i, 3) = -d5_array_(i_flip, 3); // xxxyy
          d5_array_(i, 4) = d5_array_(i_flip, 4); // xxxyz
          d5_array_(i, 5) = -d5_array_(i_flip, 5); // xxxzz
          d5_array_(i, 6) = -d5_array_(i_flip, 6); // xxyyy
          d5_array_(i, 7) = d5_array_(i_flip, 7); // xxyyz
          d5_array_(i, 8) = -d5_array_(i_flip, 8); // xxyzz
          d5_array_(i, 9) = d5_array_(i_flip, 9); // xxzzz
          d5_array_(i, 10) = -d5_array_(i_flip, 10); // xyyyy
          d5_array_(i, 11) = d5_array_(i_flip, 11); // xyyyz
          d5_array_(i, 12) = -d5_array_(i_flip, 12); // xyyzz
          d5_array_(i, 13) = d5_array_(i_flip, 13); // xyzzz
          d5_array_(i, 14) = -d5_array_(i_flip, 14); // xzzzz
          d5_array_(i, 15) = -d5_array_(i_flip, 15); // yyyyy
          d5_array_(i, 16) = d5_array_(i_flip, 16); // yyyyz
          d5_array_(i, 17) = -d5_array_(i_flip, 17); // yyyzz
          d5_array_(i, 18) = d5_array_(i_flip, 18); // yyzzz
          d5_array_(i, 19) = -d5_array_(i_flip, 19); // yzzzz
          d5_array_(i, 20) = d5_array_(i_flip, 20); // zzzzz
          
          d6_array_(i, 0) = d6_array_(i_flip, 0); // xxxxxx
          d6_array_(i, 1) = d6_array_(i_flip, 1); // xxxxxy
          d6_array_(i, 2) = -d6_array_(i_flip, 2); // xxxxxz
          d6_array_(i, 3) = d6_array_(i_flip, 3); // xxxxyy
          d6_array_(i, 4) = -d6_array_(i_flip, 4); // xxxxyz
          d6_array_(i, 5) = d6_array_(i_flip, 5); // xxxxzz
          d6_array_(i, 6) = d6_array_(i_flip, 6); // xxxyyy
          d6_array_(i, 7) = -d6_array_(i_flip, 7); // xxxyyz
          d6_array_(i, 8) = d6_array_(i_flip, 8); // xxxyzz
          d6_array_(i, 9) = -d6_array_(i_flip, 9); // xxxzzz
          d6_array_(i, 10) = d6_array_(i_flip, 10); // xxyyyy
          d6_array_(i, 11) = -d6_array_(i_flip, 11); // xxyyyz
          d6_array_(i, 12) = d6_array_(i_flip, 12); // xxyyzz
          d6_array_(i, 13) = -d6_array_(i_flip, 13); // xxyzzz
          d6_array_(i, 14) = d6_array_(i_flip, 14); // xxzzzz
          d6_array_(i, 15) = d6_array_(i_flip, 15); // xyyyyy
          d6_array_(i, 16) = -d6_array_(i_flip, 16); // xyyyyz
          d6_array_(i, 17) = d6_array_(i_flip, 17); // xyyyzz
          d6_array_(i, 18) = -d6_array_(i_flip, 18); // xyyzzz
          d6_array_(i, 19) = d6_array_(i_flip, 19); // xyzzzz
          d6_array_(i, 20) = -d6_array_(i_flip, 20); // xzzzzz
          d6_array_(i, 21) = d6_array_(i_flip, 21); // yyyyyy
          d6_array_(i, 22) = -d6_array_(i_flip, 22); // yyyyyz
          d6_array_(i, 23) = d6_array_(i_flip, 23); // yyyyzz
          d6_array_(i, 24) = -d6_array_(i_flip, 24); // yyyzzz
          d6_array_(i, 25) = d6_array_(i_flip, 25); // yyzzzz
          d6_array_(i, 26) = -d6_array_(i_flip, 26); // yzzzzz
          d6_array_(i, 27) = d6_array_(i_flip, 27); // zzzzzz
      }
    }
  }

  for (int iz = midz+1; iz < interp_zpoints_; iz++) {
    for (int iy = 0; iy <= midy; iy++) {
      for (int ix = midx+1; ix < interp_xpoints_; ix++) {
          
          int i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);
          int ix_flip = interp_xpoints_ - 1 - ix;
          int iz_flip = interp_zpoints_ - 1 - iz;
          int i_flip = ix_flip + interp_xpoints_ * (iy + iz_flip * interp_ypoints_);

          d1_array_(i, 0) = -d1_array_(i_flip, 0); // x
          d1_array_(i, 1) = d1_array_(i_flip, 1); // y
          d1_array_(i, 2) = -d1_array_(i_flip, 2); // z
          
          d2_array_(i, 0) = d2_array_(i_flip, 0); // xx
          d2_array_(i, 1) = -d2_array_(i_flip, 1); // xy
          d2_array_(i, 2) = d2_array_(i_flip, 2); // xz
          d2_array_(i, 3) = d2_array_(i_flip, 3); // yy
          d2_array_(i, 4) = -d2_array_(i_flip, 4); // yz
          d2_array_(i, 5) = d2_array_(i_flip, 5); // zz
          
          d3_array_(i, 0) = -d3_array_(i_flip, 0); // xxx
          d3_array_(i, 1) = d3_array_(i_flip, 1); // xxy
          d3_array_(i, 2) = -d3_array_(i_flip, 2); // xxz
          d3_array_(i, 3) = -d3_array_(i_flip, 3); // xyy
          d3_array_(i, 4) = d3_array_(i_flip, 4); // xyz
          d3_array_(i, 5) = -d3_array_(i_flip, 5); // xzz
          d3_array_(i, 6) = d3_array_(i_flip, 6); // yyy
          d3_array_(i, 7) = -d3_array_(i_flip, 7); // yyz
          d3_array_(i, 8) = d3_array_(i_flip, 8); // yzz
          d3_array_(i, 9) = -d3_array_(i_flip, 9); // zzz
          
          d4_array_(i, 0) = d4_array_(i_flip, 0); // xxxx
          d4_array_(i, 1) = -d4_array_(i_flip, 1); // xxxy
          d4_array_(i, 2) = d4_array_(i_flip, 2); // xxxz
          d4_array_(i, 3) = d4_array_(i_flip, 3); // xxyy
          d4_array_(i, 4) = -d4_array_(i_flip, 4); // xxyz
          d4_array_(i, 5) = d4_array_(i_flip, 5); // xxzz
          d4_array_(i, 6) = -d4_array_(i_flip, 6); // xyyy
          d4_array_(i, 7) = d4_array_(i_flip, 7); // xyyz
          d4_array_(i, 8) = -d4_array_(i_flip, 8); // xyzz
          d4_array_(i, 9) = d4_array_(i_flip, 9); // xzzz
          d4_array_(i, 10) = d4_array_(i_flip, 10); // yyyy
          d4_array_(i, 11) = -d4_array_(i_flip, 11); // yyyz
          d4_array_(i, 12) = d4_array_(i_flip, 12); // yyzz
          d4_array_(i, 13) = -d4_array_(i_flip, 13); // yzzz
          d4_array_(i, 14) = d4_array_(i_flip, 14); // zzzz
          
          d5_array_(i, 0) = -d5_array_(i_flip, 0); // xxxxx
          d5_array_(i, 1) = d5_array_(i_flip, 1); // xxxxy
          d5_array_(i, 2) = -d5_array_(i_flip, 2); // xxxxz
          d5_array_(i, 3) = -d5_array_(i_flip, 3); // xxxyy
          d5_array_(i, 4) = d5_array_(i_flip, 4); // xxxyz
          d5_array_(i, 5) = -d5_array_(i_flip, 5); // xxxzz
          d5_array_(i, 6) = d5_array_(i_flip, 6); // xxyyy
          d5_array_(i, 7) = -d5_array_(i_flip, 7); // xxyyz
          d5_array_(i, 8) = d5_array_(i_flip, 8); // xxyzz
          d5_array_(i, 9) = -d5_array_(i_flip, 9); // xxzzz
          d5_array_(i, 10) = -d5_array_(i_flip, 10); // xyyyy
          d5_array_(i, 11) = d5_array_(i_flip, 11); // xyyyz
          d5_array_(i, 12) = -d5_array_(i_flip, 12); // xyyzz
          d5_array_(i, 13) = d5_array_(i_flip, 13); // xyzzz
          d5_array_(i, 14) = -d5_array_(i_flip, 14); // xzzzz
          d5_array_(i, 15) = d5_array_(i_flip, 15); // yyyyy
          d5_array_(i, 16) = -d5_array_(i_flip, 16); // yyyyz
          d5_array_(i, 17) = d5_array_(i_flip, 17); // yyyzz
          d5_array_(i, 18) = -d5_array_(i_flip, 18); // yyzzz
          d5_array_(i, 19) = d5_array_(i_flip, 19); // yzzzz
          d5_array_(i, 20) = -d5_array_(i_flip, 20); // zzzzz
          
          d6_array_(i, 0) = d6_array_(i_flip, 0); // xxxxxx
          d6_array_(i, 1) = -d6_array_(i_flip, 1); // xxxxxy
          d6_array_(i, 2) = d6_array_(i_flip, 2); // xxxxxz
          d6_array_(i, 3) = d6_array_(i_flip, 3); // xxxxyy
          d6_array_(i, 4) = -d6_array_(i_flip, 4); // xxxxyz
          d6_array_(i, 5) = d6_array_(i_flip, 5); // xxxxzz
          d6_array_(i, 6) = -d6_array_(i_flip, 6); // xxxyyy
          d6_array_(i, 7) = d6_array_(i_flip, 7); // xxxyyz
          d6_array_(i, 8) = -d6_array_(i_flip, 8); // xxxyzz
          d6_array_(i, 9) = d6_array_(i_flip, 9); // xxxzzz
          d6_array_(i, 10) = d6_array_(i_flip, 10); // xxyyyy
          d6_array_(i, 11) = -d6_array_(i_flip, 11); // xxyyyz
          d6_array_(i, 12) = d6_array_(i_flip, 12); // xxyyzz
          d6_array_(i, 13) = -d6_array_(i_flip, 13); // xxyzzz
          d6_array_(i, 14) = d6_array_(i_flip, 14); // xxzzzz
          d6_array_(i, 15) = -d6_array_(i_flip, 15); // xyyyyy
          d6_array_(i, 16) = d6_array_(i_flip, 16); // xyyyyz
          d6_array_(i, 17) = -d6_array_(i_flip, 17); // xyyyzz
          d6_array_(i, 18) = d6_array_(i_flip, 18); // xyyzzz
          d6_array_(i, 19) = -d6_array_(i_flip, 19); // xyzzzz
          d6_array_(i, 20) = d6_array_(i_flip, 20); // xzzzzz
          d6_array_(i, 21) = d6_array_(i_flip, 21); // yyyyyy
          d6_array_(i, 22) = -d6_array_(i_flip, 22); // yyyyyz
          d6_array_(i, 23) = d6_array_(i_flip, 23); // yyyyzz
          d6_array_(i, 24) = -d6_array_(i_flip, 24); // yyyzzz
          d6_array_(i, 25) = d6_array_(i_flip, 25); // yyzzzz
          d6_array_(i, 26) = -d6_array_(i_flip, 26); // yzzzzz
          d6_array_(i, 27) = d6_array_(i_flip, 27); // zzzzzz
      }
    }
  }

  for (int iz = midz+1; iz < interp_zpoints_; iz++) {
    for (int iy = midy+1; iy < interp_ypoints_; iy++) {
      for (int ix = 0; ix <= midx; ix++) {
          
          int i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);
          int iy_flip = interp_ypoints_ - 1 - iy;
          int iz_flip = interp_zpoints_ - 1 - iz;
          int i_flip = ix + interp_xpoints_ * (iy_flip + iz_flip * interp_ypoints_);

          d1_array_(i, 0) = d1_array_(i_flip, 0); // x
          d1_array_(i, 1) = -d1_array_(i_flip, 1); // y
          d1_array_(i, 2) = -d1_array_(i_flip, 2); // z
          
          d2_array_(i, 0) = d2_array_(i_flip, 0); // xx
          d2_array_(i, 1) = -d2_array_(i_flip, 1); // xy
          d2_array_(i, 2) = -d2_array_(i_flip, 2); // xz
          d2_array_(i, 3) = d2_array_(i_flip, 3); // yy
          d2_array_(i, 4) = d2_array_(i_flip, 4); // yz
          d2_array_(i, 5) = d2_array_(i_flip, 5); // zz
          
          d3_array_(i, 0) = d3_array_(i_flip, 0); // xxx
          d3_array_(i, 1) = -d3_array_(i_flip, 1); // xxy
          d3_array_(i, 2) = -d3_array_(i_flip, 2); // xxz
          d3_array_(i, 3) = d3_array_(i_flip, 3); // xyy
          d3_array_(i, 4) = d3_array_(i_flip, 4); // xyz
          d3_array_(i, 5) = d3_array_(i_flip, 5); // xzz
          d3_array_(i, 6) = -d3_array_(i_flip, 6); // yyy
          d3_array_(i, 7) = -d3_array_(i_flip, 7); // yyz
          d3_array_(i, 8) = -d3_array_(i_flip, 8); // yzz
          d3_array_(i, 9) = -d3_array_(i_flip, 9); // zzz
          
          d4_array_(i, 0) = d4_array_(i_flip, 0); // xxxx
          d4_array_(i, 1) = -d4_array_(i_flip, 1); // xxxy
          d4_array_(i, 2) = -d4_array_(i_flip, 2); // xxxz
          d4_array_(i, 3) = d4_array_(i_flip, 3); // xxyy
          d4_array_(i, 4) = d4_array_(i_flip, 4); // xxyz
          d4_array_(i, 5) = d4_array_(i_flip, 5); // xxzz
          d4_array_(i, 6) = -d4_array_(i_flip, 6); // xyyy
          d4_array_(i, 7) = -d4_array_(i_flip, 7); // xyyz
          d4_array_(i, 8) = -d4_array_(i_flip, 8); // xyzz
          d4_array_(i, 9) = -d4_array_(i_flip, 9); // xzzz
          d4_array_(i, 10) = d4_array_(i_flip, 10); // yyyy
          d4_array_(i, 11) = d4_array_(i_flip, 11); // yyyz
          d4_array_(i, 12) = d4_array_(i_flip, 12); // yyzz
          d4_array_(i, 13) = d4_array_(i_flip, 13); // yzzz
          d4_array_(i, 14) = d4_array_(i_flip, 14); // zzzz
          
          d5_array_(i, 0) = d5_array_(i_flip, 0); // xxxxx
          d5_array_(i, 1) = -d5_array_(i_flip, 1); // xxxxy
          d5_array_(i, 2) = -d5_array_(i_flip, 2); // xxxxz
          d5_array_(i, 3) = d5_array_(i_flip, 3); // xxxyy
          d5_array_(i, 4) = d5_array_(i_flip, 4); // xxxyz
          d5_array_(i, 5) = d5_array_(i_flip, 5); // xxxzz
          d5_array_(i, 6) = -d5_array_(i_flip, 6); // xxyyy
          d5_array_(i, 7) = -d5_array_(i_flip, 7); // xxyyz
          d5_array_(i, 8) = -d5_array_(i_flip, 8); // xxyzz
          d5_array_(i, 9) = -d5_array_(i_flip, 9); // xxzzz
          d5_array_(i, 10) = d5_array_(i_flip, 10); // xyyyy
          d5_array_(i, 11) = d5_array_(i_flip, 11); // xyyyz
          d5_array_(i, 12) = d5_array_(i_flip, 12); // xyyzz
          d5_array_(i, 13) = d5_array_(i_flip, 13); // xyzzz
          d5_array_(i, 14) = d5_array_(i_flip, 14); // xzzzz
          d5_array_(i, 15) = -d5_array_(i_flip, 15); // yyyyy
          d5_array_(i, 16) = -d5_array_(i_flip, 16); // yyyyz
          d5_array_(i, 17) = -d5_array_(i_flip, 17); // yyyzz
          d5_array_(i, 18) = -d5_array_(i_flip, 18); // yyzzz
          d5_array_(i, 19) = -d5_array_(i_flip, 19); // yzzzz
          d5_array_(i, 20) = -d5_array_(i_flip, 20); // zzzzz
          
          d6_array_(i, 0) = d6_array_(i_flip, 0); // xxxxxx
          d6_array_(i, 1) = -d6_array_(i_flip, 1); // xxxxxy
          d6_array_(i, 2) = -d6_array_(i_flip, 2); // xxxxxz
          d6_array_(i, 3) = d6_array_(i_flip, 3); // xxxxyy
          d6_array_(i, 4) = d6_array_(i_flip, 4); // xxxxyz
          d6_array_(i, 5) = d6_array_(i_flip, 5); // xxxxzz
          d6_array_(i, 6) = -d6_array_(i_flip, 6); // xxxyyy
          d6_array_(i, 7) = -d6_array_(i_flip, 7); // xxxyyz
          d6_array_(i, 8) = -d6_array_(i_flip, 8); // xxxyzz
          d6_array_(i, 9) = -d6_array_(i_flip, 9); // xxxzzz
          d6_array_(i, 10) = d6_array_(i_flip, 10); // xxyyyy
          d6_array_(i, 11) = d6_array_(i_flip, 11); // xxyyyz
          d6_array_(i, 12) = d6_array_(i_flip, 12); // xxyyzz
          d6_array_(i, 13) = d6_array_(i_flip, 13); // xxyzzz
          d6_array_(i, 14) = d6_array_(i_flip, 14); // xxzzzz
          d6_array_(i, 15) = -d6_array_(i_flip, 15); // xyyyyy
          d6_array_(i, 16) = -d6_array_(i_flip, 16); // xyyyyz
          d6_array_(i, 17) = -d6_array_(i_flip, 17); // xyyyzz
          d6_array_(i, 18) = -d6_array_(i_flip, 18); // xyyzzz
          d6_array_(i, 19) = -d6_array_(i_flip, 19); // xyzzzz
          d6_array_(i, 20) = -d6_array_(i_flip, 20); // xzzzzz
          d6_array_(i, 21) = d6_array_(i_flip, 21); // yyyyyy
          d6_array_(i, 22) = d6_array_(i_flip, 22); // yyyyyz
          d6_array_(i, 23) = d6_array_(i_flip, 23); // yyyyzz
          d6_array_(i, 24) = d6_array_(i_flip, 24); // yyyzzz
          d6_array_(i, 25) = d6_array_(i_flip, 25); // yyzzzz
          d6_array_(i, 26) = d6_array_(i_flip, 26); // yzzzzz
          d6_array_(i, 27) = d6_array_(i_flip, 27); // zzzzzz
      }
    }
  }

  for (int iz = midz+1; iz < interp_zpoints_; iz++) {
    for (int iy = midy+1; iy < interp_ypoints_; iy++) {
      for (int ix = midx+1; ix < interp_xpoints_; ix++) {
          
          int i = ix + interp_xpoints_ * (iy + iz * interp_ypoints_);
          int ix_flip = interp_xpoints_ - 1 - ix;
          int iy_flip = interp_ypoints_ - 1 - iy;
          int iz_flip = interp_zpoints_ - 1 - iz;
          int i_flip = ix_flip + interp_xpoints_ * (iy_flip + iz_flip * interp_ypoints_);

          d1_array_(i, 0) = -d1_array_(i_flip, 0); // x
          d1_array_(i, 1) = -d1_array_(i_flip, 1); // y
          d1_array_(i, 2) = -d1_array_(i_flip, 2); // z
          
          d2_array_(i, 0) = d2_array_(i_flip, 0); // xx
          d2_array_(i, 1) = d2_array_(i_flip, 1); // xy
          d2_array_(i, 2) = d2_array_(i_flip, 2); // xz
          d2_array_(i, 3) = d2_array_(i_flip, 3); // yy
          d2_array_(i, 4) = d2_array_(i_flip, 4); // yz
          d2_array_(i, 5) = d2_array_(i_flip, 5); // zz
          
          d3_array_(i, 0) = -d3_array_(i_flip, 0); // xxx
          d3_array_(i, 1) = -d3_array_(i_flip, 1); // xxy
          d3_array_(i, 2) = -d3_array_(i_flip, 2); // xxz
          d3_array_(i, 3) = -d3_array_(i_flip, 3); // xyy
          d3_array_(i, 4) = -d3_array_(i_flip, 4); // xyz
          d3_array_(i, 5) = -d3_array_(i_flip, 5); // xzz
          d3_array_(i, 6) = -d3_array_(i_flip, 6); // yyy
          d3_array_(i, 7) = -d3_array_(i_flip, 7); // yyz
          d3_array_(i, 8) = -d3_array_(i_flip, 8); // yzz
          d3_array_(i, 9) = -d3_array_(i_flip, 9); // zzz
          
          d4_array_(i, 0) = d4_array_(i_flip, 0); // xxxx
          d4_array_(i, 1) = d4_array_(i_flip, 1); // xxxy
          d4_array_(i, 2) = d4_array_(i_flip, 2); // xxxz
          d4_array_(i, 3) = d4_array_(i_flip, 3); // xxyy
          d4_array_(i, 4) = d4_array_(i_flip, 4); // xxyz
          d4_array_(i, 5) = d4_array_(i_flip, 5); // xxzz
          d4_array_(i, 6) = d4_array_(i_flip, 6); // xyyy
          d4_array_(i, 7) = d4_array_(i_flip, 7); // xyyz
          d4_array_(i, 8) = d4_array_(i_flip, 8); // xyzz
          d4_array_(i, 9) = d4_array_(i_flip, 9); // xzzz
          d4_array_(i, 10) = d4_array_(i_flip, 10); // yyyy
          d4_array_(i, 11) = d4_array_(i_flip, 11); // yyyz
          d4_array_(i, 12) = d4_array_(i_flip, 12); // yyzz
          d4_array_(i, 13) = d4_array_(i_flip, 13); // yzzz
          d4_array_(i, 14) = d4_array_(i_flip, 14); // zzzz
          
          d5_array_(i, 0) = -d5_array_(i_flip, 0); // xxxxx
          d5_array_(i, 1) = -d5_array_(i_flip, 1); // xxxxy
          d5_array_(i, 2) = -d5_array_(i_flip, 2); // xxxxz
          d5_array_(i, 3) = -d5_array_(i_flip, 3); // xxxyy
          d5_array_(i, 4) = -d5_array_(i_flip, 4); // xxxyz
          d5_array_(i, 5) = -d5_array_(i_flip, 5); // xxxzz
          d5_array_(i, 6) = -d5_array_(i_flip, 6); // xxyyy
          d5_array_(i, 7) = -d5_array_(i_flip, 7); // xxyyz
          d5_array_(i, 8) = -d5_array_(i_flip, 8); // xxyzz
          d5_array_(i, 9) = -d5_array_(i_flip, 9); // xxzzz
          d5_array_(i, 10) = -d5_array_(i_flip, 10); // xyyyy
          d5_array_(i, 11) = -d5_array_(i_flip, 11); // xyyyz
          d5_array_(i, 12) = -d5_array_(i_flip, 12); // xyyzz
          d5_array_(i, 13) = -d5_array_(i_flip, 13); // xyzzz
          d5_array_(i, 14) = -d5_array_(i_flip, 14); // xzzzz
          d5_array_(i, 15) = -d5_array_(i_flip, 15); // yyyyy
          d5_array_(i, 16) = -d5_array_(i_flip, 16); // yyyyz
          d5_array_(i, 17) = -d5_array_(i_flip, 17); // yyyzz
          d5_array_(i, 18) = -d5_array_(i_flip, 18); // yyzzz
          d5_array_(i, 19) = -d5_array_(i_flip, 19); // yzzzz
          d5_array_(i, 20) = -d5_array_(i_flip, 20); // zzzzz
          
          d6_array_(i, 0) = d6_array_(i_flip, 0); // xxxxxx
          d6_array_(i, 1) = d6_array_(i_flip, 1); // xxxxxy
          d6_array_(i, 2) = d6_array_(i_flip, 2); // xxxxxz
          d6_array_(i, 3) = d6_array_(i_flip, 3); // xxxxyy
          d6_array_(i, 4) = d6_array_(i_flip, 4); // xxxxyz
          d6_array_(i, 5) = d6_array_(i_flip, 5); // xxxxzz
          d6_array_(i, 6) = d6_array_(i_flip, 6); // xxxyyy
          d6_array_(i, 7) = d6_array_(i_flip, 7); // xxxyyz
          d6_array_(i, 8) = d6_array_(i_flip, 8); // xxxyzz
          d6_array_(i, 9) = d6_array_(i_flip, 9); // xxxzzz
          d6_array_(i, 10) = d6_array_(i_flip, 10); // xxyyyy
          d6_array_(i, 11) = d6_array_(i_flip, 11); // xxyyyz
          d6_array_(i, 12) = d6_array_(i_flip, 12); // xxyyzz
          d6_array_(i, 13) = d6_array_(i_flip, 13); // xxyzzz
          d6_array_(i, 14) = d6_array_(i_flip, 14); // xxzzzz
          d6_array_(i, 15) = d6_array_(i_flip, 15); // xyyyyy
          d6_array_(i, 16) = d6_array_(i_flip, 16); // xyyyyz
          d6_array_(i, 17) = d6_array_(i_flip, 17); // xyyyzz
          d6_array_(i, 18) = d6_array_(i_flip, 18); // xyyzzz
          d6_array_(i, 19) = d6_array_(i_flip, 19); // xyzzzz
          d6_array_(i, 20) = d6_array_(i_flip, 20); // xzzzzz
          d6_array_(i, 21) = d6_array_(i_flip, 21); // yyyyyy
          d6_array_(i, 22) = d6_array_(i_flip, 22); // yyyyyz
          d6_array_(i, 23) = d6_array_(i_flip, 23); // yyyyzz
          d6_array_(i, 24) = d6_array_(i_flip, 24); // yyyzzz
          d6_array_(i, 25) = d6_array_(i_flip, 25); // yyzzzz
          d6_array_(i, 26) = d6_array_(i_flip, 26); // yzzzzz
          d6_array_(i, 27) = d6_array_(i_flip, 27); // zzzzzz
      }
    }
  }

}


/********************************************************************************************/
/*  Evaluate the Taylor series used to interpolate the Ewald derivatives */


// compute the Taylor series to interpolate derivative tensors from interpolation points to (x,y,z).


/* Note: interp_d0 is not necessary, since our Taylor series start at d1 */
// double EnzoEwald::interp_d0(double x, double y, double z) throw()
// {
//   int i;
//   double interp_x, interp_y, interp_z;
//   find_nearest_interp_point(x, y, z, &interp_x, &interp_y, &interp_z, &i);

//   std::array<double, 3> delta_r = {x - interp_x, y - interp_y, z - interp_z};
//   std::array<double, 6> delta_r2 = EnzoMethodMultipole::outer_11_(delta_r, delta_r);
//   std::array<double, 10> delta_r3 = EnzoMethodMultipole::outer_12_(delta_r, delta_r2);

//   double zeroth_term = d0_array_[i];
//   double first_term = EnzoMethodMultipole::dot_11_(delta_r, d1_array_[i]);
//   double second_term = 0.5 * EnzoMethodMultipole::dot_22_(delta_r2, d2_array_[i]);
//   double third_term = 1.0/6.0 * EnzoMethodMultipole::dot_33_(delta_r3, d3_array_[i]);

//   return zeroth_term + first_term + second_term + third_term;
  
// }

std::array<double, 3> EnzoEwald::interp_d1(double x, double y, double z) throw()
{

  int i;
  double interp_x, interp_y, interp_z;
  find_nearest_interp_point(x, y, z, &interp_x, &interp_y, &interp_z, &i);

  // CkPrintf("x, y, z: %f, %f, %f\n", x, y, z);
  // CkPrintf("interp points: %f, %f, %f\n", interp_x, interp_y, interp_z);
  // CkPrintf("i: %d\n", i);
    
  std::array<double, 3> delta_r = {x - interp_x, y - interp_y, z - interp_z};
  std::array<double, 6> delta_r2 = EnzoMethodMultipole::outer_11_(delta_r, delta_r);
  std::array<double, 10> delta_r3 = EnzoMethodMultipole::outer_12_(delta_r, delta_r2);

  // std::array<double, 3> zeroth_term = d1_array_[i]; // d1_array_.subarray(i) 
  // std::array<double, 3> first_term = EnzoMethodMultipole::dot_12_(delta_r, d2_array_[i]);
  // std::array<double, 3> second_term = EnzoMethodMultipole::dot_scalar_1_(0.5, EnzoMethodMultipole::dot_23_(delta_r2, d3_array_[i]));
  // std::array<double, 3> third_term = EnzoMethodMultipole::dot_scalar_1_(1.0/6.0, EnzoMethodMultipole::dot_34_(delta_r3, d4_array_[i]));

  // std::array<double, 3> zero_plus_one = EnzoMethodMultipole::add_11_(zeroth_term, first_term);
  // std::array<double, 3> two_plus_three = EnzoMethodMultipole::add_11_(second_term, third_term);

  std::array<double, 3> first_term;
  std::array<double, 3> second_term;
  std::array<double, 3> third_term;
  std::array<double, 3> fourth_term = {0.0, 0.0, 0.0};  // only needed if order_ >= 4

  first_term[0] = delta_r[0]*d2_array_(i,0) + delta_r[1]*d2_array_(i,1) + delta_r[2]*d2_array_(i,2);   
  first_term[1] = delta_r[0]*d2_array_(i,1) + delta_r[1]*d2_array_(i,3) + delta_r[2]*d2_array_(i,4);   
  first_term[2] = delta_r[0]*d2_array_(i,2) + delta_r[1]*d2_array_(i,4) + delta_r[2]*d2_array_(i,5);

  second_term[0] = delta_r2[0]*d3_array_(i,0) + 2*delta_r2[1]*d3_array_(i,1) + 2*delta_r2[2]*d3_array_(i,2) + delta_r2[3]*d3_array_(i,3) + 2*delta_r2[4]*d3_array_(i,4) + delta_r2[5]*d3_array_(i,5);
  second_term[1] = delta_r2[0]*d3_array_(i,1) + 2*delta_r2[1]*d3_array_(i,3) + 2*delta_r2[2]*d3_array_(i,4) + delta_r2[3]*d3_array_(i,6) + 2*delta_r2[4]*d3_array_(i,7) + delta_r2[5]*d3_array_(i,8);
  second_term[2] = delta_r2[0]*d3_array_(i,2) + 2*delta_r2[1]*d3_array_(i,4) + 2*delta_r2[2]*d3_array_(i,5) + delta_r2[3]*d3_array_(i,7) + 2*delta_r2[4]*d3_array_(i,8) + delta_r2[5]*d3_array_(i,9);

  third_term[0] = delta_r3[0]*d4_array_(i,0) + 3*delta_r3[1]*d4_array_(i,1) + 3*delta_r3[2]*d4_array_(i,2) + 3*delta_r3[3]*d4_array_(i,3) + 6*delta_r3[4]*d4_array_(i,4) + 3*delta_r3[5]*d4_array_(i,5) + delta_r3[6]*d4_array_(i,6) + 3*delta_r3[7]*d4_array_(i,7) + 3*delta_r3[8]*d4_array_(i,8) + delta_r3[9]*d4_array_(i,9);
  third_term[1] = delta_r3[0]*d4_array_(i,1) + 3*delta_r3[1]*d4_array_(i,3) + 3*delta_r3[2]*d4_array_(i,4) + 3*delta_r3[3]*d4_array_(i,6) + 6*delta_r3[4]*d4_array_(i,7) + 3*delta_r3[5]*d4_array_(i,8) + delta_r3[6]*d4_array_(i,10) + 3*delta_r3[7]*d4_array_(i,11) + 3*delta_r3[8]*d4_array_(i,12) + delta_r3[9]*d4_array_(i,13);
  third_term[2] = delta_r3[0]*d4_array_(i,2) + 3*delta_r3[1]*d4_array_(i,4) + 3*delta_r3[2]*d4_array_(i,5) + 3*delta_r3[3]*d4_array_(i,7) + 6*delta_r3[4]*d4_array_(i,8) + 3*delta_r3[5]*d4_array_(i,9) + delta_r3[6]*d4_array_(i,11) + 3*delta_r3[7]*d4_array_(i,12) + 3*delta_r3[8]*d4_array_(i,13) + delta_r3[9]*d4_array_(i,14);

  // the fourth-order term needs the 5th derivatives of the potential
  if (order_ >= 4) {
    std::array<double, 15> delta_r4 = EnzoMethodMultipole::outer_13_(delta_r, delta_r3);

    fourth_term[0] = delta_r4[0]*d5_array_(i,0) + 4*delta_r4[1]*d5_array_(i,1) + 4*delta_r4[2]*d5_array_(i,2) + 6*delta_r4[3]*d5_array_(i,3) + 12*delta_r4[4]*d5_array_(i,4) + 6*delta_r4[5]*d5_array_(i,5) + 4*delta_r4[6]*d5_array_(i,6) + 12*delta_r4[7]*d5_array_(i,7) + 12*delta_r4[8]*d5_array_(i,8) + 4*delta_r4[9]*d5_array_(i,9) + delta_r4[10]*d5_array_(i,10) + 4*delta_r4[11]*d5_array_(i,11) + 6*delta_r4[12]*d5_array_(i,12) + 4*delta_r4[13]*d5_array_(i,13) + delta_r4[14]*d5_array_(i,14);
    fourth_term[1] = delta_r4[0]*d5_array_(i,1) + 4*delta_r4[1]*d5_array_(i,3) + 4*delta_r4[2]*d5_array_(i,4) + 6*delta_r4[3]*d5_array_(i,6) + 12*delta_r4[4]*d5_array_(i,7) + 6*delta_r4[5]*d5_array_(i,8) + 4*delta_r4[6]*d5_array_(i,10) + 12*delta_r4[7]*d5_array_(i,11) + 12*delta_r4[8]*d5_array_(i,12) + 4*delta_r4[9]*d5_array_(i,13) + delta_r4[10]*d5_array_(i,15) + 4*delta_r4[11]*d5_array_(i,16) + 6*delta_r4[12]*d5_array_(i,17) + 4*delta_r4[13]*d5_array_(i,18) + delta_r4[14]*d5_array_(i,19);
    fourth_term[2] = delta_r4[0]*d5_array_(i,2) + 4*delta_r4[1]*d5_array_(i,4) + 4*delta_r4[2]*d5_array_(i,5) + 6*delta_r4[3]*d5_array_(i,7) + 12*delta_r4[4]*d5_array_(i,8) + 6*delta_r4[5]*d5_array_(i,9) + 4*delta_r4[6]*d5_array_(i,11) + 12*delta_r4[7]*d5_array_(i,12) + 12*delta_r4[8]*d5_array_(i,13) + 4*delta_r4[9]*d5_array_(i,14) + delta_r4[10]*d5_array_(i,16) + 4*delta_r4[11]*d5_array_(i,17) + 6*delta_r4[12]*d5_array_(i,18) + 4*delta_r4[13]*d5_array_(i,19) + delta_r4[14]*d5_array_(i,20);
  }

  // CkPrintf("d1 ewald zeroth_term: %f, %f, %f\n", d1_array_(i, 0), d1_array_(i, 1), d1_array_(i, 2));
  // CkPrintf("d1 ewald first_term: %f, %f, %f\n", first_term[0], first_term[1], first_term[2]);
  // CkPrintf("d1 ewald second_term: %f, %f, %f\n", 0.5*second_term[0], 0.5*second_term[1], 0.5*second_term[2]);
  // CkPrintf("d1 ewald third_term: %f, %f, %f\n", 1.0/6.0*third_term[0], 1.0/6.0*third_term[1], 1.0/6.0*third_term[2]);


  std::array<double, 3> sum; 
  sum[0] = d1_array_(i, 0) + first_term[0] + 0.5*second_term[0] + 1.0/6.0*third_term[0] + 1.0/24.0*fourth_term[0];
  sum[1] = d1_array_(i, 1) + first_term[1] + 0.5*second_term[1] + 1.0/6.0*third_term[1] + 1.0/24.0*fourth_term[1];
  sum[2] = d1_array_(i, 2) + first_term[2] + 0.5*second_term[2] + 1.0/6.0*third_term[2] + 1.0/24.0*fourth_term[2];

  return sum;
  
}

std::array<double, 6> EnzoEwald::interp_d2(double x, double y, double z) throw()
{

  int i;
  double interp_x, interp_y, interp_z;
  find_nearest_interp_point(x, y, z, &interp_x, &interp_y, &interp_z, &i);

  // CelloView<double, 1> d2_slice_ = d2(0, 0, 0);
  // CkPrintf("d2_slice: %f, %f, %f, %f, %f, %f\n", d2_slice_(0), d2_slice_(1), d2_slice_(2), d2_slice_(3), d2_slice_(4), d2_slice_(5)); 

  //CkPrintf("i, interp_x, interp_y: %d, %f, %f\n", i, interp_x, interp_y);

  std::array<double, 3> delta_r = {x - interp_x, y - interp_y, z - interp_z};
  std::array<double, 6> delta_r2 = EnzoMethodMultipole::outer_11_(delta_r, delta_r);
  std::array<double, 10> delta_r3 = EnzoMethodMultipole::outer_12_(delta_r, delta_r2);

  // std::array<double, 6> zeroth_term = d2_array_[i];
  // std::array<double, 6> first_term = EnzoMethodMultipole::dot_13_(delta_r, d3_array_[i]);
  // std::array<double, 6> second_term = EnzoMethodMultipole::dot_scalar_2_(0.5, EnzoMethodMultipole::dot_24_(delta_r2, d4_array_[i]));
  // std::array<double, 6> third_term = EnzoMethodMultipole::dot_scalar_2_(1.0/6.0, EnzoMethodMultipole::dot_35_(delta_r3, d5_array_[i]));
 
  // CkPrintf("0th term: %f\n", zeroth_term[0]);
  // CkPrintf("1st term: %f\n", first_term[0]);
  // CkPrintf("2nd term: %f\n", second_term[0]);
  // CkPrintf("3rd term: %f\n", third_term[0]);

  // std::array<double, 6> zero_plus_one = EnzoMethodMultipole::add_22_(zeroth_term, first_term);
  // std::array<double, 6> two_plus_three = EnzoMethodMultipole::add_22_(second_term, third_term);

  std::array<double, 6> first_term;
  std::array<double, 6> second_term;
  std::array<double, 6> third_term;

  first_term[0] = delta_r[0]*d3_array_(i,0) + delta_r[1]*d3_array_(i,1) + delta_r[2]*d3_array_(i,2);
  first_term[1] = delta_r[0]*d3_array_(i,1) + delta_r[1]*d3_array_(i,3) + delta_r[2]*d3_array_(i,4);
  first_term[2] = delta_r[0]*d3_array_(i,2) + delta_r[1]*d3_array_(i,4) + delta_r[2]*d3_array_(i,5);
  first_term[3] = delta_r[0]*d3_array_(i,3) + delta_r[1]*d3_array_(i,6) + delta_r[2]*d3_array_(i,7);
  first_term[4] = delta_r[0]*d3_array_(i,4) + delta_r[1]*d3_array_(i,7) + delta_r[2]*d3_array_(i,8);
  first_term[5] = delta_r[0]*d3_array_(i,5) + delta_r[1]*d3_array_(i,8) + delta_r[2]*d3_array_(i,9);

  second_term[0] = delta_r2[0]*d4_array_(i,0) + 2*delta_r2[1]*d4_array_(i,1) + 2*delta_r2[2]*d4_array_(i,2) + delta_r2[3]*d4_array_(i,3) + 2*delta_r2[4]*d4_array_(i,4) + delta_r2[5]*d4_array_(i,5);
  second_term[1] = delta_r2[0]*d4_array_(i,1) + 2*delta_r2[1]*d4_array_(i,3) + 2*delta_r2[2]*d4_array_(i,4) + delta_r2[3]*d4_array_(i,6) + 2*delta_r2[4]*d4_array_(i,7) + delta_r2[5]*d4_array_(i,8);
  second_term[2] = delta_r2[0]*d4_array_(i,2) + 2*delta_r2[1]*d4_array_(i,4) + 2*delta_r2[2]*d4_array_(i,5) + delta_r2[3]*d4_array_(i,7) + 2*delta_r2[4]*d4_array_(i,8) + delta_r2[5]*d4_array_(i,9);
  second_term[3] = delta_r2[0]*d4_array_(i,3) + 2*delta_r2[1]*d4_array_(i,6) + 2*delta_r2[2]*d4_array_(i,7) + delta_r2[3]*d4_array_(i,10) + 2*delta_r2[4]*d4_array_(i,11) + delta_r2[5]*d4_array_(i,12);
  second_term[4] = delta_r2[0]*d4_array_(i,4) + 2*delta_r2[1]*d4_array_(i,7) + 2*delta_r2[2]*d4_array_(i,8) + delta_r2[3]*d4_array_(i,11) + 2*delta_r2[4]*d4_array_(i,12) + delta_r2[5]*d4_array_(i,13);
  second_term[5] = delta_r2[0]*d4_array_(i,5) + 2*delta_r2[1]*d4_array_(i,8) + 2*delta_r2[2]*d4_array_(i,9) + delta_r2[3]*d4_array_(i,12) + 2*delta_r2[4]*d4_array_(i,13) + delta_r2[5]*d4_array_(i,14);

  third_term[0] = delta_r3[0]*d5_array_(i,0) + 3*delta_r3[1]*d5_array_(i,1) + 3*delta_r3[2]*d5_array_(i,2) + 3*delta_r3[3]*d5_array_(i,3) + 6*delta_r3[4]*d5_array_(i,4) + 3*delta_r3[5]*d5_array_(i,5) + delta_r3[6]*d5_array_(i,6) + 3*delta_r3[7]*d5_array_(i,7) + 3*delta_r3[8]*d5_array_(i,8) + delta_r3[9]*d5_array_(i,9);
  third_term[1] = delta_r3[0]*d5_array_(i,1) + 3*delta_r3[1]*d5_array_(i,3) + 3*delta_r3[2]*d5_array_(i,4) + 3*delta_r3[3]*d5_array_(i,6) + 6*delta_r3[4]*d5_array_(i,7) + 3*delta_r3[5]*d5_array_(i,8) + delta_r3[6]*d5_array_(i,10) + 3*delta_r3[7]*d5_array_(i,11) + 3*delta_r3[8]*d5_array_(i,12) + delta_r3[9]*d5_array_(i,13);
  third_term[2] = delta_r3[0]*d5_array_(i,2) + 3*delta_r3[1]*d5_array_(i,4) + 3*delta_r3[2]*d5_array_(i,5) + 3*delta_r3[3]*d5_array_(i,7) + 6*delta_r3[4]*d5_array_(i,8) + 3*delta_r3[5]*d5_array_(i,9) + delta_r3[6]*d5_array_(i,11) + 3*delta_r3[7]*d5_array_(i,12) + 3*delta_r3[8]*d5_array_(i,13) + delta_r3[9]*d5_array_(i,14);
  third_term[3] = delta_r3[0]*d5_array_(i,3) + 3*delta_r3[1]*d5_array_(i,6) + 3*delta_r3[2]*d5_array_(i,7) + 3*delta_r3[3]*d5_array_(i,10) + 6*delta_r3[4]*d5_array_(i,11) + 3*delta_r3[5]*d5_array_(i,12) + delta_r3[6]*d5_array_(i,15) + 3*delta_r3[7]*d5_array_(i,16) + 3*delta_r3[8]*d5_array_(i,17) + delta_r3[9]*d5_array_(i,18);
  third_term[4] = delta_r3[0]*d5_array_(i,4) + 3*delta_r3[1]*d5_array_(i,7) + 3*delta_r3[2]*d5_array_(i,8) + 3*delta_r3[3]*d5_array_(i,11) + 6*delta_r3[4]*d5_array_(i,12) + 3*delta_r3[5]*d5_array_(i,13) + delta_r3[6]*d5_array_(i,16) + 3*delta_r3[7]*d5_array_(i,17) + 3*delta_r3[8]*d5_array_(i,18) + delta_r3[9]*d5_array_(i,19);
  third_term[5] = delta_r3[0]*d5_array_(i,5) + 3*delta_r3[1]*d5_array_(i,8) + 3*delta_r3[2]*d5_array_(i,9) + 3*delta_r3[3]*d5_array_(i,12) + 6*delta_r3[4]*d5_array_(i,13) + 3*delta_r3[5]*d5_array_(i,14) + delta_r3[6]*d5_array_(i,17) + 3*delta_r3[7]*d5_array_(i,18) + 3*delta_r3[8]*d5_array_(i,19) + delta_r3[9]*d5_array_(i,20);

  // CkPrintf("d2 ewald zeroth_term: %f, %f, %f, %f, %f, %f\n", 
  //         d2_array_(i, 0), d2_array_(i, 1), d2_array_(i, 2), d2_array_(i, 3), d2_array_(i, 4), d2_array_(i, 5));
  // CkPrintf("d2 ewald first_term: %f, %f, %f, %f, %f, %f\n", 
  //         first_term[0], first_term[1], first_term[2], first_term[3], first_term[4], first_term[5]);
  // CkPrintf("d2 ewald second_term: %f, %f, %f, %f, %f, %f\n", 
  //         second_term[0], second_term[1], second_term[2], second_term[3], second_term[4], second_term[5]);
  // CkPrintf("d2 ewald third_term: %f, %f, %f, %f, %f, %f\n", 
  //         third_term[0], third_term[1], third_term[2], third_term[3], third_term[4], third_term[5]);

  std::array<double, 6> sum; 
  sum[0] = d2_array_(i, 0) + first_term[0] + 0.5*second_term[0] + 1.0/6.0*third_term[0];
  sum[1] = d2_array_(i, 1) + first_term[1] + 0.5*second_term[1] + 1.0/6.0*third_term[1];
  sum[2] = d2_array_(i, 2) + first_term[2] + 0.5*second_term[2] + 1.0/6.0*third_term[2];
  sum[3] = d2_array_(i, 3) + first_term[3] + 0.5*second_term[3] + 1.0/6.0*third_term[3];
  sum[4] = d2_array_(i, 4) + first_term[4] + 0.5*second_term[4] + 1.0/6.0*third_term[4];
  sum[5] = d2_array_(i, 5) + first_term[5] + 0.5*second_term[5] + 1.0/6.0*third_term[5];


  return sum;
  
}

std::array<double, 10> EnzoEwald::interp_d3(double x, double y, double z) throw()
{

  int i;
  double interp_x, interp_y, interp_z;
  find_nearest_interp_point(x, y, z, &interp_x, &interp_y, &interp_z, &i);

  std::array<double, 3> delta_r = {x - interp_x, y - interp_y, z - interp_z};
  std::array<double, 6> delta_r2 = EnzoMethodMultipole::outer_11_(delta_r, delta_r);
  std::array<double, 10> delta_r3 = EnzoMethodMultipole::outer_12_(delta_r, delta_r2);

  // std::array<double, 10> zeroth_term = d3_array_[i];
  // std::array<double, 10> first_term = EnzoMethodMultipole::dot_14_(delta_r, d4_array_[i]);
  // std::array<double, 10> second_term = EnzoMethodMultipole::dot_scalar_3_(0.5, EnzoMethodMultipole::dot_25_(delta_r2, d5_array_[i]));
  // std::array<double, 10> third_term = EnzoMethodMultipole::dot_scalar_3_(1.0/6.0, EnzoMethodMultipole::dot_36_(delta_r3, d6_array_[i]));

  // std::array<double, 10> zero_plus_one = EnzoMethodMultipole::add_33_(zeroth_term, first_term);
  // std::array<double, 10> two_plus_three = EnzoMethodMultipole::add_33_(second_term, third_term);

  std::array<double, 10> first_term;
  std::array<double, 10> second_term;
  std::array<double, 10> third_term;

  first_term[0] = delta_r[0]*d4_array_(i,0) + delta_r[1]*d4_array_(i,1) + delta_r[2]*d4_array_(i,2);
  first_term[1] = delta_r[0]*d4_array_(i,1) + delta_r[1]*d4_array_(i,3) + delta_r[2]*d4_array_(i,4);
  first_term[2] = delta_r[0]*d4_array_(i,2) + delta_r[1]*d4_array_(i,4) + delta_r[2]*d4_array_(i,5);
  first_term[3] = delta_r[0]*d4_array_(i,3) + delta_r[1]*d4_array_(i,6) + delta_r[2]*d4_array_(i,7);
  first_term[4] = delta_r[0]*d4_array_(i,4) + delta_r[1]*d4_array_(i,7) + delta_r[2]*d4_array_(i,8);
  first_term[5] = delta_r[0]*d4_array_(i,5) + delta_r[1]*d4_array_(i,8) + delta_r[2]*d4_array_(i,9);
  first_term[6] = delta_r[0]*d4_array_(i,6) + delta_r[1]*d4_array_(i,10) + delta_r[2]*d4_array_(i,11);
  first_term[7] = delta_r[0]*d4_array_(i,7) + delta_r[1]*d4_array_(i,11) + delta_r[2]*d4_array_(i,12);
  first_term[8] = delta_r[0]*d4_array_(i,8) + delta_r[1]*d4_array_(i,12) + delta_r[2]*d4_array_(i,13);
  first_term[9] = delta_r[0]*d4_array_(i,9) + delta_r[1]*d4_array_(i,13) + delta_r[2]*d4_array_(i,14);

  second_term[0] = delta_r2[0]*d5_array_(i,0) + 2*delta_r2[1]*d5_array_(i,1) + 2*delta_r2[2]*d5_array_(i,2) + delta_r2[3]*d5_array_(i,3) + 2*delta_r2[4]*d5_array_(i,4) + delta_r2[5]*d5_array_(i,5);
  second_term[1] = delta_r2[0]*d5_array_(i,1) + 2*delta_r2[1]*d5_array_(i,3) + 2*delta_r2[2]*d5_array_(i,4) + delta_r2[3]*d5_array_(i,6) + 2*delta_r2[4]*d5_array_(i,7) + delta_r2[5]*d5_array_(i,8);
  second_term[2] = delta_r2[0]*d5_array_(i,2) + 2*delta_r2[1]*d5_array_(i,4) + 2*delta_r2[2]*d5_array_(i,5) + delta_r2[3]*d5_array_(i,7) + 2*delta_r2[4]*d5_array_(i,8) + delta_r2[5]*d5_array_(i,9);
  second_term[3] = delta_r2[0]*d5_array_(i,3) + 2*delta_r2[1]*d5_array_(i,6) + 2*delta_r2[2]*d5_array_(i,7) + delta_r2[3]*d5_array_(i,10) + 2*delta_r2[4]*d5_array_(i,11) + delta_r2[5]*d5_array_(i,12);
  second_term[4] = delta_r2[0]*d5_array_(i,4) + 2*delta_r2[1]*d5_array_(i,7) + 2*delta_r2[2]*d5_array_(i,8) + delta_r2[3]*d5_array_(i,11) + 2*delta_r2[4]*d5_array_(i,12) + delta_r2[5]*d5_array_(i,13);
  second_term[5] = delta_r2[0]*d5_array_(i,5) + 2*delta_r2[1]*d5_array_(i,8) + 2*delta_r2[2]*d5_array_(i,9) + delta_r2[3]*d5_array_(i,12) + 2*delta_r2[4]*d5_array_(i,13) + delta_r2[5]*d5_array_(i,14);
  second_term[6] = delta_r2[0]*d5_array_(i,6) + 2*delta_r2[1]*d5_array_(i,10) + 2*delta_r2[2]*d5_array_(i,11) + delta_r2[3]*d5_array_(i,15) + 2*delta_r2[4]*d5_array_(i,16) + delta_r2[5]*d5_array_(i,17);
  second_term[7] = delta_r2[0]*d5_array_(i,7) + 2*delta_r2[1]*d5_array_(i,11) + 2*delta_r2[2]*d5_array_(i,12) + delta_r2[3]*d5_array_(i,16) + 2*delta_r2[4]*d5_array_(i,17) + delta_r2[5]*d5_array_(i,18);
  second_term[8] = delta_r2[0]*d5_array_(i,8) + 2*delta_r2[1]*d5_array_(i,12) + 2*delta_r2[2]*d5_array_(i,13) + delta_r2[3]*d5_array_(i,17) + 2*delta_r2[4]*d5_array_(i,18) + delta_r2[5]*d5_array_(i,19);
  second_term[9] = delta_r2[0]*d5_array_(i,9) + 2*delta_r2[1]*d5_array_(i,13) + 2*delta_r2[2]*d5_array_(i,14) + delta_r2[3]*d5_array_(i,18) + 2*delta_r2[4]*d5_array_(i,19) + delta_r2[5]*d5_array_(i,20);

  third_term[0] = delta_r3[0]*d6_array_(i,0) + 3*delta_r3[1]*d6_array_(i,1) + 3*delta_r3[2]*d6_array_(i,2) + 3*delta_r3[3]*d6_array_(i,3) + 6*delta_r3[4]*d6_array_(i,4) + 3*delta_r3[5]*d6_array_(i,5) + delta_r3[6]*d6_array_(i,6) + 3*delta_r3[7]*d6_array_(i,7) + 3*delta_r3[8]*d6_array_(i,8) + delta_r3[9]*d6_array_(i,9);
  third_term[1] = delta_r3[0]*d6_array_(i,1) + 3*delta_r3[1]*d6_array_(i,3) + 3*delta_r3[2]*d6_array_(i,4) + 3*delta_r3[3]*d6_array_(i,6) + 6*delta_r3[4]*d6_array_(i,7) + 3*delta_r3[5]*d6_array_(i,8) + delta_r3[6]*d6_array_(i,10) + 3*delta_r3[7]*d6_array_(i,11) + 3*delta_r3[8]*d6_array_(i,12) + delta_r3[9]*d6_array_(i,13);
  third_term[2] = delta_r3[0]*d6_array_(i,2) + 3*delta_r3[1]*d6_array_(i,4) + 3*delta_r3[2]*d6_array_(i,5) + 3*delta_r3[3]*d6_array_(i,7) + 6*delta_r3[4]*d6_array_(i,8) + 3*delta_r3[5]*d6_array_(i,9) + delta_r3[6]*d6_array_(i,11) + 3*delta_r3[7]*d6_array_(i,12) + 3*delta_r3[8]*d6_array_(i,13) + delta_r3[9]*d6_array_(i,14);
  third_term[3] = delta_r3[0]*d6_array_(i,3) + 3*delta_r3[1]*d6_array_(i,6) + 3*delta_r3[2]*d6_array_(i,7) + 3*delta_r3[3]*d6_array_(i,10) + 6*delta_r3[4]*d6_array_(i,11) + 3*delta_r3[5]*d6_array_(i,12) + delta_r3[6]*d6_array_(i,15) + 3*delta_r3[7]*d6_array_(i,16) + 3*delta_r3[8]*d6_array_(i,17) + delta_r3[9]*d6_array_(i,18);
  third_term[4] = delta_r3[0]*d6_array_(i,4) + 3*delta_r3[1]*d6_array_(i,7) + 3*delta_r3[2]*d6_array_(i,8) + 3*delta_r3[3]*d6_array_(i,11) + 6*delta_r3[4]*d6_array_(i,12) + 3*delta_r3[5]*d6_array_(i,13) + delta_r3[6]*d6_array_(i,16) + 3*delta_r3[7]*d6_array_(i,17) + 3*delta_r3[8]*d6_array_(i,18) + delta_r3[9]*d6_array_(i,19);
  third_term[5] = delta_r3[0]*d6_array_(i,5) + 3*delta_r3[1]*d6_array_(i,8) + 3*delta_r3[2]*d6_array_(i,9) + 3*delta_r3[3]*d6_array_(i,12) + 6*delta_r3[4]*d6_array_(i,13) + 3*delta_r3[5]*d6_array_(i,14) + delta_r3[6]*d6_array_(i,17) + 3*delta_r3[7]*d6_array_(i,18) + 3*delta_r3[8]*d6_array_(i,19) + delta_r3[9]*d6_array_(i,20);
  third_term[6] = delta_r3[0]*d6_array_(i,6) + 3*delta_r3[1]*d6_array_(i,10) + 3*delta_r3[2]*d6_array_(i,11) + 3*delta_r3[3]*d6_array_(i,15) + 6*delta_r3[4]*d6_array_(i,16) + 3*delta_r3[5]*d6_array_(i,17) + delta_r3[6]*d6_array_(i,21) + 3*delta_r3[7]*d6_array_(i,22) + 3*delta_r3[8]*d6_array_(i,23) + delta_r3[9]*d6_array_(i,24);
  third_term[7] = delta_r3[0]*d6_array_(i,7) + 3*delta_r3[1]*d6_array_(i,11) + 3*delta_r3[2]*d6_array_(i,12) + 3*delta_r3[3]*d6_array_(i,16) + 6*delta_r3[4]*d6_array_(i,17) + 3*delta_r3[5]*d6_array_(i,18) + delta_r3[6]*d6_array_(i,22) + 3*delta_r3[7]*d6_array_(i,23) + 3*delta_r3[8]*d6_array_(i,24) + delta_r3[9]*d6_array_(i,25);
  third_term[8] = delta_r3[0]*d6_array_(i,8) + 3*delta_r3[1]*d6_array_(i,12) + 3*delta_r3[2]*d6_array_(i,13) + 3*delta_r3[3]*d6_array_(i,17) + 6*delta_r3[4]*d6_array_(i,18) + 3*delta_r3[5]*d6_array_(i,19) + delta_r3[6]*d6_array_(i,23) + 3*delta_r3[7]*d6_array_(i,24) + 3*delta_r3[8]*d6_array_(i,25) + delta_r3[9]*d6_array_(i,26);
  third_term[9] = delta_r3[0]*d6_array_(i,9) + 3*delta_r3[1]*d6_array_(i,13) + 3*delta_r3[2]*d6_array_(i,14) + 3*delta_r3[3]*d6_array_(i,18) + 6*delta_r3[4]*d6_array_(i,19) + 3*delta_r3[5]*d6_array_(i,20) + delta_r3[6]*d6_array_(i,24) + 3*delta_r3[7]*d6_array_(i,25) + 3*delta_r3[8]*d6_array_(i,26) + delta_r3[9]*d6_array_(i,27);

  // CkPrintf("d3 ewald zeroth_term: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
  //         d3_array_(i, 0), d3_array_(i, 1), d3_array_(i, 2), d3_array_(i, 3), d3_array_(i, 4), d3_array_(i, 5),
  //         d3_array_(i, 6), d3_array_(i, 7), d3_array_(i, 8), d3_array_(i, 9));
  // CkPrintf("d3 ewald first_term: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
  //         first_term[0], first_term[1], first_term[2], first_term[3], first_term[4], first_term[5],
  //         first_term[6], first_term[7], first_term[8], first_term[9]);
  // CkPrintf("d3 ewald second_term: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
  //         second_term[0], second_term[1], second_term[2], second_term[3], second_term[4], second_term[5],
  //         second_term[6], second_term[7], second_term[8], second_term[9]);
  // CkPrintf("d3 ewald third_term: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
  //         third_term[0], third_term[1], third_term[2], third_term[3], third_term[4], third_term[5],
  //         third_term[6], third_term[7], third_term[8], third_term[9]);

  std::array<double, 10> sum; 
  sum[0] = d3_array_(i, 0) + first_term[0] + 0.5*second_term[0] + 1.0/6.0*third_term[0];
  sum[1] = d3_array_(i, 1) + first_term[1] + 0.5*second_term[1] + 1.0/6.0*third_term[1];
  sum[2] = d3_array_(i, 2) + first_term[2] + 0.5*second_term[2] + 1.0/6.0*third_term[2];
  sum[3] = d3_array_(i, 3) + first_term[3] + 0.5*second_term[3] + 1.0/6.0*third_term[3];
  sum[4] = d3_array_(i, 4) + first_term[4] + 0.5*second_term[4] + 1.0/6.0*third_term[4];
  sum[5] = d3_array_(i, 5) + first_term[5] + 0.5*second_term[5] + 1.0/6.0*third_term[5];
  sum[6] = d3_array_(i, 6) + first_term[6] + 0.5*second_term[6] + 1.0/6.0*third_term[6];
  sum[7] = d3_array_(i, 7) + first_term[7] + 0.5*second_term[7] + 1.0/6.0*third_term[7];
  sum[8] = d3_array_(i, 8) + first_term[8] + 0.5*second_term[8] + 1.0/6.0*third_term[8];
  sum[9] = d3_array_(i, 9) + first_term[9] + 0.5*second_term[9] + 1.0/6.0*third_term[9];

  return sum;
  
}


/********************************************************************************************/
/*  Compute the terms of the Taylor series required for interpolating Ewald derivatives */


// compute the d0 term of the Ewald sum at coordinates (x, y, z) -- only necessary if we want the potential
// double EnzoEwald::d0(double x, double y, double z) throw()
// {
//   /* d0 = g0 */

//   double d0_counter = 0;

//   Hierarchy * hierarchy = enzo::simulation()->hierarchy();
//   double lox, loy, loz; 
//   double hix, hiy, hiz;
//   hierarchy->lower(&lox, &loy, &loz);
//   hierarchy->upper(&hix, &hiy, &hiz);

//   double Lx = hix - lox;
//   double Ly = hiy - loy;
//   double Lz = hiz - loz;
//   double box_vol = Lx * Ly * Lz;

//   double alpha = 2.0 / pow(box_vol, 1./3.); 
//   double alpha2 = alpha * alpha;

//   // sum in real space 
//   for (int nz = -5; nz <= 5; nz++) {
//     for (int ny = -5; ny <= 5; ny++) {
//       for (int nx = -5; nx <= 5; nx++) {

//         double rx = x + nx * Lx;
//         double ry = y + ny * Ly;
//         double rz = z + nz * Lz;
//         double r = sqrt(rx*rx + ry*ry + rz*rz);

//         double g0;

//         if (nx != 0 || ny != 0 || nz !=0) {
//           g0 = erfc(alpha * r) / r;
//         }
//         else {
          
//           // Taylor expand if r is close to 0
//           if ((alpha * r) < 0.5) {
//             g0 = -2.0 * alpha / sqrt(M_PI) * (1.0 - pow(alpha * r, 2) / 3.0 + pow(alpha * r, 4) / 10.0 - pow(alpha * r, 6) / 42.0 
//                                             + pow(alpha * r, 8) / 216.0 - pow(alpha * r, 10) / 1320.0);
//           }
//           else { // incorporate Newtonian 1/r term
//             g0 = -erf(alpha*r)/r; // 1/r - erfc(alpha*r) / r;
//           }
//         }

//         d0_counter -= g0;
//       }
//     }
//   }

//   // sum in Fourier space
//   for (int nz = -5; nz <= 5; nz++) {
//     for (int ny = -5; ny <= 5; ny++) {
//       for (int nx = -5; nx <= 5; nx++) {

//         if (nx != 0 || ny != 0 || nz != 0) {

//           double kx = 2.0 * M_PI * nx / Lx;
//           double ky = 2.0 * M_PI * ny / Ly;
//           double kz = 2.0 * M_PI * nz / Lz;
//           double k2 = kx*kx + ky*ky + kz*kz;
//           double kdotx = kx*x + ky*y + kz*z;

//           double k_exp = exp(-k2 / (4.0*alpha2)) / k2;

//           d0_counter -= 4.0*M_PI/box_vol * k_exp * cos(kdotx);
//         } 
//       }
//     }
//   }

//   // additional pi/(alpha^2 V) term, only for d0
//   d0_counter += M_PI/(alpha2 * box_vol);


//   return d0_counter;  
// }



// compute the d1 term of the Ewald sum at coordinates (x, y, z)
CelloView<double, 1> EnzoEwald::d1(double x, double y, double z) throw()
{
  /* (d1)_i = g1 r_i */

  CelloView<double, 1> d1_counter(3);

  Hierarchy * hierarchy = enzo::simulation()->hierarchy();
  double lox, loy, loz; 
  double hix, hiy, hiz;
  hierarchy->lower(&lox, &loy, &loz);
  hierarchy->upper(&hix, &hiy, &hiz);

  double Lx = hix - lox;
  double Ly = hiy - loy;
  double Lz = hiz - loz;
  double box_vol = Lx * Ly * Lz;

  double alpha = 2.0 / pow(box_vol, 1./3.);
  double alpha2 = alpha * alpha;

  int realnxmax = (int)(8.0 / (alpha * Lx) + 0.5);
  int realnymax = (int)(8.0 / (alpha * Ly) + 0.5);
  int realnzmax = (int)(8.0 / (alpha * Lz) + 0.5);

  // sum in real space
  for (int nz = -realnzmax; nz <= realnzmax; nz++) {
    for (int ny = -realnymax; ny <= realnymax; ny++) {
      for (int nx = -realnxmax; nx <= realnxmax; nx++) {

        double rx = x + nx * Lx;
        double ry = y + ny * Ly;
        double rz = z + nz * Lz;
        double r = sqrt(rx*rx + ry*ry + rz*rz);
        // if (r==0) CkPrintf("r = 0: %d, %d, %d\n", nx, ny, nz);
        double r2 = r*r;
        double r3 = r*r2;

        double g1 = 0;

        if (nx != 0 || ny != 0 || nz != 0) {
          if (r != 0) {
            g1 = (-2.0 * alpha * r - sqrt(M_PI) * exp(alpha2*r2) * erfc(alpha*r)) * exp(-1.0 * alpha2 * r2) / (sqrt(M_PI) * r3);
          }
        }
        else {
          
          // Taylor expand if r is close to 0
          if ((alpha * r) < 0.5) {
            g1 = -4.0 * pow(alpha, 3) / sqrt(M_PI) *
                       (-1.0 / 3.0 + pow(alpha * r, 2) / 5.0 - pow(alpha * r, 4) / 14.0 + pow(alpha * r, 6) / 54.0 -
                        pow(alpha * r, 8) / 264.0 + pow(alpha * r, 10) / 1560.0);

            //if (abs(x + 0.007937) <= 1e-5 && abs(y + 0.007937) <= 1e-5 && abs(z + 0) < 1e-2) {
	    //CkPrintf("D1 Taylor: g1, rx, ry, rz: %f, %f, %f, %f\n", g1, rx, ry, rz);
            //}
          }
          else { // incorporate Newtonian 1/r term
            g1 = 1.0/r3 + ((-2.0 * alpha * r - sqrt(M_PI) * exp(alpha2*r2) * erfc(alpha*r)) 
                  * exp(-1.0 * alpha2 * r2) / (sqrt(M_PI) * r3));
          }
        }

        d1_counter(0) -= g1 * rx;
        d1_counter(1) -= g1 * ry;
        d1_counter(2) -= g1 * rz;

      }
    }
  }
  
  int knxmax = (int)(2.0 * alpha * Lx + 0.5);
  int knymax = (int)(2.0 * alpha * Ly + 0.5);
  int knzmax = (int)(2.0 * alpha * Lz + 0.5);

  // sum in Fourier space
  for (int nz = -knzmax; nz <= knzmax; nz++) {
    for (int ny = -knymax; ny <= knymax; ny++) {
      for (int nx = -knxmax; nx <= knxmax; nx++) {

        if (nx != 0 || ny != 0 || nz != 0) {

          double kx = 2.0 * M_PI * nx / Lx;
          double ky = 2.0 * M_PI * ny / Ly;
          double kz = 2.0 * M_PI * nz / Lz;
          double k2 = kx*kx + ky*ky + kz*kz;
          double kdotx = kx*x + ky*y + kz*z;

          double k_exp = 4.0*M_PI/box_vol * (exp(-k2 / (4.0*alpha2)) / k2) * sin(kdotx); 

          d1_counter(0) += k_exp * kx;
          d1_counter(1) += k_exp * ky;
          d1_counter(2) += k_exp * kz;
        } 
      }
    }
  }


  return d1_counter;

}

// compute the d2 term of the Ewald sum at coordinates (x, y, z)
CelloView<double, 1> EnzoEwald::d2(double x, double y, double z) throw()
{
  /* (d2)_ij = g1 delta_ij + g2 r_i r_j */

  CelloView<double, 1> d2_counter(6);

  Hierarchy * hierarchy = enzo::simulation()->hierarchy();
  double lox, loy, loz; 
  double hix, hiy, hiz;
  hierarchy->lower(&lox, &loy, &loz);
  hierarchy->upper(&hix, &hiy, &hiz);

  double Lx = hix - lox;
  double Ly = hiy - loy;
  double Lz = hiz - loz;
  double box_vol = Lx * Ly * Lz;

  double alpha = 2.0 / pow(box_vol, 1./3.);
  double alpha2 = alpha * alpha;
  double alpha3 = alpha2 * alpha;

  int realnxmax = (int)(8.0 / (alpha * Lx) + 0.5);
  int realnymax = (int)(8.0 / (alpha * Ly) + 0.5);
  int realnzmax = (int)(8.0 / (alpha * Lz) + 0.5);

  // sum in real space
  for (int nz = -realnzmax; nz <= realnzmax; nz++) {
    for (int ny = -realnymax; ny <= realnymax; ny++) {
      for (int nx = -realnxmax; nx <= realnxmax; nx++) {

        double rx = x + nx * Lx;
        double ry = y + ny * Ly;
        double rz = z + nz * Lz;
        double r = sqrt(rx*rx + ry*ry + rz*rz);
        double r2 = r*r;
        double r3 = r*r2;
        double r5 = r2*r3;

        double g1 = 0;
	      double g2 = 0;

        if (nx != 0 || ny != 0 || nz != 0) {
          if (r != 0) {
            // g1 = (-2.0 * alpha * r - sqrt(M_PI) * exp(alpha2*r2) * erfc(alpha*r)) * exp(-1.0 * alpha2 * r2) / (sqrt(M_PI) * r3);

            // g2 = (4.0 * sqrt(M_PI) * alpha3 * r3 + 6.0 * sqrt(M_PI) * alpha * r + 3.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r))
            //   * exp(-1.0 * alpha2 * r2) / (M_PI * r5);


            g1 = -(erfc(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) / r3;

            g2 = (3.0 * erfc(alpha * r) + (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) / r5;

                  
              // if (abs(x + 0.007937) <= 1e-5 && abs(y + 0.007937) <= 1e-5 && abs(z + 0) < 1e-2) {
              //     CkPrintf("Away: g1, g2, rx, ry, rz: %f, %f, %f, %f, %f\n", g1, g2, rx, ry, rz);
              //   }
          }
        }
        else {
          
          // Taylor expand if r is close to 0
          if ((alpha * r) < 0.5) {
	    
            g1 = -4.0 * pow(alpha, 3) / sqrt(M_PI) *
                       (-1.0 / 3.0 + pow(alpha * r, 2) / 5.0 - pow(alpha * r, 4) / 14.0 + pow(alpha * r, 6) / 54.0 -
                        pow(alpha * r, 8) / 264.0 + pow(alpha * r, 10) / 1560.0);

            g2 = -8.0 * pow(alpha, 5) / sqrt(M_PI) *
                       (1.0 / 5.0 - pow(alpha * r, 2) / 7.0 + pow(alpha * r, 4) / 18.0 - pow(alpha * r, 6) / 66.0 +
                        pow(alpha * r, 8) / 312.0 - pow(alpha * r, 10) / 1800.0);

	          // if (abs(x + 0.007937) <= 1e-5 && abs(y + 0.007937) <= 1e-5 && abs(z + 0) < 1e-2) {
            //  CkPrintf("D2 Taylor: g1, g2, rx, ry, rz: %f, %f, %f, %f, %f\n", g1, g2, rx, ry, rz);
            // }

          }
          else { // incorporate Newtonian 1/r term
            if (r != 0) {


              // g1 = 1.0/r3 + ((-2.0 * alpha * r - sqrt(M_PI) * exp(alpha2*r2) * erfc(alpha*r)) 
              //       * exp(-1.0 * alpha2 * r2) / (sqrt(M_PI) * r3));

              // g2 = -3.0/r5 + (4.0 * sqrt(M_PI) * alpha3 * r3 + 6.0 * sqrt(M_PI) * alpha * r + 3.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r))
              //     * exp(-1.0 * alpha2 * r2) / (M_PI * r5);

              g1 = -(-erf(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) / r3;

              g2 = -(3.0 * erf(alpha * r) - (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) / r5;
            }
            
            // if (abs(x + 0.007937) <= 1e-5 && abs(y + 0.007937) <= 1e-5 && abs(z + 0) < 1e-2) {
              // CkPrintf("Primary domain: g1, g2, rx, ry, rz: %f, %f, %f, %f, %f\n", g1, g2, rx, ry, rz);
	    // }
          }
        }


        d2_counter(0) -= g2*rx*rx + g1;
        d2_counter(1) -= g2*rx*ry;
        d2_counter(2) -= g2*rx*rz; 
        d2_counter(3) -= g2*ry*ry + g1;
        d2_counter(4) -= g2*ry*rz; 
        d2_counter(5) -= g2*rz*rz + g1;

      }
    }
  }

  // if (x == 0 && y == 0 && z == 0) CkPrintf("d2 counter 1: %f\n", d2_counter(0));

  //if (abs(x + 0.007937) <= 1e-5 && abs(y + 0.007937) <= 1e-5 && abs(z + 0) < 1e-2) {
  //  CkPrintf("D2_counter: %f, %f, %f, %f, %f, %f\n", d2_counter(0), d2_counter(1), d2_counter(2), d2_counter(3), d2_counter(4), d2_counter(5));
  // }
  
  int knxmax = (int)(2.0 * alpha * Lx + 0.5);
  int knymax = (int)(2.0 * alpha * Ly + 0.5);
  int knzmax = (int)(2.0 * alpha * Lz + 0.5);

  // sum in Fourier space
  for (int nz = -knzmax; nz <= knzmax; nz++) {
    for (int ny = -knymax; ny <= knymax; ny++) {
      for (int nx = -knxmax; nx <= knxmax; nx++) {

        if (nx != 0 || ny != 0 || nz != 0) {

          double kx = 2.0 * M_PI * nx / Lx;
          double ky = 2.0 * M_PI * ny / Ly;
          double kz = 2.0 * M_PI * nz / Lz;
          double k2 = kx*kx + ky*ky + kz*kz;
          double kdotx = kx*x + ky*y + kz*z;

          double k_exp = 4.0*M_PI/box_vol * (exp(-1.0*k2 / (4.0*alpha2)) / k2) * cos(kdotx);

          d2_counter(0) += k_exp * kx*kx;
          d2_counter(1) += k_exp * kx*ky;
          d2_counter(2) += k_exp * kx*kz;
          d2_counter(3) += k_exp * ky*ky;
          d2_counter(4) += k_exp * ky*kz;
          d2_counter(5) += k_exp * kz*kz;
        } 
      }
    }
  }

  // if (x == 0 && y == 0 && z == 0) CkPrintf("d2 counter 2: %f\n", d2_counter(0));

  return d2_counter;
  
}

CelloView<double, 1> EnzoEwald::d2_gadget(double x, double y, double z) throw()
{


  CelloView<double, 1> d2_counter(6);

  Hierarchy * hierarchy = enzo::simulation()->hierarchy();
  double lox, loy, loz; 
  double hix, hiy, hiz;
  hierarchy->lower(&lox, &loy, &loz);
  hierarchy->upper(&hix, &hiy, &hiz);

  double LONG_X = 1.0/(hix - lox);
  double LONG_Y = 1.0/(hiy - loy);
  double LONG_Z = 1.0/(hiz - loz);


  double leff   = pow((1.0 / LONG_X) * (1.0 / LONG_Y) * (1.0 / LONG_Z), 1.0 / 3);
  double alpha  = 2.0 / leff;
  double alpha2 = alpha * alpha;

  int qxmax = (int)(8.0 * LONG_X / alpha + 0.5);
  int qymax = (int)(8.0 * LONG_Y / alpha + 0.5);
  int qzmax = (int)(8.0 * LONG_Z / alpha + 0.5);

  int nxmax = (int)(2.0 * alpha / LONG_X + 0.5);
  int nymax = (int)(2.0 * alpha / LONG_Y + 0.5);
  int nzmax = (int)(2.0 * alpha / LONG_Z + 0.5);

 

  for(int nx = -qxmax; nx <= qxmax; nx++)
    for(int ny = -qymax; ny <= qymax; ny++)
      for(int nz = -qzmax; nz <= qzmax; nz++)
        {
          double dx = x - nx * (1.0 / LONG_X);
          double dy = y - ny * (1.0 / LONG_Y);
          double dz = z - nz * (1.0 / LONG_Z);

          double r2 = dx * dx + dy * dy + dz * dz;
          double r  = sqrt(r2);

          double rinv  = (r > 0) ? 1.0 / r : 0.0;
          double r2inv = rinv * rinv;
          double r3inv = r2inv * rinv;
          double r5inv = r3inv * r2inv;

          double g1, g2;

          if(nx != 0 || ny != 0 || nz != 0)
            {
              g1 = (erfc(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) * r3inv;

              g2 = -(3.0 * erfc(alpha * r) + (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;
            }
          else
            {
              /* we add the 1/r term here to the (0|0|0) entry, followed by differentiation, and the limit r->0 to obtain accurate
               * results at the origin
               */

              /* Note, for small r:
               *
               *   [1/- erfc(a r)]/r  =  2 a/sqrt(pi) * [ 1 - (a r)^2/3 + (a r)^4 / 10 - (a r)^6 / 42 + (a r)^8 / 216 - ...]
               *
               *   Hence for r = 0:
               *
               *   g0 =  2     * alpha   / sqrt(pi)
               *   g1 = -4/3   * alpha^3 / sqrt(pi)
               *   g2 =  8/5   * alpha^5 / sqrt(pi)
               *   g3 = -16/7  * alpha^7 / sqrt(pi)
               *   g4 =  32/9  * alpha^9 / sqrt(pi)
               *   g5 = -64/11 * alpha^11/ sqrt(pi)
               */

              if((alpha * r) < 0.5)
                {
                  g1 = 4.0 * pow(alpha, 3) / sqrt(M_PI) *
                       (-1.0 / 3.0 + pow(alpha * r, 2) / 5.0 - pow(alpha * r, 4) / 14.0 + pow(alpha * r, 6) / 54.0 -
                        pow(alpha * r, 8) / 264.0 + pow(alpha * r, 10) / 1560.0);

                  g2 = 8.0 * pow(alpha, 5) / sqrt(M_PI) *
                       (1.0 / 5.0 - pow(alpha * r, 2) / 7.0 + pow(alpha * r, 4) / 18.0 - pow(alpha * r, 6) / 66.0 +
                        pow(alpha * r, 8) / 312.0 - pow(alpha * r, 10) / 1800.0);
                }
              else
                {
                  g1 = (-erf(alpha * r) + 2.0 * alpha * r / sqrt(M_PI) * exp(-alpha2 * r2)) * r3inv;

                  g2 = (3.0 * erf(alpha * r) - (6.0 * alpha * r + 4.0 * pow(alpha * r, 3)) / sqrt(M_PI) * exp(-alpha2 * r2)) * r5inv;
                }
            }

	  d2_counter(0) += g2*dx*dx + g1;
	  d2_counter(1) += g2*dx*dy;
	  d2_counter(2) += g2*dx*dz; 
	  d2_counter(3) += g2*dy*dy + g1;
	  d2_counter(4) += g2*dy*dz; 
	  d2_counter(5) += g2*dz*dz + g1;

          //D2 += g2 * (dxyz % dxyz);
          //D2[qXX] += g1;
          //D2[qYY] += g1;
          //D2[qZZ] += g1;
        }

  for(int nx = -nxmax; nx <= nxmax; nx++)
    for(int ny = -nymax; ny <= nymax; ny++)
      for(int nz = -nzmax; nz <= nzmax; nz++)
        {
          if(nx != 0 || ny != 0 || nz != 0)
            {
              double kx = (2.0 * M_PI * LONG_X) * nx;
              double ky = (2.0 * M_PI * LONG_Y) * ny;
              double kz = (2.0 * M_PI * LONG_Z) * nz;
              double k2 = kx * kx + ky * ky + kz * kz;

              double kdotx = (x * kx + y * ky + z * kz);
              double val   = 4.0 * M_PI * (LONG_X * LONG_Y * LONG_Z) / k2 * exp(-k2 / (4.0 * alpha2)) * cos(kdotx);

              
	            d2_counter(0) += val * kx*kx;
              d2_counter(1) += val * kx*ky;
              d2_counter(2) += val * kx*kz;
              d2_counter(3) += val * ky*ky;
              d2_counter(4) += val * ky*kz;
              d2_counter(5) += val * kz*kz;
	      
              // D2 += (val * kxyz) % kxyz;
            }
        }

  // if (x == 0 && y == 0 && z == 0) CkPrintf("d2 counter gadget: %f\n", d2_counter(0));

  return d2_counter;
}



// compute the d3 term of the Ewald sum at coordinates (x, y, z)
CelloView<double, 1> EnzoEwald::d3(double x, double y, double z) throw()
{
  /* (d3)_ijk = g2 * (delta_ij * r_k + delta_jk * r_i + delta_ik * r_j) + g_3 * r_i * r_j * r_k */

  CelloView<double, 1> d3_counter(10);

  Hierarchy * hierarchy = enzo::simulation()->hierarchy();
  double lox, loy, loz; 
  double hix, hiy, hiz;
  hierarchy->lower(&lox, &loy, &loz);
  hierarchy->upper(&hix, &hiy, &hiz);

  double Lx = hix - lox;
  double Ly = hiy - loy;
  double Lz = hiz - loz;
  double box_vol = Lx * Ly * Lz;

  double alpha = 2.0 / pow(box_vol, 1./3.);
  double alpha2 = alpha * alpha;
  double alpha3 = alpha2 * alpha;
  double alpha5 = alpha2 * alpha3;

  int realnxmax = (int)(8.0 / (alpha * Lx) + 0.5);
  int realnymax = (int)(8.0 / (alpha * Ly) + 0.5);
  int realnzmax = (int)(8.0 / (alpha * Lz) + 0.5);

  // sum in real space
  for (int nz = -realnzmax; nz <= realnzmax; nz++) {
    for (int ny = -realnymax; ny <= realnymax; ny++) {
      for (int nx = -realnxmax; nx <= realnxmax; nx++) {

        double rx = x + nx * Lx;
        double ry = y + ny * Ly;
        double rz = z + nz * Lz;
        double r = sqrt(rx*rx + ry*ry + rz*rz);
        double r2 = r*r;
        double r3 = r2*r;
        double r5 = r2*r3;
        double r7 = r2*r5;

        double g2 = 0;
	      double g3 = 0;

        if (nx != 0 || ny != 0 || nz != 0) {
          if (r != 0) {
            // g2 = (4.0 * sqrt(M_PI) * alpha3 * r3 + 6.0 * sqrt(M_PI) * alpha * r + 3.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r))
            //   * exp(-1.0 * alpha2 * r2) / (M_PI * r5);

            // g3 = (-8.0 * sqrt(M_PI) * alpha5 * r5 - 20.0 * sqrt(M_PI) * alpha3 * r3 - 30.0 * sqrt(M_PI) * alpha * r
            //   - 15.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r7);

            g2 = (3.0 * erfc(alpha * r) + (6.0 * alpha * r + 4.0 * alpha3 * r3) / sqrt(M_PI) * exp(-alpha2 * r2)) / r5;

            g3 = -(15.0 * erfc(alpha * r) +
                  (30.0 * alpha * r + 20.0 * alpha3 * r3 + 8.0 * alpha5 * r5) / sqrt(M_PI) * exp(-alpha2 * r2)) / r7;

          }
        }
        else {
          
          // Taylor expand if r is close to 0
          if ((alpha * r) < 0.5) {
            g2 = -8.0 * pow(alpha, 5) / sqrt(M_PI) *
                       (1.0 / 5.0 - pow(alpha * r, 2) / 7.0 + pow(alpha * r, 4) / 18.0 - pow(alpha * r, 6) / 66.0 +
                        pow(alpha * r, 8) / 312.0 - pow(alpha * r, 10) / 1800.0);

            g3 = -16.0 * pow(alpha, 7) / sqrt(M_PI) *
                       (-1.0 / 7.0 + pow(alpha * r, 2) / 9.0 - pow(alpha * r, 4) / 22.0 + pow(alpha * r, 6) / 78.0 -
                        pow(alpha * r, 8) / 360.0 + pow(alpha * r, 10) / 2040.0);

            //if (abs(x + 0.007937) <= 1e-5 && abs(y + 0.007937) <= 1e-5 && abs(z + 0) < 1e-2) {
            //  CkPrintf("D3 Taylor: g2, g3, rx, ry, rz: %f, %f, %f, %f, %f\n", g2, g3, rx, ry, rz);
            //}

          }
          else { // incorporate Newtonian 1/r term
            // g2 = -3.0/r5 + (4.0 * sqrt(M_PI) * alpha3 * r3 + 6.0 * sqrt(M_PI) * alpha * r + 3.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r))
            //     * exp(-1.0 * alpha2 * r2) / (M_PI * r5);

            // g3 = 15.0/r7 + (-8.0 * sqrt(M_PI) * alpha5 * r5 - 20.0 * sqrt(M_PI) * alpha3 * r3 - 30.0 * sqrt(M_PI) * alpha * r
            //     - 15.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r7);

            g2 = -(3.0 * erf(alpha * r) - (6.0 * alpha * r + 4.0 * alpha3 * r3) / sqrt(M_PI) * exp(-alpha2 * r2)) / r5;

            g3 = -(-15.0 * erf(alpha * r) +
                  (30.0 * alpha * r + 20.0 * alpha3 * r3 + 8.0 * alpha5 * r5) / sqrt(M_PI) * exp(-alpha2 * r2)) / r7;

          }
        }

        d3_counter(0) -= g3*rx*rx*rx + 3*g2*rx;
        d3_counter(1) -= g3*rx*rx*ry + g2*ry;
        d3_counter(2) -= g3*rx*rx*rz + g2*rz;
        d3_counter(3) -= g3*rx*ry*ry + g2*rx;
        d3_counter(4) -= g3*rx*ry*rz;
        d3_counter(5) -= g3*rx*rz*rz + g2*rx;
        d3_counter(6) -= g3*ry*ry*ry + 3*g2*ry;
        d3_counter(7) -= g3*ry*ry*rz + g2*rz;
        d3_counter(8) -= g3*ry*rz*rz + g2*ry;
        d3_counter(9) -= g3*rz*rz*rz + 3*g2*rz;
      }
    }
  }
  
  int knxmax = (int)(2.0 * alpha * Lx + 0.5);
  int knymax = (int)(2.0 * alpha * Ly + 0.5);
  int knzmax = (int)(2.0 * alpha * Lz + 0.5);

  // sum in Fourier space
  for (int nz = -knzmax; nz <= knzmax; nz++) {
    for (int ny = -knymax; ny <= knymax; ny++) {
      for (int nx = -knxmax; nx <= knxmax; nx++) {

        if (nx != 0 || ny != 0 || nz != 0) {

          double kx = 2.0 * M_PI * nx / Lx;
          double ky = 2.0 * M_PI * ny / Ly;
          double kz = 2.0 * M_PI * nz / Lz;
          double k2 = kx*kx + ky*ky + kz*kz;
          double kdotx = kx*x + ky*y + kz*z;

          double k_exp = 4.0*M_PI/box_vol * (exp(-1.0*k2 / (4.0*alpha2)) / k2) * sin(kdotx); 

          d3_counter(0) -= k_exp * kx*kx*kx;
          d3_counter(1) -= k_exp * kx*kx*ky;
          d3_counter(2) -= k_exp * kx*kx*kz;
          d3_counter(3) -= k_exp * kx*ky*ky;
          d3_counter(4) -= k_exp * kx*ky*kz;
          d3_counter(5) -= k_exp * kx*kz*kz;
          d3_counter(6) -= k_exp * ky*ky*ky;
          d3_counter(7) -= k_exp * ky*ky*kz;
          d3_counter(8) -= k_exp * ky*kz*kz;
          d3_counter(9) -= k_exp * kz*kz*kz;
        } 
      }
    }
  }

  // if (x == 0 && y == 0 && z == 0) CkPrintf("d3 counter 2: %f\n", d3_counter(0));

  return d3_counter;
  
}

// compute the d4 term of the Ewald sum at coordinates (x, y, z)
CelloView<double, 1> EnzoEwald::d4(double x, double y, double z) throw()
{
  /* (d4)_ijkl = g2 * (delta_ij * delta_kl + delta_jk * delta_il + delta_ik * delta_jl) 
               + g3 * (delta_ij r_k r_l + delta_jk r_i r_l + delta_ik r_j r_l + delta_il r_j r_k + delta_jl r_i r_k + delta_kl r_i r_j)
               + g4 * r_i * r_j * r_k * r_l      */

  CelloView<double, 1> d4_counter(15);

  Hierarchy * hierarchy = enzo::simulation()->hierarchy();
  double lox, loy, loz; 
  double hix, hiy, hiz;
  hierarchy->lower(&lox, &loy, &loz);
  hierarchy->upper(&hix, &hiy, &hiz);

  double Lx = hix - lox;
  double Ly = hiy - loy;
  double Lz = hiz - loz;
  double box_vol = Lx * Ly * Lz;

  double alpha = 2.0 / pow(box_vol, 1./3.);
  double alpha2 = alpha * alpha;
  double alpha3 = alpha2 * alpha;
  double alpha5 = alpha2 * alpha3;
  double alpha7 = alpha2 * alpha5;

  int realnxmax = (int)(8.0 / (alpha * Lx) + 0.5);
  int realnymax = (int)(8.0 / (alpha * Ly) + 0.5);
  int realnzmax = (int)(8.0 / (alpha * Lz) + 0.5);

  // sum in real space
  for (int nz = -realnzmax; nz <= realnzmax; nz++) {
    for (int ny = -realnymax; ny <= realnymax; ny++) {
      for (int nx = -realnxmax; nx <= realnxmax; nx++) {

        double rx = x + nx * Lx;
        double ry = y + ny * Ly;
        double rz = z + nz * Lz;
        double r = sqrt(rx*rx + ry*ry + rz*rz);
        double r2 = r*r;
        double r3 = r2*r;
        double r5 = r2*r3;
        double r7 = r2*r5;
        double r9 = r2*r7;

        double g2 = 0;
        double g3 = 0;
        double g4 = 0;

        if (nx != 0 || ny != 0 || nz != 0) {
          if (r != 0) {
            g2 = (4.0 * sqrt(M_PI) * alpha3 * r3 + 6.0 * sqrt(M_PI) * alpha * r + 3.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r))
              * exp(-1.0 * alpha2 * r2) / (M_PI * r5);

            g3 = (-8.0 * sqrt(M_PI) * alpha5 * r5 - 20.0 * sqrt(M_PI) * alpha3 * r3 - 30.0 * sqrt(M_PI) * alpha * r
            - 15.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r7);

            g4 = (16.0 * sqrt(M_PI) * alpha7 * r7  + 56.0 * sqrt(M_PI) * alpha5 * r5 + 140.0 * sqrt(M_PI) * alpha3 * r3
            + 210.0 * sqrt(M_PI) * alpha * r + 105.0 * M_PI * exp(alpha2 * r2) * erfc(alpha * r)) 
              * exp(-1.0 * alpha2 * r2) / (M_PI * r9);
          }
        }
        else {
          
          // Taylor expand if r is close to 0
          if ((alpha * r) < 0.5) {
            g2 = -8.0 * pow(alpha, 5) / sqrt(M_PI) *
                       (1.0 / 5.0 - pow(alpha * r, 2) / 7.0 + pow(alpha * r, 4) / 18.0 - pow(alpha * r, 6) / 66.0 +
                        pow(alpha * r, 8) / 312.0 - pow(alpha * r, 10) / 1800.0);

            g3 = -16.0 * pow(alpha, 7) / sqrt(M_PI) *
                       (-1.0 / 7.0 + pow(alpha * r, 2) / 9.0 - pow(alpha * r, 4) / 22.0 + pow(alpha * r, 6) / 78.0 -
                        pow(alpha * r, 8) / 360.0 + pow(alpha * r, 10) / 2040.0);

            g4 = -32.0 * pow(alpha, 9) / sqrt(M_PI) *
                       (1.0 / 9.0 - pow(alpha * r, 2) / 11.0 + pow(alpha * r, 4) / 26.0 - pow(alpha * r, 6) / 90.0 +
                        pow(alpha * r, 8) / 408.0 - pow(alpha * r, 10) / 2280.0);

            //if (abs(x + 0.007937) <= 1e-5 && abs(y + 0.007937) <= 1e-5 && abs(z + 0) < 1e-2) {
            //  CkPrintf("D4 Taylor: g2, g3, g4, rx, ry, rz: %f, %f, %f, %f, %f, %f\n", g2, g3, g4, rx, ry, rz);
	    // }

          }
          else { // incorporate Newtonian 1/r term
            g2 = -3.0/r5 + (4.0 * sqrt(M_PI) * alpha3 * r3 + 6.0 * sqrt(M_PI) * alpha * r + 3.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r))
                * exp(-1.0 * alpha2 * r2) / (M_PI * r5);

            g3 = 15.0/r7 + (-8.0 * sqrt(M_PI) * alpha5 * r5 - 20.0 * sqrt(M_PI) * alpha3 * r3 - 30.0 * sqrt(M_PI) * alpha * r
                - 15.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r7);
            
            g4 = -105.0/r9 + (16.0 * sqrt(M_PI) * alpha7 * r7  + 56.0 * sqrt(M_PI) * alpha5 * r5 + 140.0 * sqrt(M_PI) * alpha3 * r3
                + 210.0 * sqrt(M_PI) * alpha * r + 105.0 * M_PI * exp(alpha2 * r2) * erfc(alpha * r)) 
                * exp(-1.0 * alpha2 * r2) / (M_PI * r9);
          }
        }

        d4_counter(0) -= g4*rx*rx*rx*rx + 6*g3*rx*rx + 3*g2;
        d4_counter(1) -= g4*rx*rx*rx*ry + 3*g3*rx*ry;
        d4_counter(2) -= g4*rx*rx*rx*rz + 3*g3*rx*rz;
        d4_counter(3) -= g4*rx*rx*ry*ry + g3*(rx*rx + ry*ry) + g2;
        d4_counter(4) -= g4*rx*rx*ry*rz + g3*ry*rz;
        d4_counter(5) -= g4*rx*rx*rz*rz + g3*(rx*rx + rz*rz) + g2;
        d4_counter(6) -= g4*rx*ry*ry*ry + 3*g3*rx*ry;
        d4_counter(7) -= g4*rx*ry*ry*rz + g3*rx*rz;
        d4_counter(8) -= g4*rx*ry*rz*rz + g3*rx*ry;
        d4_counter(9) -= g4*rx*rz*rz*rz + 3*g3*rx*rz;
        d4_counter(10) -= g4*ry*ry*ry*ry + 6*g3*ry*ry + 3*g2;
        d4_counter(11) -= g4*ry*ry*ry*rz + 3*g3*ry*rz;
        d4_counter(12) -= g4*ry*ry*rz*rz + g3*(ry*ry + rz*rz) + g2;
        d4_counter(13) -= g4*ry*rz*rz*rz + 3*g3*ry*rz;
        d4_counter(14) -= g4*rz*rz*rz*rz + 6*g3*rz*rz + 3*g2;

      }
    }
  }

  // if (x == 0 && y == 0 && z == 0) CkPrintf("d4 counter 1: %f\n", d4_counter(0));

  int knxmax = (int)(2.0 * alpha * Lx + 0.5);
  int knymax = (int)(2.0 * alpha * Ly + 0.5);
  int knzmax = (int)(2.0 * alpha * Lz + 0.5);
  
  // sum in Fourier space
  for (int nz = -knzmax; nz <= knzmax; nz++) {
    for (int ny = -knymax; ny <= knymax; ny++) {
      for (int nx = -knxmax; nx <= knxmax; nx++) {

        if (nx != 0 || ny != 0 || nz != 0) {

          double kx = 2.0 * M_PI * nx / Lx;
          double ky = 2.0 * M_PI * ny / Ly;
          double kz = 2.0 * M_PI * nz / Lz;
          double k2 = kx*kx + ky*ky + kz*kz;
          double kdotx = kx*x + ky*y + kz*z;

          double k_exp = 4.0*M_PI/box_vol * (exp(-1.0*k2 / (4.0*alpha2)) / k2) * cos(kdotx); 

          d4_counter(0) -= k_exp * kx*kx*kx*kx;
          d4_counter(1) -= k_exp * kx*kx*kx*ky;
          d4_counter(2) -= k_exp * kx*kx*kx*kz;
          d4_counter(3) -= k_exp * kx*kx*ky*ky;
          d4_counter(4) -= k_exp * kx*kx*ky*kz;
          d4_counter(5) -= k_exp * kx*kx*kz*kz;
          d4_counter(6) -= k_exp * kx*ky*ky*ky;
          d4_counter(7) -= k_exp * kx*ky*ky*kz;
          d4_counter(8) -= k_exp * kx*ky*kz*kz;
          d4_counter(9) -= k_exp * kx*kz*kz*kz;
          d4_counter(10) -= k_exp * ky*ky*ky*ky;
          d4_counter(11) -= k_exp * ky*ky*ky*kz;
          d4_counter(12) -= k_exp * ky*ky*kz*kz;
          d4_counter(13) -= k_exp * ky*kz*kz*kz;
          d4_counter(14) -= k_exp * kz*kz*kz*kz;

        } 
      }
    }
  }

  // if (x == 0 && y == 0 && z == 0) CkPrintf("d4 counter 2: %f\n", d4_counter(0));

  return d4_counter;
  
}

// compute the d5 term of the Ewald sum at coordinates (x, y, z)
CelloView<double, 1> EnzoEwald::d5(double x, double y, double z) throw()
{

  /* (d5)_ijklm = g3 [r_m (delta_ij delta_kl + delta_jk delta_il + delta_ik delta_jl) 
                      + r_l (delta_ij delta_km + delta_jk delta_im + delta_ik delta_jm) 
                      + r_k (delta_ij delta_lm + delta_jm delta_il + delta_im delta_jl)
                      + r_j (delta_im delta_kl + delta_km delta_il + delta_ik delta_lm)
                      + r_i (delta_jm delta_kl + delta_jk delta_lm + delta_km delta_jl)]
                  + g4 (delta_ij r_k r_l r_m + delta_jk r_i r_l r_m + delta_ik r_j r_l r_m
                      + delta_il r_j r_k r_m + delta_jl r_i r_k r_m + delta_kl r_i r_j r_m
                      + delta_im r_j r_k r_l + delta_jm r_i r_k r_l + delta_km r_i r_j r_l
                      + delta_lm r_i r_j r_k)
                  + g5 r_i r_j r_k r_l r_m       */

  CelloView<double, 1> d5_counter(21);

  Hierarchy * hierarchy = enzo::simulation()->hierarchy();
  double lox, loy, loz; 
  double hix, hiy, hiz;
  hierarchy->lower(&lox, &loy, &loz);
  hierarchy->upper(&hix, &hiy, &hiz);

  double Lx = hix - lox;
  double Ly = hiy - loy;
  double Lz = hiz - loz;
  double box_vol = Lx * Ly * Lz;

  double alpha = 2.0 / pow(box_vol, 1./3.);
  double alpha2 = alpha * alpha;
  double alpha3 = alpha2 * alpha;
  double alpha5 = alpha2 * alpha3;
  double alpha7 = alpha2 * alpha5;
  double alpha9 = alpha2 * alpha7;

  int realnxmax = (int)(8.0 / (alpha * Lx) + 0.5);
  int realnymax = (int)(8.0 / (alpha * Ly) + 0.5);
  int realnzmax = (int)(8.0 / (alpha * Lz) + 0.5);

  // sum in real space
  for (int nz = -realnzmax; nz <= realnzmax; nz++) {
    for (int ny = -realnymax; ny <= realnymax; ny++) {
      for (int nx = -realnxmax; nx <= realnxmax; nx++) {

        double rx = x + nx * Lx;
        double ry = y + ny * Ly;
        double rz = z + nz * Lz;
        double r = sqrt(rx*rx + ry*ry + rz*rz);
        double r2 = r*r;
        double r3 = r2*r;
        double r5 = r2*r3;
        double r7 = r2*r5;
        double r9 = r2*r7;
        double r11 = r2*r9;

        double g3 = 0;
        double g4 = 0;
        double g5 = 0;

        if (nx != 0 || ny != 0 || nz != 0) {
          if (r != 0) {
            g3 = (-8.0 * sqrt(M_PI) * alpha5 * r5 - 20.0 * sqrt(M_PI) * alpha3 * r3 - 30.0 * sqrt(M_PI) * alpha * r
            - 15.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r7);

            g4 = (16.0 * sqrt(M_PI) * alpha7 * r7  + 56.0 * sqrt(M_PI) * alpha5 * r5 + 140.0 * sqrt(M_PI) * alpha3 * r3
            + 210.0 * sqrt(M_PI) * alpha * r + 105.0 * M_PI * exp(alpha2 * r2) * erfc(alpha * r)) 
              * exp(-1.0 * alpha2 * r2) / (M_PI * r9);

            g5 = (-32.0 * sqrt(M_PI) * alpha9 * r9 - 144.0 * sqrt(M_PI) * alpha7 * r7  - 504.0 * sqrt(M_PI) * alpha5 * r5 
            - 1260.0 * sqrt(M_PI) * alpha3 * r3 - 1890.0 * sqrt(M_PI) * alpha * r - 945.0 * M_PI * exp(alpha2 * r2) 
            * erfc(alpha * r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r11);
          }
        }
        else {
          
          // Taylor expand if r is close to 0
          if ((alpha * r) < 0.5) {

            g3 = -16.0 * pow(alpha, 7) / sqrt(M_PI) *
                       (-1.0 / 7.0 + pow(alpha * r, 2) / 9.0 - pow(alpha * r, 4) / 22.0 + pow(alpha * r, 6) / 78.0 -
                        pow(alpha * r, 8) / 360.0 + pow(alpha * r, 10) / 2040.0);

            g4 = -32.0 * pow(alpha, 9) / sqrt(M_PI) *
                       (1.0 / 9.0 - pow(alpha * r, 2) / 11.0 + pow(alpha * r, 4) / 26.0 - pow(alpha * r, 6) / 90.0 +
                        pow(alpha * r, 8) / 408.0 - pow(alpha * r, 10) / 2280.0);

            g5 = -64.0 * pow(alpha, 11) / sqrt(M_PI) *
                       (-1.0 / 11.0 + pow(alpha * r, 2) / 13.0 - pow(alpha * r, 4) / 30.0 + pow(alpha * r, 6) / 102.0 -
                        pow(alpha * r, 8) / 456.0 + pow(alpha * r, 10) / 2520.0);

            //if (abs(x + 0.007937) <= 1e-5 && abs(y + 0.007937) <= 1e-5 && abs(z + 0) < 1e-2) {
            //  CkPrintf("D5 Taylor: g3, g4, g5, rx, ry, rz: %f, %f, %f, %f, %f, %f\n", g3, g4, g5, rx, ry, rz);
            //}
          }

          else { // incorporate Newtonian 1/r term
            g3 = 15.0/r7 + (-8.0 * sqrt(M_PI) * alpha5 * r5 - 20.0 * sqrt(M_PI) * alpha3 * r3 - 30.0 * sqrt(M_PI) * alpha * r
                - 15.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r7);
            
            g4 = -105.0/r9 + (16.0 * sqrt(M_PI) * alpha7 * r7  + 56.0 * sqrt(M_PI) * alpha5 * r5 + 140.0 * sqrt(M_PI) * alpha3 * r3
                + 210.0 * sqrt(M_PI) * alpha * r + 105.0 * M_PI * exp(alpha2 * r2) * erfc(alpha * r)) 
                * exp(-1.0 * alpha2 * r2) / (M_PI * r9);

            g5 = 945.0/r11 + (-32.0 * sqrt(M_PI) * alpha9 * r9 - 144.0 * sqrt(M_PI) * alpha7 * r7  - 504.0 * sqrt(M_PI) * alpha5 * r5 
                - 1260.0 * sqrt(M_PI) * alpha3 * r3 - 1890.0 * sqrt(M_PI) * alpha * r - 945.0 * M_PI * exp(alpha2 * r2) 
                * erfc(alpha * r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r11);
          }
        }

        d5_counter(0) -= g5*rx*rx*rx*rx*rx + 10*g4*rx*rx*rx + 15*g3*rx;
        d5_counter(1) -= g5*rx*rx*rx*rx*ry + 6*g4*rx*rx*ry + 3*g3*ry;
        d5_counter(2) -= g5*rx*rx*rx*rx*rz + 6*g4*rx*rx*rz + 3*g3*rz;
        d5_counter(3) -= g5*rx*rx*rx*ry*ry + g4*(3*rx*ry*ry + rx*rx*rx) + 3*g3*rx;
        d5_counter(4) -= g5*rx*rx*rx*ry*rz + 3*g4*rx*ry*rz;
        d5_counter(5) -= g5*rx*rx*rx*rz*rz + g4*(3*rx*rz*rz + rx*rx*rx) + 3*g3*rx;
        d5_counter(6) -= g5*rx*rx*ry*ry*ry + g4*(3*rx*rx*ry + ry*ry*ry) + 3*g3*ry;
        d5_counter(7) -= g5*rx*rx*ry*ry*rz + g4*(ry*ry*rz + rx*rx*rz) + g3*rz;
        d5_counter(8) -= g5*rx*rx*ry*rz*rz + g4*(ry*rz*rz + rx*rx*ry) + g3*ry;
        d5_counter(9) -= g5*rx*rx*rz*rz*rz + g4*(3*rx*rx*rz + rz*rz*rz) + 3*g3*rz;
        d5_counter(10) -= g5*rx*ry*ry*ry*ry + 6*g4*rx*ry*ry + 3*g3*rx;
        d5_counter(11) -= g5*rx*ry*ry*ry*rz + 3*g4*rx*ry*rz;
        d5_counter(12) -= g5*rx*ry*ry*rz*rz + g4*(rx*ry*ry + rx*rz*rz) + g3*rx;
        d5_counter(13) -= g5*rx*ry*rz*rz*rz + 3*g4*rx*ry*rz;
        d5_counter(14) -= g5*rx*rz*rz*rz*rz + 6*g4*rx*rz*rz + 3*g3*rx;
        d5_counter(15) -= g5*ry*ry*ry*ry*ry + 10*g4*ry*ry*ry + 15*g3*ry;
        d5_counter(16) -= g5*ry*ry*ry*ry*rz + 6*g4*ry*ry*rz + 3*g3*rz;
        d5_counter(17) -= g5*ry*ry*ry*rz*rz + g4*(3*ry*rz*rz + ry*ry*ry) + 3*g3*ry;
        d5_counter(18) -= g5*ry*ry*rz*rz*rz + g4*(3*ry*ry*rz + rz*rz*rz) + 3*g3*rz;
        d5_counter(19) -= g5*ry*rz*rz*rz*rz + 6*g4*ry*rz*rz + 3*g3*ry;
        d5_counter(20) -= g5*rz*rz*rz*rz*rz + 10*g4*rz*rz*rz + 15*g3*rz;
      }
    }
  }

  int knxmax = (int)(2.0 * alpha * Lx + 0.5);
  int knymax = (int)(2.0 * alpha * Ly + 0.5);
  int knzmax = (int)(2.0 * alpha * Lz + 0.5);
  
  // sum in Fourier space
  for (int nz = -knzmax; nz <= knzmax; nz++) {
    for (int ny = -knymax; ny <= knymax; ny++) {
      for (int nx = -knxmax; nx <= knxmax; nx++) {

        if (nx != 0 || ny != 0 || nz != 0) {

          double kx = 2.0 * M_PI * nx / Lx;
          double ky = 2.0 * M_PI * ny / Ly;
          double kz = 2.0 * M_PI * nz / Lz;
          double k2 = kx*kx + ky*ky + kz*kz;
          double kdotx = kx*x + ky*y + kz*z;

          double k_exp = 4.0*M_PI/box_vol * (exp(-1.0*k2 / (4.0*alpha2)) / k2) * sin(kdotx); 

          d5_counter(0) += k_exp * kx*kx*kx*kx*kx;
          d5_counter(1) += k_exp * kx*kx*kx*kx*ky;
          d5_counter(2) += k_exp * kx*kx*kx*kx*kz;
          d5_counter(3) += k_exp * kx*kx*kx*ky*ky;
          d5_counter(4) += k_exp * kx*kx*kx*ky*kz;
          d5_counter(5) += k_exp * kx*kx*kx*kz*kz;
          d5_counter(6) += k_exp * kx*kx*ky*ky*ky;
          d5_counter(7) += k_exp * kx*kx*ky*ky*kz;
          d5_counter(8) += k_exp * kx*kx*ky*kz*kz;
          d5_counter(9) += k_exp * kx*kx*kz*kz*kz;
          d5_counter(10) += k_exp * kx*ky*ky*ky*ky;
          d5_counter(11) += k_exp * kx*ky*ky*ky*kz;
          d5_counter(12) += k_exp * kx*ky*ky*kz*kz;
          d5_counter(13) += k_exp * kx*ky*kz*kz*kz;
          d5_counter(14) += k_exp * kx*kz*kz*kz*kz;
          d5_counter(15) += k_exp * ky*ky*ky*ky*ky;
          d5_counter(16) += k_exp * ky*ky*ky*ky*kz;
          d5_counter(17) += k_exp * ky*ky*ky*kz*kz;
          d5_counter(18) += k_exp * ky*ky*kz*kz*kz;
          d5_counter(19) += k_exp * ky*kz*kz*kz*kz;
          d5_counter(20) += k_exp * kz*kz*kz*kz*kz;

        } 
      }
    }
  }

  return d5_counter;

  
}

// compute the d6 term of the Ewald sum at coordinates (x, y, z)
CelloView<double, 1> EnzoEwald::d6(double x, double y, double z) throw()
{
  /* (d6)_ijklmn = g3 (delta_mn delta_ij delta_kl + delta_mn delta_jk delta_il + delta_mn delta_ik delta_jl
                        + delta_ln delta_ij delta_km + delta_ln delta_jk delta_im + delta_ln delta_ik delta_jm
                        + delta_kn delta_ij delta_lm + delta_kn delta_jl delta_im + delta_kn delta_il delta_jm
                        + delta_jn delta_ik delta_lm + delta_jn delta_lk delta_im + delta_jn delta_il delta_km
                        + delta_in delta_jk delta_lm + delta_in delta_lk delta_jm + delta_in delta_jl delta_km)
                +  g4 [r_i r_j (delta_kl delta_mn + delta_lm delta_kn + delta_km delta_ln)
                     + r_i r_k (delta_jl delta_mn + delta_lm delta_jn + delta_jm delta_ln)
                     + r_i r_l (delta_jk delta_mn + delta_km delta_jn + delta_jm delta_kn)
                     + r_i r_m (delta_jk delta_ln + delta_kl delta_jn + delta_jl delta_kn)
                     + r_i r_n (delta_jk delta_lm + delta_kl delta_jm + delta_jl delta_km)
                     + r_j r_k (delta_il delta_mn + delta_lm delta_in + delta_im delta_ln)
                     + r_j r_l (delta_ik delta_mn + delta_km delta_in + delta_im delta_kn)
                     + r_j r_m (delta_ik delta_ln + delta_kl delta_in + delta_il delta_kn)
                     + r_j r_n (delta_ik delta_lm + delta_kl delta_im + delta_il delta_km)
                     + r_k r_l (delta_ij delta_mn + delta_jm delta_in + delta_im delta_jn)
                     + r_k r_m (delta_ij delta_ln + delta_jl delta_in + delta_il delta_jn)
                     + r_k r_n (delta_ij delta_lm + delta_jl delta_im + delta_il delta_jm)
                     + r_l r_m (delta_ij delta_kn + delta_jk delta_in + delta_ik delta_jn)
                     + r_l r_n (delta_ij delta_km + delta_jk delta_im + delta_ik delta_jm)
                     + r_m r_n (delta_ij delta_kl + delta_jk delta_il + delta_ik delta_jl)]
                +  g5 (delta_ij r_k r_l r_m r_n + delta_ik r_j r_l r_m r_n + delta_il r_j r_k r_m r_n
                      + delta_im r_j r_k r_l r_n + delta_in r_j r_k r_l r_m + delta_jk r_i r_l r_m r_n
                      + delta_jl r_i r_k r_m r_n + delta_jm r_i r_k r_l r_n + delta_jn r_i r_k r_l r_m
                      + delta_kl r_i r_j r_m r_n + delta_km r_i r_j r_l r_n + delta_kn r_i r_j r_l r_m
                      + delta_lm r_i r_j r_k r_n + delta_ln r_i r_j r_k r_m + delta_mn r_i r_j r_k r_l)
                +  g6 r_i r_j r_k r_l r_m r_n       
                      */

  CelloView<double, 1> d6_counter(28);

  Hierarchy * hierarchy = enzo::simulation()->hierarchy();
  double lox, loy, loz; 
  double hix, hiy, hiz;
  hierarchy->lower(&lox, &loy, &loz);
  hierarchy->upper(&hix, &hiy, &hiz);

  double Lx = hix - lox;
  double Ly = hiy - loy;
  double Lz = hiz - loz;
  double box_vol = Lx * Ly * Lz;

  double alpha = 2.0 / pow(box_vol, 1./3.);
  double alpha2 = alpha * alpha;
  double alpha3 = alpha2 * alpha;
  double alpha5 = alpha2 * alpha3;
  double alpha7 = alpha2 * alpha5;
  double alpha9 = alpha2 * alpha7;
  double alpha11 = alpha2 * alpha9;

  int realnxmax = (int)(8.0 / (alpha * Lx) + 0.5);
  int realnymax = (int)(8.0 / (alpha * Ly) + 0.5);
  int realnzmax = (int)(8.0 / (alpha * Lz) + 0.5);

  // sum in real space
  for (int nz = -realnzmax; nz <= realnzmax; nz++) {
    for (int ny = -realnymax; ny <= realnymax; ny++) {
      for (int nx = -realnxmax; nx <= realnxmax; nx++) {

        double rx = x + nx * Lx;
        double ry = y + ny * Ly;
        double rz = z + nz * Lz;
        double r = sqrt(rx*rx + ry*ry + rz*rz);
        double r2 = r*r;
        double r3 = r2*r;
        double r5 = r2*r3;
        double r7 = r2*r5;
        double r9 = r2*r7;
        double r11 = r2*r9;
        double r13 = r2*r11;

        double g3 = 0;
        double g4 = 0;
        double g5 = 0;
        double g6 = 0;

        if (nx != 0 || ny != 0 || nz != 0) {
          if (r != 0) {
            g3 = (-8.0 * sqrt(M_PI) * alpha5 * r5 - 20.0 * sqrt(M_PI) * alpha3 * r3 - 30.0 * sqrt(M_PI) * alpha * r
            - 15.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r7);

            g4 = (16.0 * sqrt(M_PI) * alpha7 * r7  + 56.0 * sqrt(M_PI) * alpha5 * r5 + 140.0 * sqrt(M_PI) * alpha3 * r3
            + 210.0 * sqrt(M_PI) * alpha * r + 105.0 * M_PI * exp(alpha2 * r2) * erfc(alpha * r)) 
              * exp(-1.0 * alpha2 * r2) / (M_PI * r9);

            g5 = (-32.0 * sqrt(M_PI) * alpha9 * r9 - 144.0 * sqrt(M_PI) * alpha7 * r7  - 504.0 * sqrt(M_PI) * alpha5 * r5 
            - 1260.0 * sqrt(M_PI) * alpha3 * r3 - 1890.0 * sqrt(M_PI) * alpha * r - 945.0 * M_PI * exp(alpha2 * r2) 
            * erfc(alpha * r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r11);

            g6 = (64.0 * sqrt(M_PI) * alpha11 * r11 + 352.0 * sqrt(M_PI) * alpha9 * r9 + 1584.0 * sqrt(M_PI) * alpha7 * r7  
            + 5544.0 * sqrt(M_PI) * alpha5 * r5 + 13860.0 * sqrt(M_PI) * alpha3 * r3 + 20790.0 * sqrt(M_PI) * alpha * r 
            + 10395.0 * M_PI * exp(alpha2 * r2) * erfc(alpha * r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r13);
          }
        }
        else {
          
          // Taylor expand if r is close to 0
          if ((alpha * r) < 0.5) {
            g3 = -16.0 * pow(alpha, 7) / sqrt(M_PI) *
                       (-1.0 / 7.0 + pow(alpha * r, 2) / 9.0 - pow(alpha * r, 4) / 22.0 + pow(alpha * r, 6) / 78.0 -
                        pow(alpha * r, 8) / 360.0 + pow(alpha * r, 10) / 2040.0);

            g4 = -32.0 * pow(alpha, 9) / sqrt(M_PI) *
                       (1.0 / 9.0 - pow(alpha * r, 2) / 11.0 + pow(alpha * r, 4) / 26.0 - pow(alpha * r, 6) / 90.0 +
                        pow(alpha * r, 8) / 408.0 - pow(alpha * r, 10) / 2280.0);

            g5 = -64.0 * pow(alpha, 11) / sqrt(M_PI) *
                       (-1.0 / 11.0 + pow(alpha * r, 2) / 13.0 - pow(alpha * r, 4) / 30.0 + pow(alpha * r, 6) / 102.0 -
                        pow(alpha * r, 8) / 456.0 + pow(alpha * r, 10) / 2520.0);

            g6 = -128.0 * pow(alpha, 13) / sqrt(M_PI) *
                       (1.0 / 13.0 - pow(alpha * r, 2) / 15.0 + pow(alpha * r, 4) / 34.0 - pow(alpha * r, 6) / 114.0 +
                        pow(alpha * r, 8) / 504.0 - pow(alpha * r, 10) / 2760.0);


            //if (abs(x + 0.007937) <= 1e-5 && abs(y + 0.007937) <= 1e-5 && abs(z + 0) < 1e-2) {
            //  CkPrintf("D6 Taylor: g3, g4, g5, g6, rx, ry, rz: %f, %f, %f, %f, %f, %f, %f\n", g3, g4, g5, g6, rx, ry, rz);
            //}

          }

          else { // incorporate Newtonian 1/r term
            g3 = 15.0/r7 + (-8.0 * sqrt(M_PI) * alpha5 * r5 - 20.0 * sqrt(M_PI) * alpha3 * r3 - 30.0 * sqrt(M_PI) * alpha * r
                - 15.0 * M_PI * exp(alpha2 * r2) * erfc(alpha*r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r7);

            g4 = -105.0/r9 + (16.0 * sqrt(M_PI) * alpha7 * r7  + 56.0 * sqrt(M_PI) * alpha5 * r5 + 140.0 * sqrt(M_PI) * alpha3 * r3
                + 210.0 * sqrt(M_PI) * alpha * r + 105.0 * M_PI * exp(alpha2 * r2) * erfc(alpha * r)) 
                * exp(-1.0 * alpha2 * r2) / (M_PI * r9);

            g5 = 945.0/r11 + (-32.0 * sqrt(M_PI) * alpha9 * r9 - 144.0 * sqrt(M_PI) * alpha7 * r7  - 504.0 * sqrt(M_PI) * alpha5 * r5 
                - 1260.0 * sqrt(M_PI) * alpha3 * r3 - 1890.0 * sqrt(M_PI) * alpha * r - 945.0 * M_PI * exp(alpha2 * r2) 
                * erfc(alpha * r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r11);

            g6 = -10395.0/r13 + (64.0 * sqrt(M_PI) * alpha11 * r11 + 352.0 * sqrt(M_PI) * alpha9 * r9 + 1584.0 * sqrt(M_PI) * alpha7 * r7  
                + 5544.0 * sqrt(M_PI) * alpha5 * r5 + 13860.0 * sqrt(M_PI) * alpha3 * r3 + 20790.0 * sqrt(M_PI) * alpha * r 
                + 10395.0 * M_PI * exp(alpha2 * r2) * erfc(alpha * r)) * exp(-1.0 * alpha2 * r2) / (M_PI * r13);
          }
        }

        d6_counter(0) -= g6*rx*rx*rx*rx*rx*rx + 15*g5*rx*rx*rx*rx + 45*g4*rx*rx + 15*g3;
        d6_counter(1) -= g6*rx*rx*rx*rx*rx*ry + 10*g5*rx*rx*rx*ry + 15*g4*rx*ry; 
        d6_counter(2) -= g6*rx*rx*rx*rx*rx*rz + 10*g5*rx*rx*rx*rz + 15*g4*rx*rz;
        d6_counter(3) -= g6*rx*rx*rx*rx*ry*ry + g5*(6*rx*rx*ry*ry + rx*rx*rx*rx) + g4*(6*rx*rx + 3*ry*ry) + 3*g3; 
        d6_counter(4) -= g6*rx*rx*rx*rx*ry*rz + 6*g5*rx*rx*ry*rz + 3*g4*ry*rz;
        d6_counter(5) -= g6*rx*rx*rx*rx*rz*rz + g5*(6*rx*rx*rz*rz + rx*rx*rx*rx) + g4*(6*rx*rx + 3*rz*rz) + 3*g3;
        d6_counter(6) -= g6*rx*rx*rx*ry*ry*ry + g5*(3*rx*ry*ry*ry + 3*rx*rx*rx*ry) + 9*g4*rx*ry;
        d6_counter(7) -= g6*rx*rx*rx*ry*ry*rz + g5*(3*rx*ry*ry*rz + rx*rx*rx*rz) + 3*g4*rx*rz;
        d6_counter(8) -= g6*rx*rx*rx*ry*rz*rz + g5*(3*rx*ry*rz*rz + rx*rx*rx*ry) + 3*g4*rx*ry;
        d6_counter(9) -= g6*rx*rx*rx*rz*rz*rz + g5*(3*rx*rz*rz*rz + 3*rx*rx*rx*rz) + 9*g4*rx*rz;
        d6_counter(10) -= g6*rx*rx*ry*ry*ry*ry + g5*(6*rx*rx*ry*ry + ry*ry*ry*ry) + g4*(6*ry*ry + 3*rx*rx) + 3*g3;
        d6_counter(11) -= g6*rx*rx*ry*ry*ry*rz + g5*(3*rx*rx*ry*rz + ry*ry*ry*rz) + 3*g4*ry*rz;
        d6_counter(12) -= g6*rx*rx*ry*ry*rz*rz + g5*(ry*ry*rz*rz + rx*rx*rz*rz + rx*rx*ry*ry) + g4*(rx*rx + ry*ry + rz*rz) + g3;
        d6_counter(13) -= g6*rx*rx*ry*rz*rz*rz + g5*(3*rx*rx*ry*rz + ry*rz*rz*rz) + 3*g4*ry*rz;
        d6_counter(14) -= g6*rx*rx*rz*rz*rz*rz + g5*(6*rx*rx*rz*rz + rz*rz*rz*rz) + g4*(6*rz*rz + 3*rx*rx) + 3*g3;
        d6_counter(15) -= g6*rx*ry*ry*ry*ry*ry + 10*g5*rx*ry*ry*ry + 15*g4*rx*ry;
        d6_counter(16) -= g6*rx*ry*ry*ry*ry*rz + 6*g5*rx*ry*ry*rz + 3*g4*rx*rz;
        d6_counter(17) -= g6*rx*ry*ry*ry*rz*rz + g5*(3*rx*ry*rz*rz + rx*ry*ry*ry) + 3*g4*rx*ry;
        d6_counter(18) -= g6*rx*ry*ry*rz*rz*rz + g5*(3*rx*ry*ry*rz + rx*rz*rz*rz) + 3*g4*rx*rz;
        d6_counter(19) -= g6*rx*ry*rz*rz*rz*rz + 6*g5*rx*ry*rz*rz + 3*g4*rx*ry;
        d6_counter(20) -= g6*rx*rz*rz*rz*rz*rz + 10*g5*rx*rz*rz*rz + 15*g4*rx*rz;
        d6_counter(21) -= g6*ry*ry*ry*ry*ry*ry + 15*g5*ry*ry*ry*ry + 45*g4*ry*ry + 15*g3;
        d6_counter(22) -= g6*ry*ry*ry*ry*ry*rz + 10*g5*ry*ry*ry*rz + 15*g4*ry*rz;
        d6_counter(23) -= g6*ry*ry*ry*ry*rz*rz + g5*(6*ry*ry*rz*rz + ry*ry*ry*ry) + g4*(6*ry*ry + 3*rz*rz) + 3*g3;
        d6_counter(24) -= g6*ry*ry*ry*rz*rz*rz + g5*(3*ry*rz*rz*rz + 3*ry*ry*ry*rz) + 9*g4*ry*rz;
        d6_counter(25) -= g6*ry*ry*rz*rz*rz*rz + g5*(6*ry*ry*rz*rz + rz*rz*rz*rz) + g4*(6*rz*rz + 3*ry*ry) + 3*g3;
        d6_counter(26) -= g6*ry*rz*rz*rz*rz*rz + 10*g5*ry*rz*rz*rz + 15*g4*ry*rz;
        d6_counter(27) -= g6*rz*rz*rz*rz*rz*rz + 15*g5*rz*rz*rz*rz + 45*g4*rz*rz + 15*g3;
        
      }
    }
  }

  int knxmax = (int)(2.0 * alpha * Lx + 0.5);
  int knymax = (int)(2.0 * alpha * Ly + 0.5);
  int knzmax = (int)(2.0 * alpha * Lz + 0.5);
  
  // sum in Fourier space
  for (int nz = -knzmax; nz <= knzmax; nz++) {
    for (int ny = -knymax; ny <= knymax; ny++) {
      for (int nx = -knxmax; nx <= knxmax; nx++) {

        if (nx != 0 || ny != 0 || nz != 0) {

          double kx = 2.0 * M_PI * nx / Lx;
          double ky = 2.0 * M_PI * ny / Ly;
          double kz = 2.0 * M_PI * nz / Lz;
          double k2 = kx*kx + ky*ky + kz*kz;
          double kdotx = kx*x + ky*y + kz*z;

          double k_exp = 4.0*M_PI/box_vol * (exp(-1.0*k2 / (4.0*alpha2)) / k2) * cos(kdotx); 

          d6_counter(0) += k_exp * kx*kx*kx*kx*kx*kx;
          d6_counter(1) += k_exp * kx*kx*kx*kx*kx*ky;
          d6_counter(2) += k_exp * kx*kx*kx*kx*kx*kz;
          d6_counter(3) += k_exp * kx*kx*kx*kx*ky*ky;
          d6_counter(4) += k_exp * kx*kx*kx*kx*ky*kz;
          d6_counter(5) += k_exp * kx*kx*kx*kx*kz*kz;
          d6_counter(6) += k_exp * kx*kx*kx*ky*ky*ky;
          d6_counter(7) += k_exp * kx*kx*kx*ky*ky*kz;
          d6_counter(8) += k_exp * kx*kx*kx*ky*kz*kz;
          d6_counter(9) += k_exp * kx*kx*kx*kz*kz*kz;
          d6_counter(10) += k_exp * kx*kx*ky*ky*ky*ky;
          d6_counter(11) += k_exp * kx*kx*ky*ky*ky*kz;
          d6_counter(12) += k_exp * kx*kx*ky*ky*kz*kz;
          d6_counter(13) += k_exp * kx*kx*ky*kz*kz*kz;
          d6_counter(14) += k_exp * kx*kx*kz*kz*kz*kz;
          d6_counter(15) += k_exp * kx*ky*ky*ky*ky*ky;
          d6_counter(16) += k_exp * kx*ky*ky*ky*ky*kz;
          d6_counter(17) += k_exp * kx*ky*ky*ky*kz*kz;
          d6_counter(18) += k_exp * kx*ky*ky*kz*kz*kz;
          d6_counter(19) += k_exp * kx*ky*kz*kz*kz*kz;
          d6_counter(20) += k_exp * kx*kz*kz*kz*kz*kz;
          d6_counter(21) += k_exp * ky*ky*ky*ky*ky*ky;
          d6_counter(22) += k_exp * ky*ky*ky*ky*ky*kz;
          d6_counter(23) += k_exp * ky*ky*ky*ky*kz*kz;
          d6_counter(24) += k_exp * ky*ky*ky*kz*kz*kz;
          d6_counter(25) += k_exp * ky*ky*kz*kz*kz*kz;
          d6_counter(26) += k_exp * ky*kz*kz*kz*kz*kz;
          d6_counter(27) += k_exp * kz*kz*kz*kz*kz*kz;
        } 
      }
    }
  }

  return d6_counter;
  
}

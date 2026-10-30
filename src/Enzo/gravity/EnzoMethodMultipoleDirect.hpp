// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMethodMultipoleDirect.hpp
/// @brief    [\ref Enzo] Direct (leaf-leaf) gravitational interaction kernels
///           used by EnzoMethodMultipole::interact_direct_.
///
/// This header is deliberately free of Charm++/Cello dependencies so that the
/// kernels can be tested and benchmarked standalone.
///
/// All "accumulate" functions ADD  grav_const * sum_j m_j * K(x_j - x_i)
/// to the output accelerations, where K is the softened Newtonian kernel plus
/// (optionally) the Ewald correction for periodic boundaries.
///
/// The cell-cell sums are the original brute-force loops.

#ifndef ENZO_METHOD_MULTIPOLE_DIRECT_HPP
#define ENZO_METHOD_MULTIPOLE_DIRECT_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace fmm_direct {

//----------------------------------------------------------------------

struct Params {
  double eps0 = 0.0;           // gravitational softening length
  double grav_const = 1.0;     // includes 1/(4 pi a^2) in cosmology
  int    rank = 3;
  bool   periodic[3] = {false, false, false};
  double dom_width[3] = {1.0, 1.0, 1.0};
};

/// Ewald functor placeholder for non-periodic runs
struct NoEwald {
  static constexpr bool enabled = false;
  std::array<double,3> operator() (double, double, double) const
  { return {0.0, 0.0, 0.0}; }
  std::array<double,6>  d2 (double, double, double) const { return {}; }
  std::array<double,10> d3 (double, double, double) const { return {}; }
  const void * id () const { return nullptr; }
};

/// Geometry and density of the cells of one block. `dens` is the full
/// padded array (including ghost zones); nx,ny,nz count interior cells.
struct CellGrid {
  int nx, ny, nz;
  int mx, my, mz;
  int gx, gy, gz;
  double h[3];
  double lo[3];
  const double * dens;

  double cell_volume() const { return h[0]*h[1]*h[2]; }
  double pos(int axis, int i) const { return lo[axis] + (i + 0.5)*h[axis]; }
  const double * row (int j, int k) const
  { return dens + gx + mx*((j+gy) + my*(size_t)(k+gz)); }
};

//----------------------------------------------------------------------
// Scalar helpers

/// Softening kernel (gradient of potential; Springel, Yoshida & White 2001)
inline double W2 (double u)
{
  if (u < 0.5)
    return 32./3 * u - 192./5 * u*u*u + 32. * u*u*u*u;
  else if (u < 1.0)
    return -1./15 / (u*u) + 64./3 * u - 48. * u*u + 192./5 * u*u*u
           - 32./3 * u*u*u*u;
  else
    return 1./(u*u);
}

/// Acceleration at the origin due to a unit mass at displacement d
/// (softened Newtonian, zero at d = 0)
inline void newton_unit (double h, double dx, double dy, double dz,
                         double & ax, double & ay, double & az)
{
  const double r2 = dx*dx + dy*dy + dz*dz;
  if (r2 == 0.0) { ax = ay = az = 0.0; return; }
  const double r = std::sqrt(r2);
  double s;                      // acceleration magnitude / r
  if (r < h) s = W2(r/h) / (h*h) / r;
  else       s = 1.0 / (r2*r);
  ax = s*dx; ay = s*dy; az = s*dz;
}

/// Component of the displacement (image of b nearest to a) - a
inline double wrapped_disp (double xb, double xa, double width, bool periodic)
{
  if (!periodic) return xb - xa;
  const double q = (xb - xa) / width;
  const double n = std::floor(q);
  return (((q - n) < 0.5) ? (xb - n*width) : (xb - (n + 1.0)*width)) - xa;
}

/// Displacement from a to b. With Ewald (periodic) the nearest periodic
/// image is used and unused dimensions are zeroed; otherwise it is the raw
/// difference (this mirrors the original code).
template <bool EW>
inline void displacement (const Params & p, const double * a, const double * b,
                          double * d)
{
  if (EW) {
    d[0] = wrapped_disp(b[0], a[0], p.dom_width[0], p.periodic[0]);
    d[1] = (p.rank < 2) ? 0.0 :
           wrapped_disp(b[1], a[1], p.dom_width[1], p.periodic[1]);
    d[2] = (p.rank < 3) ? 0.0 :
           wrapped_disp(b[2], a[2], p.dom_width[2], p.periodic[2]);
  } else {
    d[0] = b[0] - a[0]; d[1] = b[1] - a[1]; d[2] = b[2] - a[2];
  }
}

/// Length h of the gravitational softening: the force is exactly Newtonian
/// beyond distance h = 2.8 eps0 (Springel, Yoshida & White 2001)
inline double softening_length (const Params & p) { return 2.8*p.eps0; }

/// Kernel per unit source mass: Newton(d) + Ewald(d)
template <class Ew>
inline void kernel (const Params & p, const Ew & ew, const double * d,
                    double * k)
{
  const double h = softening_length(p);
  newton_unit(h, d[0], d[1], d[2], k[0], k[1], k[2]);
  if (Ew::enabled) {
    const std::array<double,3> e = ew(d[0], d[1], d[2]);
    k[0] += e[0]; k[1] += e[1]; k[2] += e[2];
  }
}

//----------------------------------------------------------------------
// Point-point accumulation (particles, or cells treated as points)

/// out[t] += grav_const * sum_s sm[s] * K(s - t)
/// If `exclude` is given, source exclude[t] (if >= 0) is skipped for target t
/// (used to leave out a particle's interaction with itself).
template <class Ew>
void accumulate_points (const Params & p, const Ew & ew,
                        int nt, const double * tx, const double * ty,
                        const double * tz,
                        int ns, const double * sx, const double * sy,
                        const double * sz, const double * sm,
                        double * ox, double * oy, double * oz,
                        const int * exclude = nullptr)
{
  for (int t = 0; t < nt; t++) {
    const double a[3] = {tx[t], ty[t], tz[t]};
    const int skip = exclude ? exclude[t] : -1;
    double sumx = 0, sumy = 0, sumz = 0;
    for (int s = 0; s < ns; s++) {
      if (s == skip) continue;
      const double b[3] = {sx[s], sy[s], sz[s]};
      double d[3], k[3];
      displacement<Ew::enabled>(p, a, b, d);
      kernel(p, ew, d, k);
      sumx += sm[s]*k[0]; sumy += sm[s]*k[1]; sumz += sm[s]*k[2];
    }
    ox[t] += p.grav_const*sumx;
    oy[t] += p.grav_const*sumy;
    oz[t] += p.grav_const*sumz;
  }
}

/// Cell centers and masses as point lists (interior cells only, x fastest)
inline void cells_to_points (const CellGrid & g,
                             std::vector<double> & x, std::vector<double> & y,
                             std::vector<double> & z, std::vector<double> & m)
{
  const size_t n = (size_t)g.nx*g.ny*g.nz;
  x.resize(n); y.resize(n); z.resize(n); m.resize(n);
  const double vol = g.cell_volume();
  size_t c = 0;
  for (int k = 0; k < g.nz; k++)
    for (int j = 0; j < g.ny; j++) {
      const double * d = g.dens ? g.row(j,k) : nullptr;   // null: positions only
      for (int i = 0; i < g.nx; i++, c++) {
        x[c] = g.pos(0,i); y[c] = g.pos(1,j); z[c] = g.pos(2,k);
        m[c] = d ? d[i]*vol : 0.0;
      }
    }
}

//----------------------------------------------------------------------
// Cell-cell accumulation

/// The original brute-force six-deep loop. accel arrays are padded like
/// `a.dens`. If `self` is true, a and b are the same block and each cell's
/// interaction with itself is left out.
template <class Ew, class real>
void accumulate_cells_cells (const Params & p, const Ew & ew,
                      const CellGrid & a, const CellGrid & b,
                      real * ax, real * ay, real * az, bool self = false)
{
  const double vol_b = b.cell_volume();
  for (int k = 0; k < a.nz; k++)
  for (int j = 0; j < a.ny; j++)
  for (int i = 0; i < a.nx; i++) {
    const size_t ia = (i + a.gx) + a.mx*((size_t)(j + a.gy) + (size_t)a.my*(k + a.gz));
    const double pa[3] = {a.pos(0,i), a.pos(1,j), a.pos(2,k)};
    for (int k2 = 0; k2 < b.nz; k2++)
    for (int j2 = 0; j2 < b.ny; j2++)
    for (int i2 = 0; i2 < b.nx; i2++) {
      if (self && i2 == i && j2 == j && k2 == k) continue;
      const double m = b.row(j2,k2)[i2]*vol_b;
      const double pb[3] = {b.pos(0,i2), b.pos(1,j2), b.pos(2,k2)};
      double d[3], kv[3];
      displacement<Ew::enabled>(p, pa, pb, d);
      kernel(p, ew, d, kv);
      ax[ia] += p.grav_const*m*kv[0];
      ay[ia] += p.grav_const*m*kv[1];
      az[ia] += p.grav_const*m*kv[2];
    }
  }
}

//----------------------------------------------------------------------
// Whole block-pair sums (cells and particles)

/// A list of points (particles); `m` is not used for target lists
struct PointSet {
  int n = 0;
  const double * x = nullptr;
  const double * y = nullptr;
  const double * z = nullptr;
  const double * m = nullptr;
};

/// Add the acceleration on all cells of block `a` and on the target points
/// `tgt` (output tox,toy,toz, one entry per target) due to the cells of
/// block `b` and the source points `src`. If `self`, a and b are the same
/// block and each cell/particle's interaction with itself is omitted
/// (`exclude[t]` is the source index of target particle t, or -1).
/// All outputs are scaled by p.grav_const.
template <class Ew, class real>
void accumulate_block_pair_direct
(const Params & p, const Ew & ew, const CellGrid & a, const CellGrid & b,
 bool self, real * ax, real * ay, real * az,
 const PointSet & tgt, double * tox, double * toy, double * toz,
 const PointSet & src, const int * exclude)
{
  // cells <- cells
  accumulate_cells_cells(p, ew, a, b, ax, ay, az, self);

  // cells <- particles
  if (src.n > 0) {
    std::vector<double> cx, cy, cz, cm;
    cells_to_points(a, cx, cy, cz, cm);
    const size_t nc = cx.size();
    std::vector<double> ox(nc,0.0), oy(nc,0.0), oz(nc,0.0);
    accumulate_points(p, ew, (int)nc, cx.data(), cy.data(), cz.data(),
                      src.n, src.x, src.y, src.z, src.m,
                      ox.data(), oy.data(), oz.data());
    size_t c = 0;
    for (int k = 0; k < a.nz; k++)
    for (int j = 0; j < a.ny; j++)
    for (int i = 0; i < a.nx; i++, c++) {
      const size_t ia = (i + a.gx) + a.mx*((size_t)(j + a.gy) + (size_t)a.my*(k + a.gz));
      ax[ia] += ox[c]; ay[ia] += oy[c]; az[ia] += oz[c];
    }
  }

  if (tgt.n == 0) return;

  // particles <- cells
  {
    std::vector<double> cx, cy, cz, cm;
    cells_to_points(b, cx, cy, cz, cm);
    accumulate_points(p, ew, tgt.n, tgt.x, tgt.y, tgt.z,
                      (int)cx.size(), cx.data(), cy.data(), cz.data(), cm.data(),
                      tox, toy, toz);
  }

  // particles <- particles
  if (src.n > 0) {
    accumulate_points(p, ew, tgt.n, tgt.x, tgt.y, tgt.z,
                      src.n, src.x, src.y, src.z, src.m,
                      tox, toy, toz, self ? exclude : nullptr);
  }
}

} // namespace fmm_direct

#endif /* ENZO_METHOD_MULTIPOLE_DIRECT_HPP */

// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMethodMultipoleSubtree.hpp
/// @brief    [\ref Enzo] Tree-based (FMM-style) evaluation of the leaf-leaf
///           interaction, done locally within a pair of leaf Blocks.
///
/// Each leaf Block is divided into chunks of `chunk`^3 cells (e.g. 2x2x2),
/// which form the leaves of a small local octree (a Block of n^3 cells has
/// log2(n/chunk) levels below it). Cells *and* particles belong to the chunk
/// containing them, so both contribute to that chunk's multipole moments.
///
/// For a pair of Blocks (A receives acceleration from B), a dual tree walk
/// over tree(A) x tree(B) either
///   - accepts a pair of nodes as well separated (opening angle `theta`) and
///     adds the multipole-to-local (M2L) interaction to the local expansion of
///     the node in A, using exactly the same expansion (order 3 in the
///     acceleration, quadrupole moments) as the global FMM, including the Ewald
///     derivative tensors when periodic; or
///   - reaches a pair of chunks that are too close, which interact directly
///     (cells and particles, using the same kernel as the direct sums).
/// The local expansions are then shifted down the tree (L2L) and evaluated at
/// every cell and particle of A. No communication is needed.
///
/// With theta = 0 nothing is ever well separated and the result is identical
/// to the direct sums (up to rounding); this checks the bookkeeping and the direct
/// part, the error of the expansions is checked by decreasing theta.

#ifndef ENZO_METHOD_MULTIPOLE_SUBTREE_HPP
#define ENZO_METHOD_MULTIPOLE_SUBTREE_HPP

#include "EnzoMethodMultipoleDirect.hpp"
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <unordered_map>

namespace fmm_direct {

struct SubtreeOpts {
  int chunk = 2;         // cells per chunk along each axis (leaf of local tree)
  double theta = 0.5;    // opening angle: (ra+rb) < theta * distance
  bool timing = false;   // accumulate per-phase timings in subtree::stats()
};

namespace subtree {

//----------------------------------------------------------------------
// Packed symmetric tensor indexing (same ordering as EnzoMethodMultipole:
// 00 01 02 11 12 22 and 000 001 002 011 012 022 111 112 122 222)

struct Idx {
  int i2[3][3];
  int i3[3][3][3];
  Idx () {
    int m3[3][3][3] = {};
    int cnt = 0;
    for (int a = 0; a < 3; a++) for (int b = a; b < 3; b++) for (int c = b; c < 3; c++)
      m3[a][b][c] = cnt++;
    const int m2[3][3] = {{0,1,2},{1,3,4},{2,4,5}};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
      i2[i][j] = m2[i][j];
      for (int k = 0; k < 3; k++) {
        int t[3] = {i,j,k};
        if (t[0] > t[1]) std::swap(t[0],t[1]);
        if (t[1] > t[2]) std::swap(t[1],t[2]);
        if (t[0] > t[1]) std::swap(t[0],t[1]);
        i3[i][j][k] = m3[t[0]][t[1]][t[2]];
      }
    }
  }
};
inline const Idx & idx () { static const Idx t; return t; }

/// Derivative tensors of 1/r at displacement r (d1 = -r/r^3, ...)
inline void newton_tensors (const double * r, double * d1, double * d2, double * d3)
{
  const double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
  const double r3 = r2*std::sqrt(r2), r5 = r2*r3, r7 = r2*r5;
  const double x = r[0], y = r[1], z = r[2];
  d1[0] = -x/r3; d1[1] = -y/r3; d1[2] = -z/r3;
  d2[0] = 3.0/r5*x*x - 1.0/r3;  d2[1] = 3.0/r5*x*y;  d2[2] = 3.0/r5*x*z;
  d2[3] = 3.0/r5*y*y - 1.0/r3;  d2[4] = 3.0/r5*y*z;  d2[5] = 3.0/r5*z*z - 1.0/r3;
  d3[0] = -15.0/r7*x*x*x + 9.0/r5*x;  d3[1] = -15.0/r7*x*x*y + 3.0/r5*y;
  d3[2] = -15.0/r7*x*x*z + 3.0/r5*z;  d3[3] = -15.0/r7*x*y*y + 3.0/r5*x;
  d3[4] = -15.0/r7*x*y*z;             d3[5] = -15.0/r7*x*z*z + 3.0/r5*x;
  d3[6] = -15.0/r7*y*y*y + 9.0/r5*y;  d3[7] = -15.0/r7*y*y*z + 3.0/r5*z;
  d3[8] = -15.0/r7*y*z*z + 3.0/r5*y;  d3[9] = -15.0/r7*z*z*z + 9.0/r5*z;
}

//----------------------------------------------------------------------

/// Node of the local tree. The moments describe the mass inside the node
/// about the node's geometric center s: total mass m, dipole P = sum m (x-s)
/// and second moment Q = sum m (x-s)(x-s). Expanding about the geometric
/// center (rather than the center of mass) makes the displacement between two
/// nodes a lattice vector, so the M2L derivative tensors can be cached. The
/// coefficients c1,c2,c3 are those of the local expansion about the center
/// of the node: the acceleration at center + a is  -( c1 - c2.a + 1/2 c3:aa ).
struct TNode {
  double m = 0;
  double P[3] = {0,0,0};
  double Q[6] = {0,0,0,0,0,0};
  double c1[3] = {0,0,0};
  double c2[6] = {0,0,0,0,0,0};
  double c3[10] = {0,0,0,0,0,0,0,0,0,0};
};

/// Derivative tensors (d1,d2,d3) of the full periodic Green's function at a
/// lattice displacement, cached by (level, wrapped integer offset)
struct M2LTensors { double d1[3], d2[6], d3[10]; };

struct M2LKey {
  unsigned long long scale, pk;
  bool operator== (const M2LKey & o) const { return scale == o.scale && pk == o.pk; }
};
struct M2LKeyHash {
  size_t operator() (const M2LKey & k) const {
    unsigned long long h = k.scale * 0x9E3779B97F4A7C15ULL ^ (k.pk + 0x7F4A7C15ULL + (k.scale << 6));
    h ^= h >> 29; h *= 0xBF58476D1CE4E5B9ULL; h ^= h >> 32;
    return (size_t) h;
  }
};

struct M2LCache {
  const void * id = nullptr;
  std::unordered_map<M2LKey, M2LTensors, M2LKeyHash> map;
};
inline M2LCache & m2l_cache () { static thread_local M2LCache c; return c; }

/// Kernel between cells of blocks of different resolution, tabulated by the
/// wrapped integer offset in units of half the fine cell width. Such pairs are
/// rare (refinement interfaces), so a hash map is adequate.
struct MixedKernelCache {
  const void * id = nullptr;
  double eps0 = -1, u[3] = {0,0,0}, width[3] = {0,0,0};
  bool periodic[3] = {false,false,false};
  std::unordered_map<unsigned long long, std::array<double,3> > map;
  bool matches (const void * id_, const Params & p, const double * u_) const {
    if (id_ != id || p.eps0 != eps0) return false;
    for (int d = 0; d < 3; d++)
      if (u_[d] != u[d] || p.periodic[d] != periodic[d] || p.dom_width[d] != width[d]) return false;
    return true;
  }
  void reset (const void * id_, const Params & p, const double * u_) {
    id = id_; eps0 = p.eps0; map.clear();
    for (int d = 0; d < 3; d++) { u[d] = u_[d]; periodic[d] = p.periodic[d]; width[d] = p.dom_width[d]; }
  }
};
inline MixedKernelCache & mixed_kernel_cache () { static thread_local MixedKernelCache c; return c; }

/// Kernel (Newton + Ewald, per unit mass) between two cells, tabulated by the
/// wrapped integer offset between them. The near field only ever needs
/// small offsets, so a small dense table suffices however large the domain.
struct NearKernelTable {
  static constexpr int R = 40, S = 2*R + 1;
  const void * id = nullptr;
  double eps0 = -1, h[3] = {0,0,0};
  bool periodic[3] = {false,false,false};
  double width[3] = {0,0,0};
  std::vector<double> K;      // 3 per entry
  std::vector<char> done;
  bool matches (const void * id_, const Params & p, const double * h_) const {
    if (id_ != id || p.eps0 != eps0 || done.empty()) return false;
    for (int d = 0; d < 3; d++)
      if (h_[d] != h[d] || p.periodic[d] != periodic[d] || p.dom_width[d] != width[d])
        return false;
    return true;
  }
  void reset (const void * id_, const Params & p, const double * h_) {
    id = id_; eps0 = p.eps0;
    for (int d = 0; d < 3; d++) { h[d] = h_[d]; periodic[d] = p.periodic[d]; width[d] = p.dom_width[d]; }
    K.assign(3*(size_t)S*S*S, 0.0); done.assign((size_t)S*S*S, 0);
  }
};
inline NearKernelTable & near_table () { static thread_local NearKernelTable t; return t; }

/// Accumulated timings and counts of the tree solver (diagnostics)
struct Stats {
  double build = 0, walk = 0, down = 0, eval = 0;       // phases of run()
  double p2p_cc = 0, p2p_cp = 0, p2p_pc = 0, p2p_pp = 0; // within the walk
  long long m2l = 0, p2p_leaf = 0, calls = 0;
  void clear () { *this = Stats(); }
};
inline Stats & stats () { static thread_local Stats s; return s; }

inline double wall_seconds ()
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

struct Tree {
  int K = 0;                      // leaf level (root is level 0)
  double lo[3], width[3];         // lower corner and extent of the Block
  std::vector<size_t> base;       // first node of each level
  std::vector<TNode> nodes;

  void init (int K_, const double * lo_, const double * width_) {
    K = K_;
    for (int d = 0; d < 3; d++) { lo[d] = lo_[d]; width[d] = width_[d]; }
    base.assign(K+2, 0);
    for (int k = 0; k <= K; k++) { const size_t g = size_t(1) << k; base[k+1] = base[k] + g*g*g; }
    nodes.assign(base[K+1], TNode());
  }
  size_t id (int k, int i, int j, int l) const
  { const size_t g = size_t(1) << k; return base[k] + i + g*(j + g*l); }
  void center (int k, int i, int j, int l, double * c) const {
    const int t[3] = {i,j,l};
    for (int d = 0; d < 3; d++) c[d] = lo[d] + (t[d] + 0.5)*width[d]/double(1 << k);
  }
  /// half the diagonal of a node at level k
  double radius (int k) const {
    double s = 0;
    for (int d = 0; d < 3; d++) { const double w = width[d]/double(1 << k); s += w*w; }
    return 0.5*std::sqrt(s);
  }
};

/// Points bucketed by the leaf chunk containing them
struct Buckets {
  std::vector<int> start, items;
  void build (const PointSet & pts, const Tree & t) {
    const int g = 1 << t.K;
    const size_t nl = size_t(g)*g*g;
    start.assign(nl + 1, 0);
    std::vector<int> leaf(pts.n);
    for (int q = 0; q < pts.n; q++) {
      const double x[3] = {pts.x[q], pts.y[q], pts.z[q]};
      int c[3];
      for (int d = 0; d < 3; d++) {
        c[d] = (int) std::floor((x[d] - t.lo[d]) / (t.width[d]/g));
        c[d] = std::min(std::max(c[d], 0), g-1);
      }
      leaf[q] = c[0] + g*(c[1] + g*c[2]);
      start[leaf[q] + 1]++;
    }
    for (size_t l = 0; l < nl; l++) start[l+1] += start[l];
    items.resize(pts.n);
    std::vector<int> fill(start.begin(), start.end()-1);
    for (int q = 0; q < pts.n; q++) items[fill[leaf[q]]++] = q;
  }
};

//----------------------------------------------------------------------

template <class Ew, class real>
class Solver {

public:

  Solver (const Params & p, const Ew & ew, const SubtreeOpts & o,
          const CellGrid & a, const CellGrid & b, bool self,
          const PointSet & tgt, const PointSet & src, const int * exclude)
    : p_(p), ew_(ew), o_(o), a_(a), b_(b), self_(self), tgt_(tgt), src_(src),
      exclude_(exclude), hsoft_(softening_length(p)) {}

  /// returns false (having changed nothing) if this pair of Blocks is not
  /// supported, in which case the caller should use the direct sums; why()
  /// then says why
  bool run (real * ax, real * ay, real * az,
            double * tox, double * toy, double * toz)
  {
    if (!supported_()) return false;

    // --- local trees
    double wa[3], wb[3];
    for (int d = 0; d < 3; d++) {
      wa[d] = n_[d]*a_.h[d]; wb[d] = n_[d]*b_.h[d];
      wa_[d] = wa[d]; wb_[d] = wb[d];
    }
    ta_.init(K_, a_.lo, wa);
    tb_.init(K_, b_.lo, wb);
    const double t0 = wall_seconds();
    bsrc_.build(src_, tb_);
    btgt_.build(tgt_, ta_);
    build_multipoles_();

    // --- kernel lookup for cell-cell near field
    setup_lookup_();
    setup_m2l_cache_(wa, wb);

    // --- accumulators (before scaling by grav_const)
    ncell_ = (size_t)n_[0]*n_[1]*n_[2];
    cacc_.assign(3*ncell_, 0.0);
    pacc_.assign(3*(size_t)tgt_.n, 0.0);

    const double t1 = wall_seconds();
    walk_(0,0,0,0, 0,0,0,0);
    const double t2 = wall_seconds();
    push_down_();
    const double t3 = wall_seconds();
    evaluate_leaves_();
    const double t4 = wall_seconds();
    Stats & st = stats();
    st.build += t1 - t0; st.walk += t2 - t1; st.down += t3 - t2; st.eval += t4 - t3;
    st.m2l += num_m2l; st.p2p_leaf += num_p2p_leaf; st.calls++;

    // --- write out
    const double G = p_.grav_const;
    size_t c = 0;
    for (int k = 0; k < n_[2]; k++)
    for (int j = 0; j < n_[1]; j++)
    for (int i = 0; i < n_[0]; i++, c++) {
      const size_t ia = (i + a_.gx) + a_.mx*((size_t)(j + a_.gy) + (size_t)a_.my*(k + a_.gz));
      ax[ia] += G*cacc_[3*c]; ay[ia] += G*cacc_[3*c+1]; az[ia] += G*cacc_[3*c+2];
    }
    for (int t = 0; t < tgt_.n; t++) {
      tox[t] += G*pacc_[3*t]; toy[t] += G*pacc_[3*t+1]; toz[t] += G*pacc_[3*t+2];
    }
    return true;
  }

  const char * why () const { return why_; }

  // statistics (diagnostics)
  long long num_m2l = 0, num_p2p_leaf = 0;

private:

  //-------------------------------------------------- setup

  /// records why the tree cannot be used for this pair of Blocks
  bool unsupported_ (const char * why) { why_ = why; return false; }

  bool supported_ () {
    if (p_.rank != 3) return unsupported_("the run is not 3D");
    n_[0] = a_.nx; n_[1] = a_.ny; n_[2] = a_.nz;
    if (b_.nx != n_[0] || b_.ny != n_[1] || b_.nz != n_[2])
      return unsupported_("the two Blocks do not have the same number of cells");
    if (n_[0] != n_[1] || n_[1] != n_[2])
      return unsupported_("the Blocks are not cubic");
    const int c = o_.chunk;
    if (c < 1 || c > 8 || n_[0] % c != 0)
      return unsupported_("leaf_chunk is not in 1..8 or does not divide the number of cells per Block");
    int g = n_[0] / c;
    if (g < 2 || (g & (g-1)) != 0 || g > 32)   // power of two, 2 to 6 levels
      return unsupported_("the number of cells per Block divided by leaf_chunk is not a power "
                          "of two between 2 and 32");
    K_ = 0; while ((1 << K_) < g) K_++;

    mixed_ = !same_spacing(a_, b_);
    if (mixed_) {
      // blocks of different resolution: cell widths differ by a power of two
      // and all cell centers lie on a lattice of half the fine cell width
      mg_ = mixed_geometry<Ew>(p_, a_, b_);
      if (!mg_.ok)
        return unsupported_("Blocks of different resolution do not have cell widths "
                            "related by a power of two on a common lattice");
      Ra_ = mg_.a_fine ? 1 : mg_.R;
      Rb_ = mg_.a_fine ? mg_.R : 1;
      for (int d = 0; d < 3; d++) wrap_[d] = mg_.wrap[d];
      return true;
    }

    // blocks must sit on a common lattice of cells
    for (int d = 0; d < 3; d++) {
      const double t = (b_.lo[d] - a_.lo[d]) / a_.h[d];
      O_[d] = std::lround(t);
      if (std::fabs(t - O_[d]) > 1e-6)
        return unsupported_("the Blocks are not on a common lattice of cells");
      if (Ew::enabled && p_.periodic[d]) {      // periodic images on the cell lattice
        const double nn = p_.dom_width[d] / a_.h[d];
        if (std::fabs(nn - std::lround(nn)) > 1e-6 || std::lround(nn) < 1)
          return unsupported_("the periodic domain is not a whole number of cells wide");
      }
    }
    return true;
  }

  /// width of the nodes of A at level la and of B at level lb, and the
  /// lattice unit for the displacement between their centers
  void m2l_unit_ (int la, int lb, const double * wa, const double * wb, double * un) const {
    for (int d = 0; d < 3; d++)
      un[d] = 0.5*std::min(wa[d]/double(1 << la), wb[d]/double(1 << lb));
  }

  void setup_m2l_cache_ (const double * wa, const double * wb) {
    m2l_cache_ok_ = false;
    if (!(Ew::enabled && p_.periodic[0] && p_.periodic[1] && p_.periodic[2])) return;
    for (int la = 0; la <= K_; la++) for (int lb = 0; lb <= K_; lb++) {
      double un[3];
      m2l_unit_(la, lb, wa, wb, un);
      m2l_ok_[la][lb] = true;
      for (int d = 0; d < 3; d++) {
        const double nn = p_.dom_width[d] / un[d];
        Ne_m2l_[la][lb][d] = std::lround(nn);
        if (std::fabs(nn - Ne_m2l_[la][lb][d]) > 1e-6 || Ne_m2l_[la][lb][d] < 1 ||
            Ne_m2l_[la][lb][d] >= (1L << 17)) m2l_ok_[la][lb] = false;
      }
    }
    // the cache is shared by all pairs; the physical scale and the resolution
    // ratio are part of the key
    m2l_scale_ = 0;
    for (int d = 0; d < 3; d++) {
      const double wf = mixed_ ? (mg_.a_fine ? wa[d] : wb[d]) : wa[d];
      const double dom = p_.dom_width[d];
      unsigned long long b1, b2;
      std::memcpy(&b1, &wf, 8); std::memcpy(&b2, &dom, 8);
      m2l_scale_ = (m2l_scale_ ^ b1 ^ (b2 * 0x9E3779B97F4A7C15ULL)) * 0xFF51AFD7ED558CCDULL;
    }
    m2l_ratio_code_ = 0;
    if (mixed_) {
      int l2 = 0; while ((1L << l2) < mg_.R) l2++;
      m2l_ratio_code_ = (mg_.a_fine ? 0 : 16) + l2;
    }
    M2LCache & mc = m2l_cache();
    if (mc.id != ew_.id() || mc.map.size() > (size_t(1) << 22)) {
      mc.map.clear(); mc.id = ew_.id();
    }
    m2l_cache_ok_ = true;
  }

  void setup_lookup_ () {
    if (mixed_) {
      MixedKernelCache & mk = mixed_kernel_cache();
      if (!mk.matches(ew_.id(), p_, mg_.u)) mk.reset(ew_.id(), p_, mg_.u);
      return;
    }
    for (int d = 0; d < 3; d++) {
      wrap_[d] = Ew::enabled && p_.periodic[d];
      Ne_[d] = 0;
      if (wrap_[d]) Ne_[d] = std::lround(p_.dom_width[d] / a_.h[d]);
    }
    NearKernelTable & nt = near_table();
    if (!nt.matches(ew_.id(), p_, a_.h)) nt.reset(ew_.id(), p_, a_.h);
  }

  /// kernel for a source cell at wrapped integer offset w (in units of half
  /// the fine cell width) from the target cell; blocks of different resolution
  void kernel_u_ (long wx, long wy, long wz, double * k) {
    MixedKernelCache & mk = mixed_kernel_cache();
    if (std::labs(wx) >= (1L << 20) || std::labs(wy) >= (1L << 20) || std::labs(wz) >= (1L << 20)) {
      // too large for the key (never the case in practice): do not cache
      const double d[3] = {wx*mg_.u[0], wy*mg_.u[1], wz*mg_.u[2]};
      kernel(p_, ew_, d, k);
      return;
    }
    const unsigned long long key = ((unsigned long long)(wx + (1L << 20)) << 42)
      | ((unsigned long long)(wy + (1L << 20)) << 21) | (unsigned long long)(wz + (1L << 20));
    auto f = mk.map.find(key);
    if (f == mk.map.end()) {
      const double d[3] = {wx*mg_.u[0], wy*mg_.u[1], wz*mg_.u[2]};
      std::array<double,3> kv;
      kernel(p_, ew_, d, kv.data());
      f = mk.map.emplace(key, kv).first;
    }
    k[0] = f->second[0]; k[1] = f->second[1]; k[2] = f->second[2];
  }

  /// integer cell offset o along axis d, wrapped to the nearest periodic image
  /// (same tie-breaking as wrapped_disp) if that axis is periodic
  long wrap_offset_ (long o, int d) const {
    if (!wrap_[d]) return o;
    long r = o % Ne_[d]; if (r < 0) r += Ne_[d];
    return (2*r < Ne_[d]) ? r : r - Ne_[d];
  }

  /// kernel (per unit mass) for a source cell at wrapped integer offset w
  void kernel_w_ (long wx, long wy, long wz, double * k) {
    const long R = NearKernelTable::R, S = NearKernelTable::S;
    if (std::labs(wx) <= R && std::labs(wy) <= R && std::labs(wz) <= R) {
      NearKernelTable & nt = near_table();
      const size_t g = (wx + R) + (size_t)S*((wy + R) + (size_t)S*(wz + R));
      if (!nt.done[g]) {
        const double d[3] = {wx*a_.h[0], wy*a_.h[1], wz*a_.h[2]};
        double kv[3];
        kernel(p_, ew_, d, kv);
        nt.K[3*g] = kv[0]; nt.K[3*g+1] = kv[1]; nt.K[3*g+2] = kv[2];
        nt.done[g] = 1;
      }
      k[0] = nt.K[3*g]; k[1] = nt.K[3*g+1]; k[2] = nt.K[3*g+2];
    } else {
      const double d[3] = {wx*a_.h[0], wy*a_.h[1], wz*a_.h[2]};
      kernel(p_, ew_, d, k);
    }
  }

  //-------------------------------------------------- multipoles

  void build_multipoles_ () {
    const int g = 1 << K_, c = o_.chunk;
    const double vol = b_.cell_volume();

    // leaves: moments about the leaf center from cells and particles
    for (int lz = 0; lz < g; lz++)
    for (int ly = 0; ly < g; ly++)
    for (int lx = 0; lx < g; lx++) {
      TNode & nd = tb_.nodes[tb_.id(K_, lx, ly, lz)];
      double ctr[3];
      tb_.center(K_, lx, ly, lz, ctr);
      for (int k = 0; k < c; k++)
      for (int j = 0; j < c; j++) {
        const double * dens = b_.row(c*ly + j, c*lz + k);
        for (int i = 0; i < c; i++) {
          const double cm = dens[c*lx + i]*vol;
          const double u[3] = {b_.pos(0, c*lx + i) - ctr[0], b_.pos(1, c*ly + j) - ctr[1],
                               b_.pos(2, c*lz + k) - ctr[2]};
          add_point_(cm, u, nd);
        }
      }
      const size_t lid = lx + g*(ly + (size_t)g*lz);
      for (int e = bsrc_.start[lid]; e < bsrc_.start[lid+1]; e++) {
        const int pi = bsrc_.items[e];
        const double u[3] = {src_.x[pi] - ctr[0], src_.y[pi] - ctr[1], src_.z[pi] - ctr[2]};
        add_point_(src_.m[pi], u, nd);
      }
    }

    // upward pass: shift the moments of the children to the parent's center
    const Idx & I = idx();
    for (int k = K_-1; k >= 0; k--) {
      const int gk = 1 << k;
      for (int lz = 0; lz < gk; lz++)
      for (int ly = 0; ly < gk; ly++)
      for (int lx = 0; lx < gk; lx++) {
        TNode & par = tb_.nodes[tb_.id(k, lx, ly, lz)];
        double cp[3];
        tb_.center(k, lx, ly, lz, cp);
        for (int ch = 0; ch < 8; ch++) {
          const int ci = 2*lx + (ch&1), cj = 2*ly + ((ch>>1)&1), ck = 2*lz + (ch>>2);
          const TNode & c1 = tb_.nodes[tb_.id(k+1, ci, cj, ck)];
          if (c1.m == 0) continue;   // (densities, hence masses, are not negative)
          double cc[3], dl[3];
          tb_.center(k+1, ci, cj, ck, cc);
          for (int d = 0; d < 3; d++) dl[d] = cc[d] - cp[d];
          par.m += c1.m;
          for (int i = 0; i < 3; i++) par.P[i] += c1.P[i] + c1.m*dl[i];
          for (int i = 0; i < 3; i++) for (int j = i; j < 3; j++)
            par.Q[I.i2[i][j]] += c1.Q[I.i2[i][j]] + c1.P[i]*dl[j] + dl[i]*c1.P[j]
                                 + c1.m*dl[i]*dl[j];
        }
      }
    }
  }

  /// add a point mass at offset u from the node center to its moments
  static void add_point_ (double m, const double * u, TNode & nd) {
    nd.m += m;
    nd.P[0] += m*u[0]; nd.P[1] += m*u[1]; nd.P[2] += m*u[2];
    nd.Q[0] += m*u[0]*u[0]; nd.Q[1] += m*u[0]*u[1]; nd.Q[2] += m*u[0]*u[2];
    nd.Q[3] += m*u[1]*u[1]; nd.Q[4] += m*u[1]*u[2]; nd.Q[5] += m*u[2]*u[2];
  }

  //-------------------------------------------------- dual tree walk

  void walk_ (int la, int ia, int ja, int ka, int lb, int ib, int jb, int kb)
  {
    const TNode & B = tb_.nodes[tb_.id(lb, ib, jb, kb)];
    if (B.m == 0) return;   // (densities, hence masses, are not negative)

    double ca[3], cb[3], r[3];
    ta_.center(la, ia, ja, ka, ca);
    tb_.center(lb, ib, jb, kb, cb);
    displacement<Ew::enabled>(p_, ca, cb, r);
    const double dist = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    const double ra = ta_.radius(la), rb = tb_.radius(lb);

    // well separated, and not (even partly) inside the softening length
    if ((ra + rb) < o_.theta*dist && dist > ra + rb + hsoft_) {
      m2l_(la, lb, ta_.nodes[ta_.id(la, ia, ja, ka)], B, r);
      num_m2l++;
      return;
    }

    if (la == K_ && lb == K_) {
      p2p_leaf_(ia, ja, ka, ib, jb, kb);
      num_p2p_leaf++;
      return;
    }

    const bool split_a = (lb == K_) || (la < K_ && ra >= rb);
    if (split_a) {
      for (int ch = 0; ch < 8; ch++)
        walk_(la+1, 2*ia + (ch&1), 2*ja + ((ch>>1)&1), 2*ka + (ch>>2), lb, ib, jb, kb);
    } else {
      for (int ch = 0; ch < 8; ch++)
        walk_(la, ia, ja, ka, lb+1, 2*ib + (ch&1), 2*jb + ((ch>>1)&1), 2*kb + (ch>>2));
    }
  }

  /// tensors of the Green's function (Newton and, if periodic, the Ewald
  /// correction) at displacement r
  void compute_tensors_ (const double * r, M2LTensors & T) const {
    newton_tensors(r, T.d1, T.d2, T.d3);
    if (Ew::enabled) {
      // The Ewald tables (EnzoEwald) hold minus the derivative tensors of the
      // periodic correction (they are in the sign convention of the
      // acceleration, which is what the direct sums add), whereas d1,d2,d3
      // here are the tensors of 1/r, so the tables are subtracted.
      const std::array<double,3> e1 = ew_(r[0], r[1], r[2]);
      const std::array<double,6> e2 = ew_.d2(r[0], r[1], r[2]);
      const std::array<double,10> e3 = ew_.d3(r[0], r[1], r[2]);
      for (int i = 0; i < 3; i++) T.d1[i] -= e1[i];
      for (int i = 0; i < 6; i++) T.d2[i] -= e2[i];
      for (int i = 0; i < 10; i++) T.d3[i] -= e3[i];
    }
  }

  /// For periodic problems the displacement between the centers of two nodes
  /// is a lattice vector (the unit is half the smaller node width), and there
  /// are only a few hundred distinct ones per pair of levels, so the tensors
  /// are cached by (levels, resolution ratio, wrapped offset).
  const M2LTensors & tensors_ (int la, int lb, const double * r, M2LTensors & scratch) {
    if (m2l_cache_ok_ && m2l_ok_[la][lb]) {
      double un[3]; long o[3], g[3];
      m2l_unit_(la, lb, wa_, wb_, un);
      bool lattice = true;
      for (int d = 0; d < 3; d++) {
        const double t = r[d]/un[d];
        o[d] = std::lround(t);
        if (std::fabs(t - o[d]) > 1e-6) lattice = false;
        g[d] = wrapped_index(o[d], Ne_m2l_[la][lb][d]);
      }
      if (lattice) {
        M2LKey key;
        key.scale = m2l_scale_;
        key.pk = ((((((unsigned long long)(la << 4 | lb) << 5 | (unsigned long long)m2l_ratio_code_) << 17
                   | (unsigned long long)g[0]) << 17) | (unsigned long long)g[1]) << 17)
                 | (unsigned long long)g[2];
        M2LCache & mc = m2l_cache();
        auto f = mc.map.find(key);
        if (f != mc.map.end()) return f->second;
        const double d[3] = {(g[0] - Ne_m2l_[la][lb][0]/2)*un[0],
                             (g[1] - Ne_m2l_[la][lb][1]/2)*un[1],
                             (g[2] - Ne_m2l_[la][lb][2]/2)*un[2]};
        M2LTensors T;
        compute_tensors_(d, T);
        return mc.map.emplace(key, T).first->second;
      }
    }
    compute_tensors_(r, scratch);
    return scratch;
  }

  void m2l_ (int la, int lb, TNode & A, const TNode & B, const double * r) {
    M2LTensors scratch;
    const M2LTensors & T = tensors_(la, lb, r, scratch);
    const Idx & I = idx();
    // local coefficients from the source moments about its center:
    //   c1 = M d1 + P.d2 + 1/2 Q:d3,   c2 = M d2 + P.d3,   c3 = M d3
    for (int i = 0; i < 3; i++) {
      double s = B.m*T.d1[i];
      for (int j = 0; j < 3; j++) {
        s += B.P[j]*T.d2[I.i2[i][j]];
        for (int k = 0; k < 3; k++) s += 0.5*B.Q[I.i2[j][k]]*T.d3[I.i3[i][j][k]];
      }
      A.c1[i] += s;
    }
    for (int i = 0; i < 3; i++) for (int j = i; j < 3; j++) {
      double s = B.m*T.d2[I.i2[i][j]];
      for (int k = 0; k < 3; k++) s += B.P[k]*T.d3[I.i3[i][j][k]];
      A.c2[I.i2[i][j]] += s;
    }
    for (int i = 0; i < 10; i++) A.c3[i] += B.m*T.d3[i];
  }

  /// direct interaction between leaf chunk (ia,ja,ka) of A and (ib,jb,kb) of B
  void p2p_leaf_ (int ia, int ja, int ka, int ib, int jb, int kb)
  {
    const int c = o_.chunk, g = 1 << K_;
    const double vol_b = b_.cell_volume();
    const bool same_leaf = self_ && ia == ib && ja == jb && ka == kb;
    const bool tm = o_.timing;
    double tc0 = tm ? wall_seconds() : 0;

    if (mixed_) {
      // blocks of different resolution: the offset between a cell of A and a
      // cell of B, in units of half the fine cell width, is
      //   D = LO + (2 j + 1) Rb - (2 i + 1) Ra
      long W[3][8][8];
      const int la3[3] = {ia, ja, ka}, lb3[3] = {ib, jb, kb};
      for (int d = 0; d < 3; d++)
        for (int i = 0; i < c; i++) for (int j = 0; j < c; j++) {
          long Dv = mg_.LO[d] + (2*(c*lb3[d] + j) + 1)*Rb_ - (2*(c*la3[d] + i) + 1)*Ra_;
          if (mg_.wrap[d]) Dv = wrap_units(Dv, mg_.Ne[d]);
          W[d][i][j] = Dv;
        }
      for (int k = 0; k < c; k++) for (int j = 0; j < c; j++) for (int i = 0; i < c; i++) {
        const int ci = c*ia + i, cj = c*ja + j, ck = c*ka + k;
        double sx = 0, sy = 0, sz = 0;
        for (int k2 = 0; k2 < c; k2++) for (int j2 = 0; j2 < c; j2++) {
          const double * dens = b_.row(c*jb + j2, c*kb + k2) + c*ib;
          for (int i2 = 0; i2 < c; i2++) {
            double kv[3];
            kernel_u_(W[0][i][i2], W[1][j][j2], W[2][k][k2], kv);
            const double m = dens[i2]*vol_b;
            sx += m*kv[0]; sy += m*kv[1]; sz += m*kv[2];
          }
        }
        const size_t cid = ci + (size_t)n_[0]*(cj + (size_t)n_[1]*ck);
        cacc_[3*cid] += sx; cacc_[3*cid+1] += sy; cacc_[3*cid+2] += sz;
      }
    } else {
    // cells <- cells: a pair of cells is characterized by the integer offset
    // between them; there are only (2c-1)^3 distinct offsets for the c^3 x c^3
    // pairs of this leaf pair, so evaluate the kernel once per offset
    const int D = 2*c - 1;
    kloc_.resize(3*(size_t)D*D*D);
    {
      const long bx = O_[0] + (long)c*(ib - ia), by = O_[1] + (long)c*(jb - ja),
                 bz = O_[2] + (long)c*(kb - ka);
      long wx[16], wy[16], wz[16];
      for (int o = 0; o < D; o++) {
        wx[o] = wrap_offset_(bx + o - (c-1), 0);
        wy[o] = wrap_offset_(by + o - (c-1), 1);
        wz[o] = wrap_offset_(bz + o - (c-1), 2);
      }
      size_t e = 0;
      for (int oz = 0; oz < D; oz++) for (int oy = 0; oy < D; oy++) for (int ox = 0; ox < D; ox++, e++)
        kernel_w_(wx[ox], wy[oy], wz[oz], &kloc_[3*e]);
    }
    for (int k = 0; k < c; k++) for (int j = 0; j < c; j++) for (int i = 0; i < c; i++) {
      const int ci = c*ia + i, cj = c*ja + j, ck = c*ka + k;
      double sx = 0, sy = 0, sz = 0;
      for (int k2 = 0; k2 < c; k2++) for (int j2 = 0; j2 < c; j2++) {
        const double * dens = b_.row(c*jb + j2, c*kb + k2) + c*ib;
        const size_t row = (size_t)D*((j2 - j + c-1) + (size_t)D*(k2 - k + c-1));
        for (int i2 = 0; i2 < c; i2++) {
          if (same_leaf && i2 == i && j2 == j && k2 == k) continue;
          const double * kv = &kloc_[3*(row + (i2 - i + c-1))];
          const double m = dens[i2]*vol_b;
          sx += m*kv[0]; sy += m*kv[1]; sz += m*kv[2];
        }
      }
      const size_t cid = ci + (size_t)n_[0]*(cj + (size_t)n_[1]*ck);
      cacc_[3*cid] += sx; cacc_[3*cid+1] += sy; cacc_[3*cid+2] += sz;
    }

    }

    double tc1 = tm ? wall_seconds() : 0;
    if (tm) stats().p2p_cc += tc1 - tc0;
    const size_t lidA = ia + g*(ja + (size_t)g*ka), lidB = ib + g*(jb + (size_t)g*kb);
    const int s0 = bsrc_.start[lidB], s1 = bsrc_.start[lidB+1];
    const int t0 = btgt_.start[lidA], t1 = btgt_.start[lidA+1];

    // cells(A) <- particles(B)
    if (s1 > s0) {
      for (int k = 0; k < c; k++) for (int j = 0; j < c; j++) for (int i = 0; i < c; i++) {
        const int ci = c*ia + i, cj = c*ja + j, ck = c*ka + k;
        const double pa[3] = {a_.pos(0,ci), a_.pos(1,cj), a_.pos(2,ck)};
        double sx = 0, sy = 0, sz = 0;
        for (int e = s0; e < s1; e++) {
          const int q = bsrc_.items[e];
          const double pb[3] = {src_.x[q], src_.y[q], src_.z[q]};
          double d[3], kv[3];
          displacement<Ew::enabled>(p_, pa, pb, d);
          kernel(p_, ew_, d, kv);
          sx += src_.m[q]*kv[0]; sy += src_.m[q]*kv[1]; sz += src_.m[q]*kv[2];
        }
        const size_t cid = ci + (size_t)n_[0]*(cj + (size_t)n_[1]*ck);
        cacc_[3*cid] += sx; cacc_[3*cid+1] += sy; cacc_[3*cid+2] += sz;
      }
    }

    if (tm) stats().p2p_cp += wall_seconds() - tc1;

    // particles(A) <- cells(B), particles(B)
    double t_pc = 0, t_pp = 0;
    for (int e = t0; e < t1; e++) {
      double tp0 = tm ? wall_seconds() : 0;
      const int t = btgt_.items[e];
      const double pa[3] = {tgt_.x[t], tgt_.y[t], tgt_.z[t]};
      double sx = 0, sy = 0, sz = 0;
      for (int k2 = 0; k2 < c; k2++) for (int j2 = 0; j2 < c; j2++) {
        const double * dens = b_.row(c*jb + j2, c*kb + k2);
        for (int i2 = 0; i2 < c; i2++) {
          const double pb[3] = {b_.pos(0, c*ib + i2), b_.pos(1, c*jb + j2), b_.pos(2, c*kb + k2)};
          double d[3], kv[3];
          displacement<Ew::enabled>(p_, pa, pb, d);
          kernel(p_, ew_, d, kv);
          const double m = dens[c*ib + i2]*vol_b;
          sx += m*kv[0]; sy += m*kv[1]; sz += m*kv[2];
        }
      }
      double tp1 = tm ? wall_seconds() : 0;
      t_pc += tp1 - tp0;
      const int skip = (self_ && exclude_) ? exclude_[t] : -1;
      for (int f = s0; f < s1; f++) {
        const int q = bsrc_.items[f];
        if (q == skip) continue;
        const double pb[3] = {src_.x[q], src_.y[q], src_.z[q]};
        double d[3], kv[3];
        displacement<Ew::enabled>(p_, pa, pb, d);
        kernel(p_, ew_, d, kv);
        sx += src_.m[q]*kv[0]; sy += src_.m[q]*kv[1]; sz += src_.m[q]*kv[2];
      }
      pacc_[3*t] += sx; pacc_[3*t+1] += sy; pacc_[3*t+2] += sz;
      if (tm) t_pp += wall_seconds() - tp1;
    }
    if (tm) { stats().p2p_pc += t_pc; stats().p2p_pp += t_pp; }
  }

  //-------------------------------------------------- downward pass

  /// shift the coefficients of `par` (center cp) to a child (center cc) and add
  static void l2l_ (const TNode & par, const double * cp, TNode & ch, const double * cc) {
    const Idx & I = idx();
    double dl[3];
    for (int d = 0; d < 3; d++) dl[d] = cc[d] - cp[d];
    for (int i = 0; i < 3; i++) {
      double s = par.c1[i];
      for (int j = 0; j < 3; j++) {
        s -= par.c2[I.i2[i][j]]*dl[j];
        double t = 0;
        for (int k = 0; k < 3; k++) t += par.c3[I.i3[i][j][k]]*dl[k];
        s += 0.5*t*dl[j];
      }
      ch.c1[i] += s;
    }
    for (int i = 0; i < 3; i++) for (int j = i; j < 3; j++) {
      double s = par.c2[I.i2[i][j]];
      for (int k = 0; k < 3; k++) s -= par.c3[I.i3[i][j][k]]*dl[k];
      ch.c2[I.i2[i][j]] += s;
    }
    for (int i = 0; i < 10; i++) ch.c3[i] += par.c3[i];
  }

  void push_down_ () {
    for (int k = 0; k < K_; k++) {
      const int gk = 1 << k;
      for (int lz = 0; lz < gk; lz++)
      for (int ly = 0; ly < gk; ly++)
      for (int lx = 0; lx < gk; lx++) {
        const TNode & par = ta_.nodes[ta_.id(k, lx, ly, lz)];
        double cp[3], cc[3];
        ta_.center(k, lx, ly, lz, cp);
        for (int ch = 0; ch < 8; ch++) {
          const int ci = 2*lx + (ch&1), cj = 2*ly + ((ch>>1)&1), ck = 2*lz + (ch>>2);
          ta_.center(k+1, ci, cj, ck, cc);
          l2l_(par, cp, ta_.nodes[ta_.id(k+1, ci, cj, ck)], cc);
        }
      }
    }
  }

  static void eval_local_ (const TNode & n, const double * a, double * acc) {
    const Idx & I = idx();
    for (int i = 0; i < 3; i++) {
      double s = n.c1[i];
      for (int j = 0; j < 3; j++) {
        s -= n.c2[I.i2[i][j]]*a[j];
        double t = 0;
        for (int k = 0; k < 3; k++) t += n.c3[I.i3[i][j][k]]*a[k];
        s += 0.5*t*a[j];
      }
      acc[i] = -s;
    }
  }

  void evaluate_leaves_ () {
    const int g = 1 << K_, c = o_.chunk;
    for (int lz = 0; lz < g; lz++)
    for (int ly = 0; ly < g; ly++)
    for (int lx = 0; lx < g; lx++) {
      const TNode & nd = ta_.nodes[ta_.id(K_, lx, ly, lz)];
      double ctr[3];
      ta_.center(K_, lx, ly, lz, ctr);
      for (int k = 0; k < c; k++) for (int j = 0; j < c; j++) for (int i = 0; i < c; i++) {
        const int ci = c*lx + i, cj = c*ly + j, ck = c*lz + k;
        const double a[3] = {a_.pos(0,ci) - ctr[0], a_.pos(1,cj) - ctr[1], a_.pos(2,ck) - ctr[2]};
        double acc[3];
        eval_local_(nd, a, acc);
        const size_t cid = ci + (size_t)n_[0]*(cj + (size_t)n_[1]*ck);
        cacc_[3*cid] += acc[0]; cacc_[3*cid+1] += acc[1]; cacc_[3*cid+2] += acc[2];
      }
      const size_t lid = lx + g*(ly + (size_t)g*lz);
      for (int e = btgt_.start[lid]; e < btgt_.start[lid+1]; e++) {
        const int t = btgt_.items[e];
        const double a[3] = {tgt_.x[t] - ctr[0], tgt_.y[t] - ctr[1], tgt_.z[t] - ctr[2]};
        double acc[3];
        eval_local_(nd, a, acc);
        pacc_[3*t] += acc[0]; pacc_[3*t+1] += acc[1]; pacc_[3*t+2] += acc[2];
      }
    }
  }

  //-------------------------------------------------- data

  const Params & p_;
  const Ew & ew_;
  SubtreeOpts o_;
  const CellGrid & a_;
  const CellGrid & b_;
  bool self_;
  const char * why_ = nullptr;
  PointSet tgt_, src_;
  const int * exclude_;
  double hsoft_;

  int n_[3] = {0,0,0};
  int K_ = 0;
  long O_[3] = {0,0,0};
  bool wrap_[3] = {false,false,false};
  bool m2l_cache_ok_ = false;
  bool m2l_ok_[16][16];
  long Ne_m2l_[16][16][3];
  unsigned long long m2l_scale_ = 0;
  int m2l_ratio_code_ = 0;
  double wa_[3] = {0,0,0}, wb_[3] = {0,0,0};     // block extents of A and B
  bool mixed_ = false;                           // blocks of different resolution
  MixedGeom mg_;
  long Ra_ = 1, Rb_ = 1;
  long Ne_[3] = {0,0,0};
  size_t ncell_ = 0;

  Tree ta_, tb_;
  Buckets bsrc_, btgt_;
  std::vector<double> cacc_, pacc_, kloc_;
};

} // namespace subtree

/// Tree-based version of accumulate_block_pair_direct (same arguments plus
/// options). Returns false, leaving the outputs untouched, if the pair of
/// Blocks does not have the structure the tree needs (equal cubic Blocks of
/// chunk*2^K cells on a common lattice); use the direct sums then. If `why`
/// is given it is then set to the reason.
template <class Ew, class real>
bool accumulate_block_pair_tree
(const Params & p, const Ew & ew, const SubtreeOpts & opts,
 const CellGrid & a, const CellGrid & b, bool self,
 real * ax, real * ay, real * az,
 const PointSet & tgt, double * tox, double * toy, double * toz,
 const PointSet & src, const int * exclude, const char ** why = nullptr)
{
  subtree::Solver<Ew, real> solver(p, ew, opts, a, b, self, tgt, src, exclude);
  const bool ok = solver.run(ax, ay, az, tox, toy, toz);
  if (!ok && why) *why = solver.why();
  return ok;
}

} // namespace fmm_direct

#endif /* ENZO_METHOD_MULTIPOLE_SUBTREE_HPP */

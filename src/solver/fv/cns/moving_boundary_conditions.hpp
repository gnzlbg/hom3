#ifndef HOM3_SOLVERS_FV_CNS_MOVING_BOUNDARY_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_CNS_MOVING_BOUNDARY_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the CNS moving boundary conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/boundary_condition.hpp"
#include "interpolation/rbf.hpp"
#include "geometry/algorithms.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace cns {
////////////////////////////////////////////////////////////////////////////////
namespace bc {
////////////////////////////////////////////////////////////////////////////////

#define NAIVE
/// \brief Moving wall boundary conditions
namespace moving_wall {

/// \brief Dirichlet boundary condition for the temperature
template<class Solver>
struct AdiabaticNoSlip : fv::bc::Condition<AdiabaticNoSlip<Solver>> {
  static const SInd nd = Solver::nd;
  static const SInd nvars = Solver::nvars;
  using V = Indices<nd>;

  template<class G, class V, class A>
      AdiabaticNoSlip(Solver& solver, G&& geometry, V&& vel, A&& acc) noexcept
    : s(solver)
    , g([=](const NumA<Solver::nd> x) { return (*geometry)(x); })
    , velocity(vel), acceleration(acc) {}


  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, const CellIdx ghostIdx) const noexcept {
    NumA<nvars> gCellVars;  /// PV of ghost cell:

    auto u_d = [&](const CellIdx i, const SInd d) {
      const auto pvars = s.pv(_(), i);
      return pvars(V::u(d));
    };

    for (auto d : s.grid().dimensions()) {  /// velocity_srfc = 0
      gCellVars(V::u(d)) = this->mb_dirichlet([&](CellIdx i) { return u_d(i, d); }, s,
                                              g, ghostIdx, [&](CellIdx) {
          return velocity(d);
      });
    }

    /// pressure gradient normal to surface = 0
    auto p = [&](const CellIdx i) {
      const auto pvars = s.pv(_(), i);
      return pvars(V::p());
    };
    Num p_n_srfc = 0.;
    Num p_srfc = 0.;
    Num rho_srfc = 0.;
    {  // compute p_n_srfc
#ifdef NAIVE
      std::vector<CellIdx> nghbrs; nghbrs.reserve(6);
      for (auto nghbr_pos : s.grid().neighbor_positions()) {
        const auto nghbrIdx = s.cells().neighbors(ghostIdx, nghbr_pos);
        if (is_valid(nghbrIdx) && !s.is_ghost_cell(nghbrIdx)
            && !is_valid(s.cut_by_which_moving_boundary(nghbrIdx))) {
          nghbrs.push_back(nghbrIdx);
        }
      }

      auto cs = s.grid().cut_surface(s.node_idx(ghostIdx), g);
      auto n = cs.normal;

      /// normal velocity
      NumA<nd> u_body = NumA<nd>::Constant(0.);
      for(SInd d = 0; d < Solver::nd; ++d) {
        u_body(d) = velocity(d);
      }
      Num u_n = u_body.dot(n);

      /// normal acceleration
      NumA<nd> a_body = NumA<nd>::Constant(0.);
      for(SInd d = 0; d < Solver::nd; ++d) {
        a_body(d) = acceleration(d);
      }
      Num a_n = a_body.dot(n);

      Num dx = 0;
      Num u_mp_n = 0;
      Num rho = 0;
      Num p_bndry = 0;
      if (!(nghbrs.size() == 0)) {
        if (nghbrs.size() == 1) {
          auto bndryIdx = nghbrs[0];
          dx = s.cell_dx(bndryIdx, ghostIdx);
          NumA<nd> u_bndry = NumA<nd>::Constant(0.);
          for(SInd d = 0; d < Solver::nd; ++d) {
            u_bndry(d) = u_d(bndryIdx, d);
          }
          u_mp_n = u_bndry.dot(n);
          rho = s.pv(_(), bndryIdx)(V::rho());
          p_bndry = s.pv(_(), bndryIdx)(V::p());
        } else if (nghbrs.size() == 2) {
          auto bndryIdx0 = nghbrs[0];
          auto bndryIdx1 = nghbrs[1];
          auto dx0 = s.cell_dx(bndryIdx0, ghostIdx);
          auto dx1 = s.cell_dx(bndryIdx1, ghostIdx);
          NumA<nd> u_bndry0 = NumA<nd>::Constant(0.);
          for(SInd d = 0; d < Solver::nd; ++d) {
            u_bndry0(d) = u_d(bndryIdx0, d);
          }
          NumA<nd> u_bndry1 = NumA<nd>::Constant(0.);
          for(SInd d = 0; d < Solver::nd; ++d) {
            u_bndry1(d) = u_d(bndryIdx1, d);
          }
          auto u0 = u_bndry0.dot(n);
          auto u1 = u_bndry1.dot(n);
          u_mp_n = 0.5 * (u0 + u1);
          dx = 0.5 * (dx0 + dx1);
          auto rho0 = s.pv(_(), bndryIdx0)(V::rho());
          auto rho1 = s.pv(_(), bndryIdx1)(V::rho());
          rho = 0.5 * (rho0 + rho1);
          auto p0 = s.pv(_(), bndryIdx0)(V::p());
          auto p1 = s.pv(_(), bndryIdx1)(V::p());
          p_bndry = 0.5 * (p0 + p1);
        }

        /// (u_n)_n
        auto u_n_n = dx > 1e-5? (u_mp_n - u_n) / dx : 0.;
        rho_srfc = rho;
        p_srfc = p_bndry;
        Num p_n = - (rho * a_n + 2 * rho * u_n * u_n_n);
        p_n_srfc = p_n;
        // std::cerr << "gIdx: " << ghostIdx << "  p_n: " << p_n
        //           << " u_n: " << u_n << " us_mp_n: " << u_mp_n

        //           << " a_n: " << a_n << " u_n_n: " << u_n_n
        //           << " rho: " << rho << " rho_s: " << rho_srfc
        //           << " p_s: " << p_srfc
        //           << " rho_n: " <<  p_n_srfc * rho_srfc / p_srfc
        //           << " dx: " << dx << "\n";
      }
#endif
#ifdef ACCURATE
      bool has_direct_non_mb_nghbrs = false;
      for(auto npos : s.grid().neighbor_positions()) {
        auto id = s.cells().neighbors(ghostIdx, npos);
        if (is_valid(id) && !is_valid(s.cut_by_which_moving_boundary(id))) {
          has_direct_non_mb_nghbrs = true;
        }
      }

      if(has_direct_non_mb_nghbrs) {

      /// find normal vector and mirror point
      auto cs = s.grid().cut_surface(s.node_idx(ghostIdx), g);
      auto n = cs.normal;

      /// In plane surface vector
      auto srfc = Eigen::ParametrizedLine<Num, Solver::nd>::Through
                  (cs.cut_points[0], cs.cut_points[1]);

      /// Line perpendicular to surface crossing the cell center
      NumA<Solver::nd> x_gc = s.cells().x_center.row(ghostIdx);
      auto line = Eigen::ParametrizedLine<Num, Solver::nd>(x_gc, cs.normal);
      auto x_ip = line.intersectionPoint
                  (Eigen::Hyperplane<Num, Solver::nd>(srfc));

      /// Compute the mirror point
      auto sgc = math::signum(g(x_gc));
      auto dist = (x_ip - x_gc).norm();
      auto x_mp = sgc > 0.  // Is always inside the domain
                  ? line.pointAt(dist)  // note: line starts at gc
                  : line.pointAt(2. * dist);

      /// Create interpolation stencil for mirror point
      auto ip_nghbrs = s.find_interpolation_neighbors(ghostIdx);
      ASSERT(ip_nghbrs.size() != 0, "no nghbrs found?!");
      std::vector<CellIdx> nghbrs; nghbrs.reserve(10);
      for(auto in : ip_nghbrs) {  /// set of non-mb interpolation neighbors
        if(s.is_ghost_cell(in) && is_valid(s.cut_by_which_moving_boundary(in))) {
          continue;
        } else {
          nghbrs.push_back(in);
        }
      }
      ASSERT(nghbrs.size() > 0, "cell: " << ghostIdx
             << " has no neighbors for mb neumann condition");

      /// interpolation stencil for velocity
      std::vector<NumA<Solver::nd>> xs;
      std::vector<std::vector<Num>> u_vs;
      xs.reserve(nghbrs.size() + 3); u_vs.reserve(nghbrs.size() + 3);
      for (auto in : nghbrs) {
        xs.push_back(s.cells().x_center.row(in));
        u_vs.push_back(std::vector<Num>{ u_d(in, 0), u_d(in, 1) });
      }
      xs.push_back(x_ip);  // add surface vel to stencil
      u_vs.push_back(std::vector<Num>{velocity(0), velocity(1)});

      /// interpolate values to mirror point
      std::vector<Num> vs_mp = s.interpolate(x_mp, xs, u_vs);

      /// interpolate rho to x_gc
      xs.pop_back();   // remove surface from stencil
      std::vector<Num> rho_vs;
      for (auto in : nghbrs) {
        rho_vs.push_back(s.pv(_(), in)(V::rho()));
      }
      Num rho = s.interpolate(x_gc, xs, rho_vs);

      /// interpolate rho to x_srfc_centroid
      rho_srfc = s.interpolate(cs.centroid, xs, rho_vs);

      /// interpolate p to x_srfc_centroid
      SInd index = 0;
      for (auto in : nghbrs) {
        rho_vs[index++] = s.pv(_(), in)(V::p());
      }
      p_srfc = s.interpolate(cs.centroid, xs, rho_vs);

      /// vel vector at mirror point
      NumA<nd> us_mp = NumA<nd>::Constant(0.);
      for(SInd d = 0; d < Solver::nd; ++d) {
        us_mp(d) = vs_mp[d];
      }
      /// normal vel vector at mirror point
      auto us_mp_n = us_mp.dot(n);

      /// compute (u_n)_n
      NumA<nd> u_body = NumA<nd>::Constant(0.);
      for(SInd d = 0; d < Solver::nd; ++d) {
        u_body(d) = velocity(d);
      }
      Num u_n = u_body.dot(n);
      Num dx = (x_ip - x_mp).norm();
      auto u_n_n = dx > 1e-5? (us_mp_n - u_n) / dx : 0.;

      /// normal acceleration
      NumA<nd> a_body = NumA<nd>::Constant(0.);
      for(SInd d = 0; d < Solver::nd; ++d) {
        a_body(d) = acceleration(d);
      }
      Num a_n = a_body.dot(n);

      Num rho_t = 0.;  ///< Assume small density changes in time for now

      /// p_n - (rho_t * u_n + rho * a + 2 rho u_n u_n_n)
      //Num p_n = - (rho_t * u_n + rho * a_n + 2 * rho * u_n * u_n_n);
      Num p_n = - rho * a_n;
      p_n_srfc = p_n;
      std::cerr << "gIdx: " << ghostIdx << "  p_n: " << p_n
                << " u_n: " << u_n << " us_mp_n: " << us_mp_n
                << " us_mp: " << us_mp(0) << ", " << us_mp(1)
                << " a_n: " << a_n << " u_n_n: " << u_n_n
                << " rho: " << rho << " rho_s: " << rho_srfc
                << " p_s: " << p_srfc
                << " rho_n: " <<  p_n_srfc * rho_srfc / p_srfc
                << " dx: " << dx << "\n";
      }
      #endif
    }
    gCellVars(V::p()) = this->mb_neumann(p, s, g, ghostIdx, [&](CellIdx) { return p_n_srfc; });


    /// density gradient normal to surface = 0
    auto rho = [&](const CellIdx i) {
      const auto pvars = s.pv(_(), i);
      return pvars(V::rho());
    };
    Num rho_n = 0.;
    {
      rho_n = p_n_srfc * rho_srfc / p_srfc;
    }
    gCellVars(V::rho()) = this->mb_neumann(rho, s, g, ghostIdx, [&](CellIdx){ return rho_n;});


    s.Q(_(), ghostIdx) = s.cv(gCellVars);
  }

  bool is_moving() const noexcept { return true; }

  Solver& s;
  const std::function<Num(NumA<Solver::nd>)> g;
  const std::function<Num(SInd d)> velocity;
  const std::function<Num(SInd d)> acceleration;
};


}  // namespace moving_wall

////////////////////////////////////////////////////////////////////////////////
}  // namespace bc
}  // namespace cns
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif

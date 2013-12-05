#ifndef HOM3_SOLVERS_FV_BOUNDARY_CONDITION_HPP_
#define HOM3_SOLVERS_FV_BOUNDARY_CONDITION_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include "globals.hpp"
#include "grid/boundary.hpp"
#include "tags.hpp"
/// Options:
#define ENABLE_DBG_ 0
#include "misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv {
////////////////////////////////////////////////////////////////////////////////

namespace bc {

/// \brief Implements the Finite Volume boundary condition concept
///
/// Refines the grid::boundary concept.
template<SInd nd> struct Interface : grid::boundary::Interface<nd> {
  /// \brief Construct a boundary condition
  template<class Solver, class Geometry, class Condition>
  Interface(const String name, Geometry&& geometry,
            const Solver& solver, Condition&& condition) noexcept
      : grid::boundary::Interface<nd>(name, geometry, solver)
      , apply_lhs(condition), apply_rhs(condition)
      , slope_lhs([=](const CellIdx cIdx, const SInd v, const SInd dir) {
          return condition.template slope<lhs_tag>(cIdx, v, dir); })
      , slope_rhs([=](const CellIdx cIdx, const SInd v, const SInd dir) {
          return condition.template slope<rhs_tag>(cIdx, v, dir); })
      , is_moving([=]() -> bool { return condition.is_moving(); })
      , velocity([=](SInd d) -> Num { return condition.velocity(d); })
  {}

  /// Applies the boundary condition to the lhs of a GhostCellRange
  void apply(lhs_tag, const CellIdx ghostIdx) const noexcept
  { apply_lhs(lhs, ghostIdx); }

  /// Applies the boundary condition to the rhs of a GhostCellRange
  void apply(rhs_tag, const CellIdx ghostIdx) const noexcept
  { apply_rhs(rhs, ghostIdx); }

  Num slope(lhs_tag, const CellIdx cIdx, const SInd v,
            const SInd dir) const noexcept
  { return slope_lhs(cIdx, v, dir); }

  Num slope(rhs_tag, const CellIdx cIdx, const SInd v,
            const SInd dir) const noexcept
  { return slope_rhs(cIdx, v, dir); }

 private:
  const std::function<void(lhs_tag, const CellIdx)> apply_lhs;
  const std::function<void(rhs_tag, const CellIdx)> apply_rhs;
  const std::function<Num(CellIdx, SInd, SInd)> slope_lhs;
  const std::function<Num(CellIdx, SInd, SInd)> slope_rhs;
 public:
  const std::function<bool()> is_moving;
  const std::function<Num(SInd)> velocity;
};

/// \brief Contains general boundary condition functionality
template<class BC> struct Condition {
  /// \brief Imposes boundary cell's tangential slope on the boundary surface
  template<class _>
  Num slope(const CellIdx bndryIdx, const SInd v,
            const SInd dir) const noexcept {
    return static_cast<const BC*>(this)->s.template slope<_>(bndryIdx, v, dir);
  }

  /// \brief Ghost cell value for dirichlet boundary condition
  inline Num dirichlet(const Num boundaryCellValue,
                       const Num surfaceValue = 0) const noexcept
  { return 2.0 * surfaceValue - boundaryCellValue; }

  /// \brief Ghost cell value for dirichlet boundary condition
  inline Num neumann(const Num boundaryCellValue, const Num surfaceFlux = 0,
                     const Num cell_distance = 0) const noexcept
  { return boundaryCellValue - cell_distance * surfaceFlux; }

  inline bool is_moving() const noexcept { return false; }
  inline Num velocity(SInd) const noexcept { return 0.0; }

  template<class CV, class S, class G, class V>
  inline Num mb_dirichlet(CV&& cell_values, S&& s, G&& g, const CellIdx ghostIdx,
                          V&& value_srfc) const noexcept {
    /// cut surface:
    auto cs = s.grid().cut_surface(s.node_idx(ghostIdx), g);

    /// get interpolation neighbors
    auto nghbrs = s.find_interpolation_neighbors(ghostIdx);
    ASSERT(nghbrs.size() != 0, "no nghbrs found?!");
    std::size_t mbNghbrs = 0;
    for(auto n : nghbrs) {
      if(s.is_ghost_cell(n)
         && is_valid(s.cut_by_which_moving_boundary(n))) {
        ++mbNghbrs;
      }
    }

    if(!(nghbrs.size() - mbNghbrs > 1)) {
      if(nghbrs.size() != mbNghbrs) {
        /// only one neighbor is not a mb neighbor
        // note: this cell won't be part of any stencil
        ASSERT(nghbrs.size() > 0, "no neighbors found!!");
        CellIdx nghbr = invalid<CellIdx>();
        for(auto n : nghbrs) {
          if(s.is_ghost_cell(n) && is_valid(s.cut_by_which_moving_boundary(n))) {
            continue;
          }
          nghbr = n;
          break;
        }
        ASSERT(is_valid(nghbr), "WTF: ghostIdx: " << ghostIdx << " has: "
               << nghbrs.size() << " interpolation neighbors, of which "
               << mbNghbrs << " are mb ghost cells");
        auto x_n = NumA<std::decay_t<S>::nd>{s.cells().x_center.row(nghbr)};
        auto points = [&](SInd p) {
          return p == 0? cs.centroid : x_n;
        };
        auto values = [&](SInd p) {
          return p == 0? value_srfc(ghostIdx) : cell_values(nghbr);
        };
        auto x_gc = NumA<std::decay_t<S>::nd>{s.cells().x_center.row(ghostIdx)};
        return geometry::algorithm::linear_interpolation(x_gc, points, values);
      } else {
        /// all neighbors are moving boundary cells
        /// the right value will be interpolated at the ghost cell centroid
        /// automatically
      }
    }

    std::vector<NumA<std::decay_t<S>::nd>> xs;
    std::vector<Num> vs;
    xs.reserve(nghbrs.size() + 3); vs.reserve(nghbrs.size() + 3);

    xs.push_back(cs.centroid);      vs.push_back(value_srfc(ghostIdx));
    xs.push_back(cs.cut_points[0]); vs.push_back(value_srfc(ghostIdx));
    xs.push_back(cs.cut_points[1]); vs.push_back(value_srfc(ghostIdx));

    /// add neighbors to stencil
    for (auto n : nghbrs) {
      if(is_valid(s.cut_by_which_moving_boundary(n))) {
        // if the nghbr is a mb cell, add its surface centroid to stencil
        xs.push_back(s.grid().cut_surface_centroid(s.node_idx(n), g));
        vs.push_back(value_srfc(n));
      } else {
        xs.push_back(s.cells().x_center.row(n));
        vs.push_back(cell_values(n));
      }
    }

    NumA<std::decay_t<S>::nd> x_gc = s.cells().x_center.row(ghostIdx);
    auto v_gc = s.interpolate(x_gc, xs, vs);
    return v_gc;
  }

  template<class CV, class S, class G, class V>
  inline Num mb_neumann(CV&& cell_values, S&& s, G&& g, const CellIdx ghostIdx,
                        V&& flux_srfc) const noexcept {
    #define NAIVE
    //#define SIMPLE
    // #define ITERATIVE
    // #define ACCURATE

    #ifdef NAIVE
    std::vector<CellIdx> nghbrs; nghbrs.reserve(6);
    for (auto nghbr_pos : s.grid().neighbor_positions()) {
      const auto nghbrIdx = s.cells().neighbors(ghostIdx, nghbr_pos);
      if (is_valid(nghbrIdx) && !s.is_ghost_cell(nghbrIdx)
          && !is_valid(s.cut_by_which_moving_boundary(nghbrIdx))) {
        nghbrs.push_back(nghbrIdx);
      }
    }

    if (nghbrs.size() == 0) {
      auto average_v = 0.;
      auto count = 0;
      for (auto nghbr_pos : s.grid().neighbor_positions()) {
        const auto nghbrIdx = s.cells().neighbors(ghostIdx, nghbr_pos);
        if (is_valid(nghbrIdx)) {
          average_v += cell_values(nghbrIdx);
          count++;
        }
      }
      average_v /= count;
      return average_v;
    } else if (nghbrs.size() == 1) {
      auto bndryIdx = nghbrs[0];
      return this->neumann(cell_values(bndryIdx), flux_srfc(ghostIdx),
                           s.cell_dx(bndryIdx, ghostIdx));
    } else if (nghbrs.size() == 2) {
      auto bndryIdx0 = nghbrs[0];
      auto bndryIdx1 = nghbrs[1];
      auto dx0 = s.cell_dx(bndryIdx0, ghostIdx);
      auto dx1 = s.cell_dx(bndryIdx1, ghostIdx);
      auto T0 = cell_values(bndryIdx0);
      auto T1 = cell_values(bndryIdx1);
      auto T = 0.5 * (T0 + T1);
      auto dx = 0.5 * (dx0 + dx1);
      return this->neumann(T, flux_srfc(ghostIdx), dx);
    } else {
      TERMINATE("UNKNOWN NO NGHBRS!");
    }
    #endif
    ////////////////////////////////////////////////////////////////////////////
    #ifdef SIMPLE
    auto gradient = flux_srfc(ghostIdx);

    /// cut surface:
    auto cs = s.grid().cut_surface(s.node_idx(ghostIdx), g);

    /// In plane surface vector
    auto srfc = Eigen::ParametrizedLine<Num, std::decay_t<S>::nd>::Through
                (cs.cut_points[0], cs.cut_points[1]);

    /// Line perpendicular to surface crossing the cell center
    NumA<std::decay_t<S>::nd> x_gc = s.cells().x_center.row(ghostIdx);
    auto line = Eigen::ParametrizedLine<Num, std::decay_t<S>::nd>(x_gc, cs.normal);
    auto x_ip = line.intersectionPoint
                (Eigen::Hyperplane<Num, std::decay_t<S>::nd>(srfc));

    DBGV((ghostIdx)(x_gc)(s.cells().length(ghostIdx))
            (cs.centroid)(cs.cut_points[0])(cs.cut_points[1])(cs.normal)(x_ip));

    auto sgc = math::signum(g(x_gc));

    /// Interpolation neighbors
    auto ip_nghbrs = s.find_interpolation_neighbors(ghostIdx);

    ASSERT(ip_nghbrs.size() != 0, "no nghbrs found?!");
    std::vector<CellIdx> nghbrs; nghbrs.reserve(10);
    std::size_t mbNghbrs = 0;
    for(auto n : ip_nghbrs) {
      if(s.is_ghost_cell(n) && is_valid(s.cut_by_which_moving_boundary(n))) {
        ++mbNghbrs;
      } else {
        nghbrs.push_back(n);
      }
    }
    auto noNghbrs = nghbrs.size();
    auto noIpNghbrs = ip_nghbrs.size();
    if (!(ip_nghbrs.size() - mbNghbrs > 1)) {
      if (ip_nghbrs.size() == mbNghbrs) {
        /// all neighbors are moving boundary cells
        auto average_v = 0.;
        for(auto n : ip_nghbrs) {
          average_v += cell_values(n);
        }
        average_v /= ip_nghbrs.size();
        return average_v;
      } else {
        ASSERT(nghbrs.size() > 0, "cell: " << ghostIdx
               << " has no neighbors for mb neumann condition");

        /// only one neighbor is not a mb neighbor
        CellIdx nghbr = invalid<CellIdx>();
        for(auto n : nghbrs) {
          if(s.is_ghost_cell(n) && is_valid(s.cut_by_which_moving_boundary(n))) {
            continue;
          }
          nghbr = n;
          break;
        }
        ASSERT(is_valid(nghbrs), "WTF");
        auto bndryIdx = nghbr;
        return this->neumann(cell_values(bndryIdx), flux_srfc(ghostIdx),
                             s.cell_dx(bndryIdx, ghostIdx));
      }
    }

    ASSERT(nghbrs.size() > 0, "cell: " << ghostIdx
           << " has no neighbors for mb neumann condition");

    DBGV((ghostIdx)(x_gc)(sgc)(s.cells().length(ghostIdx))(cs.centroid)
         (cs.cut_points[0])(cs.cut_points[1])(cs.normal)(x_ip)
         (noNghbrs)(mbNghbrs)(noIpNghbrs));
    std::cerr << "nghbrIds: ";
    for (auto n : nghbrs) { std:: cerr << n << " "; }
    std::cerr << std::endl;

    std::vector<NumA<std::decay_t<S>::nd>> xs;
    std::vector<Num> vs;
    xs.reserve(nghbrs.size() + 3); vs.reserve(nghbrs.size() + 3);
    for (auto n : nghbrs) {
      xs.push_back(s.cells().x_center.row(n));
      vs.push_back(cell_values(n));
    }
    ASSERT(xs.size() > 0, "error: empty stencil!");

    // Interpolate value at intersection point with surface:
    auto v_ip = s.interpolate(x_ip, xs, vs);
    Num v_gc = 0;
    if(sgc > 0.) {
      v_gc = v_ip + gradient * (x_gc - x_ip).norm();
    } else {
      v_gc = v_ip - gradient * (x_gc - x_ip).norm();
    }

    return v_gc;
    #endif
    ////////////////////////////////////////////////////////////////////////////
    #ifdef ITERATIVE
    auto gradient = flux_srfc(ghostIdx);

    /// cut surface:
    auto cs = s.grid().cut_surface(s.node_idx(ghostIdx), g);

    /// In plane surface vector
    auto srfc = Eigen::ParametrizedLine<Num, std::decay_t<S>::nd>::Through
                (cs.cut_points[0], cs.cut_points[1]);

    /// Line perpendicular to surface crossing the cell center
    NumA<std::decay_t<S>::nd> x_gc = s.cells().x_center.row(ghostIdx);
    auto line = Eigen::ParametrizedLine<Num, std::decay_t<S>::nd>(x_gc, cs.normal);
    auto x_ip = line.intersectionPoint
                (Eigen::Hyperplane<Num, std::decay_t<S>::nd>(srfc));

    DBGV((ghostIdx)(x_gc)(s.cells().length(ghostIdx))
            (cs.centroid)(cs.cut_points[0])(cs.cut_points[1])(cs.normal)(x_ip));


    /// Interpolation neighbors
    auto ip_nghbrs = s.find_interpolation_neighbors(ghostIdx);

    ASSERT(ip_nghbrs.size() != 0, "no nghbrs found?!");
    std::vector<CellIdx> nghbrs; nghbrs.reserve(10);
    std::size_t mbNghbrs = 0;
    for(auto n : ip_nghbrs) {
      if(s.is_ghost_cell(n) && is_valid(s.cut_by_which_moving_boundary(n))) {
        ++mbNghbrs;
      } else {
        nghbrs.push_back(n);
      }
    }
    auto noNghbrs = nghbrs.size();
    auto noIpNghbrs = ip_nghbrs.size();
    if (!(ip_nghbrs.size() - mbNghbrs > 1)) {
      if (ip_nghbrs.size() == mbNghbrs) {
        /// all neighbors are moving boundary cells
        auto average_v = 0.;
        for(auto n : ip_nghbrs) {
          average_v += cell_values(n);
        }
        average_v /= ip_nghbrs.size();
        return average_v;
      } else {
        ASSERT(nghbrs.size() > 0, "cell: " << ghostIdx
               << " has no neighbors for mb neumann condition");

        /// only one neighbor is not a mb neighbor
        CellIdx nghbr = invalid<CellIdx>();
        for(auto n : nghbrs) {
          if(s.is_ghost_cell(n) && is_valid(s.cut_by_which_moving_boundary(n))) {
            continue;
          }
          nghbr = n;
          break;
        }
        ASSERT(is_valid(nghbrs), "WTF");
        auto bndryIdx = nghbr;
        return this->neumann(cell_values(bndryIdx), flux_srfc(ghostIdx),
                             s.cell_dx(bndryIdx, ghostIdx));
        // note: this cell won't be part of any stencil
        // ASSERT(nghbrs.size() > 0, "no neighbors found!!");
        // CellIdx nghbr = invalid<CellIdx>();
        // for(auto n : nghbrs) {
        //   if(s.is_ghost_cell(n) && is_valid(s.cut_by_which_moving_boundary(n))) {
        //     continue;
        //   }
        //   nghbr = n;
        //   break;
        // }
        // ASSERT(is_valid(nghbr), "WTF: ghostIdx: " << ghostIdx << " has: "
        //        << nghbrs.size() << " interpolation neighbors, of which "
        //        << mbNghbrs << " are mb ghost cells");
        // auto x_n = NumA<std::decay_t<S>::nd>{s.cells().x_center.row(nghbr)};

        // auto distance = (x_n - x_gc).norm();
        // auto v_n = cell_values(nhgbr);


        // auto points = [&](SInd p) {
        //   return p == 0? cs.centroid : x_n;
        // };
        // auto values = [&](SInd p) {
        //   return p == 0? value_srfc(ghostIdx) : cell_values(nghbr);
        // };
        // auto x_gc = NumA<std::decay_t<S>::nd>{s.cells().x_center.row(ghostIdx)};
        // return geometry::algorithm::linear_interpolation(x_gc, points, values);
      }
    }

    ASSERT(nghbrs.size() > 0, "cell: " << ghostIdx
           << " has no neighbors for mb neumann condition");

    auto sgc = math::signum(g(x_gc));
    DBGV((ghostIdx)(x_gc)(sgc)(s.cells().length(ghostIdx))
            (cs.centroid)(cs.cut_points[0])(cs.cut_points[1])(cs.normal)(x_ip)
            (noNghbrs)(mbNghbrs)(noIpNghbrs));
    std::cerr << "nghbrIds: ";
    for (auto n : nghbrs) { std:: cerr << n << " "; }
    std::cerr << std::endl;

    std::vector<NumA<std::decay_t<S>::nd>> xs;
    std::vector<Num> vs;
    xs.reserve(nghbrs.size() + 3); vs.reserve(nghbrs.size() + 3);
    for (auto n : nghbrs) {
      xs.push_back(s.cells().x_center.row(n));
      vs.push_back(cell_values(n));
    }
    ASSERT(xs.size() > 0, "error: empty stencil!");

    /// Compute gradient
    auto vi_ip = s.interpolate(x_ip, xs, vs);
    auto vi_gc = s.interpolate(x_gc, xs, vs);
    auto dx = (x_gc - x_ip).norm();
    Num gi = 0.;
    if(sgc > 0.) {
      gi = !math::approx(dx, 0.)? (vi_gc - vi_ip) / dx : 0.;
    } else {
      gi = !math::approx(dx, 0.)? -(vi_gc - vi_ip) /dx : 0.;
    }

    auto xs_ip = xs; xs_ip.push_back(x_gc);
    auto vs_ip = vs; vs_ip.push_back(vi_gc);

    Num error = std::abs(gradient - gi);
    const Num surface_gradient_accuracy = 1e-5;
    SInd no_iterations = 0;
    const SInd limit_iterations = 200;
    while (error > surface_gradient_accuracy) {
      if (sgc > 0.) {
        vi_gc = (gradient - gi) * dx + vi_ip;
      } else {
        vi_gc = -(gradient - gi) * dx + vi_ip;
      }
      vs_ip.back() = vi_gc;
      vi_ip = s.interpolate(x_ip, xs_ip, vs_ip);  // compute new surface value
      if(sgc > 0.) {
        gi = !math::approx(dx, 0.)? (vi_gc - vi_ip) / dx : 0.;
      } else {
        gi = !math::approx(dx, 0.)? -(vi_gc - vi_ip) /dx : 0.;
      }
      error = std::abs(gradient - gi);
      ++no_iterations;
      DBGV_ON((no_iterations)(vi_gc)(vi_ip)(gradient)(gi)(error)(dx));
      ASSERT(no_iterations < limit_iterations, "Limit iterations achieved!"
             << " Error is: " << std::abs(error));
    }

    return vi_gc;

    // // Interpolate value at intersection point with surface:
    // auto v_ip = s.interpolate(x_ip, xs, vs);

    // Num v_gc = 0;
    // if(sgc > 0.) {
    //   v_gc = v_ip + gradient * dx;
    // } else {
    //   v_gc = v_ip - gradient * dx;
    // }


    // auto xs_ip = xs; xs_ip.push_back(x_gc);
    // auto vs_ip = vs; vs_ip.push_back(v_gc);

    // // add ghost cell to the stencil of the ip
    // auto xs_ip = xs; xs_ip.push_back(x_gc);
    // auto vs_ip = vs; vs_ip.push_back(v_gc);

    // v_ip = s.interpolate(x_ip, xs_ip, vs_ip);

    // Num error = std::abs(dx * gradient - (v_gc - v_ip));

    // const Num surface_gradient_accuracy = 1e-5;
    // SInd no_iterations = 0;
    // const SInd limit_iterations = 200;
    // while (error > surface_gradient_accuracy) {
    //   if(sgc > 0.) {
    //     v_gc = v_ip + gradient * dx;
    //   } else {
    //     v_gc = v_ip - gradient * dx;
    //   }

    //   vs_ip.back() = v_gc;  // Add to ip stencil
    //   v_ip = s.interpolate(x_ip, xs_ip, vs_ip);  // compute new surface value
    //   error = std::abs(dx * gradient - (v_gc - v_ip));  // update error
    //   ++no_iterations;

    //   DBGV_ON((no_iterations)(v_gc)(v_ip)(gradient)(error)(dx));
    // //   ASSERT(no_iterations < limit_iterations, "Limit iterations achieved!"
    // //          << " Error is: " << std::abs(error));
    // }

    // return v_gc;
    #endif

    ////////////////////////////////////////////////////////////////////////////
    #ifdef ACCURATE
    auto gradient = flux_srfc(ghostIdx);

    /// cut surface:
    auto cs = s.grid().cut_surface(s.node_idx(ghostIdx), g);

    /// In plane surface vector
    auto srfc = Eigen::ParametrizedLine<Num, std::decay_t<S>::nd>::Through
                (cs.cut_points[0], cs.cut_points[1]);

    /// Line perpendicular to surface crossing the cell center
    NumA<std::decay_t<S>::nd> x_gc = s.cells().x_center.row(ghostIdx);
    auto line = Eigen::ParametrizedLine<Num, std::decay_t<S>::nd>(x_gc, cs.normal);
    auto x_ip = line.intersectionPoint
                (Eigen::Hyperplane<Num, std::decay_t<S>::nd>(srfc));

    /// Compute the mirror point
    auto sgc = math::signum(g(x_gc));
    auto dist = (x_ip - x_gc).norm();
    auto x_mp = sgc > 0.  // Is always inside the domain
                ? line.pointAt(dist)  // note: line starts at gc
                : line.pointAt(2. * dist);

    DBGV((ghostIdx)(x_gc)(s.cells().length(ghostIdx))
         (cs.centroid)(cs.cut_points[0])(cs.cut_points[1])(cs.normal)
         (x_ip)(x_mp));

    /// Interpolation neighbors
    auto ip_nghbrs = s.find_interpolation_neighbors(ghostIdx);

    ASSERT(ip_nghbrs.size() != 0, "no nghbrs found?!");
    std::vector<CellIdx> nghbrs; nghbrs.reserve(10);
    std::size_t mbNghbrs = 0;
    for(auto n : ip_nghbrs) {
      if(s.is_ghost_cell(n) && is_valid(s.cut_by_which_moving_boundary(n))) {
        ++mbNghbrs;
      } else {
        nghbrs.push_back(n);
      }
    }
    auto noNghbrs = nghbrs.size();
    auto noIpNghbrs = ip_nghbrs.size();
    if (!(ip_nghbrs.size() - mbNghbrs > 1)) {
      if (ip_nghbrs.size() == mbNghbrs) {
        /// all neighbors are moving boundary cells -> do nothing meaningful
        auto average_v = 0.;
        for(auto n : ip_nghbrs) {
          average_v += cell_values(n);
        }
        average_v /= ip_nghbrs.size();
        return average_v;
      }
    }

    ASSERT(nghbrs.size() > 0, "cell: " << ghostIdx
           << " has no neighbors for mb neumann condition");

    DBGV((ghostIdx)(x_gc)(sgc)(s.cells().length(ghostIdx))
            (cs.centroid)(cs.cut_points[0])(cs.cut_points[1])(cs.normal)(x_ip)
            (noNghbrs)(mbNghbrs)(noIpNghbrs));
    std::cerr << "nghbrIds: ";
    for (auto n : nghbrs) { std:: cerr << n << " "; }
    std::cerr << std::endl;

    std::vector<NumA<std::decay_t<S>::nd>> xs;
    std::vector<Num> vs;
    xs.reserve(nghbrs.size() + 3); vs.reserve(nghbrs.size() + 3);
    for (auto n : nghbrs) {
      xs.push_back(s.cells().x_center.row(n));
      vs.push_back(cell_values(n));
    }
    ASSERT(xs.size() > 0, "error: empty stencil!");

    /// interpolate values to mirror point
    auto v_mp = s.interpolate(x_mp, xs, vs);
    Num v_gc = 0.;
    /// set ghost cell value
    if(sgc > 0.) {
      /// impose gradient in the ghost cell center
      v_gc = v_mp - gradient * dist;
    } else {
      /// impose gradient in the intersection point
      v_gc = v_mp - gradient * 2. * dist;
    }

    return v_gc;
    #endif
  }
};

}  // namespace bc

////////////////////////////////////////////////////////////////////////////////
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif

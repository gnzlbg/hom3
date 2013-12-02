#ifndef HOM3_SOLVERS_FV_HEAT_MOVING_BOUNDARY_CONDITIONS_HPP_
#define HOM3_SOLVERS_FV_HEAT_MOVING_BOUNDARY_CONDITIONS_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief Implements the Heat Equantion's moving boundary conditions.
////////////////////////////////////////////////////////////////////////////////
#include "solver/fv/boundary_condition.hpp"
#include "interpolation/rbf.hpp"
#include "geometry/algorithms.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace solver { namespace fv { namespace heat {
////////////////////////////////////////////////////////////////////////////////
namespace bc {
////////////////////////////////////////////////////////////////////////////////

/// \brief Moving boundary conditions
namespace mb {

/// \brief Dirichlet boundary condition for the temperature
template<class Solver> struct Dirichlet : fv::bc::Condition<Dirichlet<Solver>> {
  /// \brief Imposes a \p temperature distribution
  template<class F, class G> Dirichlet(Solver& solver, G geometry,
                                       F&& temperature) noexcept
      : T_srfc(temperature), s(solver),
      g([=](const NumA<Solver::nd> x) { return (*geometry)(x); }) {}

  /// \brief Imposes a constant surface \p temperature
  template<class G>
      Dirichlet(Solver& solver, G geometry, const Num temperature) noexcept
      : T_srfc([=](CellIdx) { return temperature; }), s(solver),
      g([=](const NumA<Solver::nd> x) { return (*geometry)(x); }) {}

  /// \brief Applies the boundary condition to a range of ghost cells
  template<class _>
  void operator()(_, const CellIdx ghostIdx) const noexcept {
     /// add surface centroid/cut_points to stencil
    auto cs = s.grid().cut_surface(s.node_idx(ghostIdx), g);

    auto nghbrs = s.find_interpolation_neighbors2(ghostIdx);
    SInd mbNghbrs = 0;
    for(auto n : nghbrs) {
      if(s.is_ghost_cell(n)
         && is_valid(s.cut_by_which_moving_boundary(n))) {
        ++mbNghbrs;
      }
    }

    if(!(nghbrs.size() - mbNghbrs > 1)) {
      /// only one neighbor is not a mb neighbor
      // note: this cell won't be part of any stencil
      ASSERT(nghbrs.size() > 0, "no neighbors found!!");
      CellIdx nghbr;
      for(auto n : nghbrs) {
        if(s.is_ghost_cell(n) && is_valid(s.cut_by_which_moving_boundary(n))) {
          continue;
        }
        nghbr = n;
        break;
      }
      auto points = [&](SInd p) {
        return p == 0? cs.centroid
        : NumA<Solver::nd>(s.cells().x_center.row(nghbr));
      };
      auto values = [&](SInd p) {
        return p == 0? T_srfc(ghostIdx) : s.T(_(), nghbr);
      };
      s.T(_(), ghostIdx)
          = geometry::algorithm::linear_interpolation
          (NumA<Solver::nd>(s.cells().x_center.row(ghostIdx)),
           points, values);
      return;

    }

    std::vector<NumA<Solver::nd>> xs;
    std::vector<Num> vs;
    xs.reserve(nghbrs.size() + 3); vs.reserve(nghbrs.size() + 3);

    /// add surface centroid/cut_points to stencil
        //auto cs = s.grid().cut_surface(s.node_idx(ghostIdx), g);
    xs.push_back(cs.centroid);      vs.push_back(T_srfc(ghostIdx));
    xs.push_back(cs.cut_points[0]); vs.push_back(T_srfc(ghostIdx));
    xs.push_back(cs.cut_points[1]); vs.push_back(T_srfc(ghostIdx));

    /// add neighbors to stencil
    for (auto n : nghbrs) {
      if(is_valid(s.cut_by_which_moving_boundary(n))) {
        // if the nghbr is a mb cell, add its surface centroid to stencil
        xs.push_back(s.grid().cut_surface_centroid(s.node_idx(n), g));
        vs.push_back(T_srfc(n));
      } else {
        xs.push_back(s.cells().x_center.row(n));
        vs.push_back(s.T(_(), n));
      }
    }

    namespace rip = interpolation::rbf;
    // const auto kernel
    //     = rip::kernel::ThinPlate();
    // NumA<Solver::nd> x_gc = s.cells().x_center.row(ghostIdx);
    // auto v_gc = rip::interpolate(x_gc, xs, rip::build_weights(xs, vs, kernel), kernel);

    const auto kernel
        = rip::kernel::InverseMultiquadric(1.0);
    NumA<Solver::nd> x_gc = s.cells().x_center.row(ghostIdx);
    auto v_gc = rip::interpolate(x_gc, xs, rip::build_weights(xs, vs, kernel), kernel);

    // const auto minmax = std::minmax(std::begin(vs), std::end(vs));

    // if (v_gc < (*minmax.first)
    //     && v_gc < T_srfc(ghostIdx)) { v_gc = std::min((*minmax.first),
    //                                                   T_srfc(ghostIdx)); }
    // if (v_gc > (*minmax.second)
    //     && v_gc > T_srfc(ghostIdx)) { v_gc = std::max((*minmax.second),
    //                                                   T_srfc(ghostIdx)); }
    // if(!s.is_in_positive_lsv(ghostIdx)) { v_gc = T_srfc(ghostIdx); }

    s.T(_(), ghostIdx) = v_gc;
  }

  bool is_moving() const noexcept { return true; }

  /// Surface temperature distribution
  const std::function<Num(CellIdx)> T_srfc;
  Solver& s;
  const std::function<Num(NumA<Solver::nd>)> g;
};


}  // namespace mb

////////////////////////////////////////////////////////////////////////////////
}  // namespace bc
}  // namespace heat
}  // namespace fv
}  // namespace solver
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#endif

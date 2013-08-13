#ifndef HOM3_GRID_GENERATION_HPP_
#define HOM3_GRID_GENERATION_HPP_
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file implements mesh generation back-ends.
////////////////////////////////////////////////////////////////////////////////
/// Options:
#define ENABLE_DBG_ 0
#include "../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace grid {
////////////////////////////////////////////////////////////////////////////////

/// \brief Contains mesh generators
namespace generation {

/// \brief Mesh generator interface
///
/// A mesh generator must implement:
///   - a generic operator()(Grid&)
///
/// It is good practice to pass mesh generator properties via an io::Properties
/// in its constructor but this is not required.
///
template<class Implementation> struct Interface {
  template<class Grid>
  void operator()(Grid& g) { imp()->generate_mesh(g); }
 private:
        Implementation* imp()       { return static_cast<Implementation*>(this); }
  const Implementation* imp() const { return static_cast<Implementation*>(this); }
};

/// \brief Min-Level mesh generator: refines the grid up to a specified minimum
/// refinement \p level (see constructor)
///
struct MinLevel : Interface<MinLevel> {
  const Ind minDesiredLvl;
  MinLevel(io::Properties input) : minDesiredLvl(io::read<Ind>(input,"level"))
  { TRACE_IN_(); TRACE_OUT(); }

  template<class Grid> void generate_mesh(Grid& g) {
    TRACE_IN_();
    bool done = false;
    std::cerr << "minDesLvl: " << minDesiredLvl << std::endl;
    while(!done) {
      done = true;
      for(auto&& nIdx : g.nodes().nodes()) {
        if(g.nodes().is_leaf(nIdx) && g.nodes().level(nIdx) != minDesiredLvl) {
          g.refine_node(nIdx);
          done = false;
        }
      }
    }
    TRACE_OUT();
  }
};

////////////////////////////////////////////////////////////////////////////////
}}} // hom3::grid::generation namespace
////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif

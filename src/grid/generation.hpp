#ifndef HOM3_GRID_GENERATION_HPP_
#define HOM3_GRID_GENERATION_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Options:
#define ENABLE_DBG_ 0
#include "../misc/dbg.hpp"
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
namespace grid {
////////////////////////////////////////////////////////////////////////////////
namespace generation {
////////////////////////////////////////////////////////////////////////////////

template<class Implementation>
struct Interface {
  template<class Grid>
  void operator()(Grid* g) { imp()->generate_mesh(g); }
 private:
        Implementation* imp()       { return static_cast<Implementation*>(this); }
  const Implementation* imp() const { return static_cast<Implementation*>(this); }
};

struct MinLevel : Interface<MinLevel> {
  const Ind minDesiredLvl;
  MinLevel(io::Properties input) : minDesiredLvl(io::read<Ind>(input,"level"))
  { TRACE_IN_(); TRACE_OUT(); }

  template<class Grid> void generate_mesh(Grid* g) {
    TRACE_IN_();
    bool done = false;
    std::cerr << "minDesLvl: " << minDesiredLvl << std::endl;
    while(!done) {
      done = true;
      for(auto&& nId : g->nodes().nodes()) {
        if(g->nodes().is_leaf(nId) && g->nodes().level(nId) != minDesiredLvl) {
          g->refine_cell(nId);
          done = false;
        }
      }
    }
    TRACE_OUT();
  }
};

////////////////////////////////////////////////////////////////////////////////
} // namespace generation
////////////////////////////////////////////////////////////////////////////////
} // namespace grid
////////////////////////////////////////////////////////////////////////////////


/// Grid generation:
// namespace default_meshgen_ {

// template<class C> struct grid_traits {
//   typedef typename C::UNKNOWN_GRID_CONNECTIVITY_TYPE grid_tag;
// };
// struct structured_tag{};
// template<> struct grid_traits<structured_connectivity>{typedef structured_tag type;};
// }


// struct DefaultMeshGeneration {
//   template<class GridT> static void createMesh_(GridT* grid) {
//     typedef typename default_meshgen_::grid_traits<typename GridT::ConnectivityT::type>::type grid_tag;
//     createMeshImpl_(grid, grid_tag());
//   }


//   template<class GridT> static void createMeshImpl_(GridT* grid, default_meshgen_::structured_tag ){
//     std::cout << "CREATE STRUCTURED MESH" << std::endl;
//     typedef typename Traits<typename GridT::ElementT,Structured>::node_type ConnectivityNode;
//     typedef NodeContainer<ConnectivityNode> ConnectivityNodeContainer;

//     IndT noGhostElements = 2;
//     IndT noElements = io::cast<IndT>(grid->properties_["noElements"]) + noGhostElements;
//     IndT noInternalElements = noElements - noGhostElements;
//     NumT dx = ( grid->BBox_.xMax[0] - grid->BBox_.xMin[0] ) / noElements;

//     //! Get grid connectivity, element container, and accessor:
//     ConnectivityNodeContainer* connectivityNodes = grid->connectivity_.getNodeContainer();
//     ElementContainer<typename GridT::ElementT>* elements;


//     //! Create grid elements:
//     typename GridT::ElementT element;
//     element.order = io::cast<IntT>(grid->properties_["elementOrder"]);
//     elements->push_back(noElements,element);


//     //! Set connectivity and coordinates:
//     NumT xi = grid->BBox_.xMin[0] - 0.5*dx;
//     std::for_each(elements->begin(),elements->end(),
//                   [&](typename GridT::ElementT& e){
//                     int nNodes = e.nNodes();
//                     ASSERT(nNodes > 1,"ERROR: nNodes <= 1");
//                     NumT dxLocal = dx/(nNodes-1);
//                     for(int n = 0; n < nNodes; ++n)
//                       e.data[ grid->x(n,xi+n*dxLocal) ];
//                     xi += dx;
//                     connectivityNodes->push_back(ConnectivityNode(&e));
//                   });
//   }
// };

////////////////////////////////////////////////////////////////////////////////
#undef ENABLE_DBG_
////////////////////////////////////////////////////////////////////////////////
#endif

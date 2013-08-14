#ifndef HOM3_OUTPUT_
#define HOM3_OUTPUT_
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include "../globals.hpp"
#include "../misc/helpers.hpp"
////////////////////////////////////////////////////////////////////////////////
namespace hom3 {

namespace grid { template<SInd nd> struct CartesianHSP; } // grid namespace

namespace io {
////////////////////////////////////////////////////////////////////////////////

/// Some notes: the io classes are able to write "streamable" grids with a
/// variable amount of data.
/// Any object that satisfies the "streamable" grid interface can be written out.
/// First the output is initialized as follows:

/// Here fName is the file name, and grid_information is an streamable grid
/// interface. The Output format ascii/binary can be customized as well as
/// the precision of the floating point data.
/// io::Vtk out(fName, io::grid_information(this), format = io::format::ascii,
///                                        precision = io::precision::default);
///
/// To write data fields you just stream them to the output:
/// out << io::stream(CV->V) << io::stream(m_coordinates) << io::stream(cells());
/// The output class will count for you how many variables you want to stream,
/// and stream everything to disk when the destructor is called.


/// VTK IO Class, this should be abstracted to a general IO class.

////////////////////////////////////////////////////////////////////////////////
/// IO Tags: output format, precision, output field types...
namespace format     { struct binary   {}; struct ascii  {}; }
namespace precision  { struct standard {}; struct low    {}; }
namespace data_types { struct scalar   {}; struct vector {}; }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Streamable Variable Interface:
///
/// Requirements on T:
///   value_type: underlying value type
///   name: stream name
///   dimensions: stream dimension
///   operator(Ind i, SInd d): value at element i for dimension d
///
/// Provides: name, dimension, ()(i,d), and value_type() which returns
/// a Streamable::value_type enum
struct StreamableVariable {
  /// \brief Copy constructor
  StreamableVariable(const StreamableVariable& o)
      : value_type_(o.value_type_), as_int(o.as_int), as_num(o.as_num),
        name_(o.name_), noDims_(o.noDims_) {}
  /// \brief Constructor for Ind variables
  template<class T> StreamableVariable
  (const T& t, Ind, std::function<std::string(const SInd)> n, SInd nd)
      : value_type_(StreamableVariable::value_types::Ind),
        as_int([=](const Ind i, const SInd d){ return t(i,d); }),
        as_num([](const Ind, const SInd) { return invalid<Num>(); }),
        name_(n), noDims_(nd) {}
  /// \brief Constructor for SInd variables
  template<class T> StreamableVariable
  (const T& t, SInd, std::function<std::string(const SInd)> n, SInd nd)
      : value_type_(StreamableVariable::value_types::Ind),
        as_int([=](const Ind i, const SInd d){ return t(i,d); }),
        as_num([](const Ind, const SInd) { return invalid<Num>(); }),
        name_(n), noDims_(nd) {}
  /// \brief Constructor for Num variables
  template<class T> StreamableVariable
  (const T& t, Num, std::function<std::string(const SInd)> n, SInd nd)
      : value_type_(StreamableVariable::value_types::Num),
        as_int([](const Ind, const SInd){ return invalid<Ind>(); }),
        as_num([=](const Ind i, const SInd d){ return t(i,d); }),
        name_(n), noDims_(nd) {}
  /// \brief Constructor for types that satisfy the Streamable concept
  template<class T> StreamableVariable(const T& t)
      : StreamableVariable(t, typename T::value_type(), t.name(), t.dimensions()) {}
  /// \brief Constructor from function object
  template<class F> StreamableVariable(std::string n, SInd nd, F&& f)
      : StreamableVariable(f, decltype(f(0,0))(), make_name(n,nd), nd) {}
  template<class F> StreamableVariable
  (std::function<std::string(const SInd)> n, SInd nd, F&& f)
      : StreamableVariable(f, decltype(f(0,0))(), n, nd) {}

  std::string name(SInd d = 0) const { return name_(d); }
  SInd no_dimensions() const { return noDims_; }
  SInd value_type() const { return value_type_; }
  struct value_types { enum { Ind = 0, Num = 1 }; };
 private:
  const SInd value_type_;
 public:
  const std::function<Ind(const Ind, const SInd)> as_int;
  const std::function<Num(const Ind, const SInd)> as_num;
 private:
  //const std::string name_;
  const std::function<std::string(const SInd)> name_;
  const SInd noDims_;
  static std::function<std::string(const SInd)> make_name(std::string n, SInd nd) {
    if(nd > 1) {
      return [=](const SInd d) { return n + "_" + std::to_string(d); };
    } else {
      return [=](const SInd) { return n; };
    }
  }
};

/// "make_stream" helper function:
template<class Variable> StreamableVariable stream
(const Variable& v) { return StreamableVariable(v); }

template<class Functor> StreamableVariable stream
(std::string name, SInd nd, Functor&& f) {
  return StreamableVariable(name, nd, f);
}

template<class Functor> StreamableVariable stream
(std::function<std::string(const SInd)> name, SInd nd, Functor&& f) {
  return StreamableVariable(name, nd, f);
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// "Streamable" grid interface:
template<SInd nd> struct StreamableDomain {

  /// \todo this function should be removed eventually, output should be
  /// exclusively range-based
  StreamableDomain(const grid::CartesianHSP<nd>* g)
      : cells([=](){ return g->nodes().leaf_nodes()
              | boost::adaptors::transformed(
                  [&](const NodeIdx i){
                  return i();
                })
              ; }),
        cell_vertices([=](const Ind nIdx){
            return g->compute_cell_vertices(NodeIdx{nIdx}); })
  {}

  template<class CellRange, class CellVerticesRange>
  StreamableDomain(CellRange&& cells_, CellVerticesRange&& cell_vertices_)
      : cells(std::forward<CellRange>(cells_)),
        cell_vertices(std::forward<CellVerticesRange>(cell_vertices_)) {}

  const std::function<AnyRange<Ind>(void)> cells;
  const std::function<typename grid::CartesianHSP<nd>::CellVertices(Ind)> cell_vertices;

  inline auto dimensions() const -> decltype(Range<SInd>(SInd{0},nd)) {
    return Range<SInd>(SInd{0},nd);
  }
  static constexpr SInd no_dimensions() { return nd; }
};
////////////////////////////////////////////////////////////////////////////////

template<template <SInd> class Grid, SInd nd>
StreamableDomain<nd> stream_grid(const Grid<nd>& g) {
  return StreamableDomain<nd>(g);
}


////////////////////////////////////////////////////////////////////////////////
/// VTK Output class:
template<SInd nd, class Format = io::format::ascii> struct Vtk {

  /// The constructor configures the output and prepares the header information
  template<class Precision = io::precision::standard>
  Vtk(StreamableDomain<nd> d, const std::string fN, const Precision p = Precision())
      :  streams_(), fileName(fN + ".vtk"), domain(d), isFileOpen(false),
         isHeaderWritten(false) {
    set_precision(p);

    /// VTK Header Information:
    const std::string version_  = "# vtk DataFile Version 2.0\n";
    const std::string title_    = "HOM3 legacy vtk Output \n";
    const std::string format_   = header_format_string();
    const std::string meshType_ = "DATASET UNSTRUCTURED_GRID\n";

    header_ << version_ << title_ << format_ << meshType_;
  }

  ~Vtk() {

    /// Output log:
    std::cout << "#of grid dimensions: " << domain.no_dimensions() << "\n";
    for(auto stream : streams_) {
      std::cerr << "name: " << stream.name() << " dim: " << stream.no_dimensions() << "\n";
    }

    initializeOutput();

    write_to_file();

    finalizeOutput();
  }

  inline void write_stream(const StreamableVariable& stream) {
    for(SInd d = 0; d < stream.no_dimensions(); ++d) {
      write_scalar_stream(stream,stream.name(d), d);
    }
  }

  inline void write_scalar_stream(const StreamableVariable& stream, const std::string& name, SInd dim) {
    std::string type;
    switch (stream.value_type()) {
      case StreamableVariable::value_types::Ind: {
        type = "int";
        break;
      }
      case StreamableVariable::value_types::Num: {
        type = "float";
        break;
      }
      default: { TERMINATE("unknown stream value_type"); }
    }
    write_header("SCALARS " + name + " " + type + "\n");
    write_header("LOOKUP_TABLE default\n");
    for(auto cell : domain.cells()) {
      switch (stream.value_type()) {
        case StreamableVariable::value_types::Ind: {
          auto tmp = stream.as_int(cell,dim);
          write_formated_scalar(tmp);// << "\n";
          break;
        }
        case StreamableVariable::value_types::Num: {
          write_formated_scalar(static_cast<float>(stream.as_num(cell,dim)));// << "\n";
          break;
        }
        default: { TERMINATE("unknown stream value_type (with name: " + name + ")"); }
      }
    }
    write_new_line();
  }

  void write_to_file() {

    /// Cell types: 2d-pixel, 3d-voxel
    static constexpr SInd cellType = (nd == 2) ? 8 : 11;

    /// complexity: O(N)
    const Ind noCells = boost::distance(domain.cells());
    static constexpr Ind noCellVertices = math::ct::ipow(2u,nd);

    auto cell_positions = [&](){ return boost::counting_range(Ind(0),noCells); };
    auto vertex_positions = [&](){ return boost::counting_range(Ind(0),noCellVertices); };

    /// Store each vertex of each cell into a vector
    /// The vertices of each cell are stored in counter-clockwise order
    /// memory: O(N*M), where N = #cells, M=#vertices/cell
    /// complexity: O(N) = counting (O(N)) + copying (O(N))
    std::vector<Vertex> cell_vertices;

    cell_vertices.reserve(noCells * noCellVertices); // m: O(N*M)

    for(auto cellId : domain.cells()) {
      boost::copy(domain.cell_vertices(cellId)(), std::back_inserter(cell_vertices));
    }


    /// Create a unique cell vertices vector
    /// memory: O(N*M) memory | cumulative: 2*O(N+M)
    /// complexity: O((N*M)log(N*M)) | cumulative: O((N*M)log(N*M))
    auto unique_cell_vertices = cell_vertices;

    /// LessThanComparable: is a < b ?
    auto vertex_cmp = [](const Vertex& a, const Vertex& b) {
      if(a(0) > b(0)) { // \todo optimize?
        return false;
      } else if(math::apprx(a(0),b(0))) {
        if(a(1) > b(1)) {
          return false;
        } else if(math::apprx(a(1),b(1))) {
          if(nd == 2) {
            return false;
          } else if(nd == 3 && (a(2) > b(2) || math::apprx(a(2),b(2)))) {
            return false;
          }
        }
      }
      return true;
    };

    /// EqualityComparable: is a == b ?
    auto vertex_eq = [](const Vertex& a, const Vertex& b) { return a.isApprox(b, math::eps); };

    /// Removes the duplicated elements: erase(unique(sort(rng)))
    /// Note: resulting unique range remains sorted according to vertex_cmp
    boost::erase(unique_cell_vertices,
                 boost::unique<boost::return_found_end>(
                     boost::sort(unique_cell_vertices, vertex_cmp), vertex_eq));

    /// Write unique vertices (points) to file
    write_header("POINTS " + std::to_string(unique_cell_vertices.size())+ " FLOAT\n");
    for(auto&& vertex : unique_cell_vertices ) {
      for(auto d : domain.dimensions()) {
        write_formated_scalar(vertex(d));
      }
      if(nd == 2) { write_formated_scalar(0.0); }
    }
    write_new_line();

    /// Write cell-node incidence table The vertices in cell_vertices are
    /// already stored in the required counter-clockwise order (see vtk elements
    /// 8 and 11)
    write_header("CELLS " + std::to_string(noCells)
                 + " " + std::to_string(noCells*(noCellVertices + 1)) + "\n");

    for(const auto& cellPos : cell_positions()) {
      write_formated_scalar(noCellVertices);
      for(const auto& vertexPos : vertex_positions()) {
        auto vertexId = cellPos * noCellVertices + vertexPos;
        auto vertex = cell_vertices[vertexId];

        // (O(logN)):
        auto it = boost::lower_bound(unique_cell_vertices, vertex, vertex_cmp);
        const Ind uniqueVertexId = it - std::begin(unique_cell_vertices);
        write_formated_scalar(uniqueVertexId);

        ASSERT(it != std::end(unique_cell_vertices), "vertex not found!");
      }
    }
    write_new_line();

    /// Write cell-types:
    write_header("CELL_TYPES " +std::to_string(noCells) + "\n");
    boost::for_each(domain.cells(), [&](const Ind&) { write_formated_scalar(cellType); });
    write_new_line();

    /// Write cell data:
    write_header("CELL_DATA " + std::to_string(noCells) + "\n");
    for(const auto& stream : streams_) {
      write_stream(stream);
    }
  }

  inline void finalizeOutput() { os.close(); }

  Vtk& operator<<(const StreamableVariable& s) {
    streams_.push_back(s);
    return *this;
  }

  static const std::string write_dim_type(io::data_types::scalar) { return "SCALARS "; }
  static const std::string write_dim_type(io::data_types::vector) { return  "VECTOR "; }
  static const std::string write_val_type(int)   { return " int"; }
  static const std::string write_val_type(float) { return " float"; }

 private:
  void write_header(const std::string& i) { write_header(i,Format()); }
  void write_header(const std::string& i, io::format::binary) {
    if(os.is_open()) { os.close(); }
    os.open(fileName, std::ios_base::app); // reopen in ascii mode
    os << i;
    os.close();
    os.open(fileName, std::ios_base::app|std::ios_base::binary);
  }
  void write_header(const std::string& i, io::format::ascii) { os << i; }

  template<class T> int interpret_id(T&& t) { return is_valid<Ind>(t) ? static_cast<int>(t) : -1; }
  template<class T> float interpret_num(T&& t) { return static_cast<float>(t); }
  template<class T> inline void write_formated_scalar(T&& t) {
    using type = typename std::remove_reference<T>::type;
    static_assert( std::is_integral<type>::value || std::is_floating_point<type>::value,
                   "t is either an integral or a floating point value!");
    if(std::is_integral<type>::value) {
      write_formated_id(std::forward<T>(t),Format());
    } else if (std::is_floating_point<type>::value) {
      write_formated_num(std::forward<T>(t),Format());
    }
  }
  template<class T> inline void write_formated_num(T&& t, io::format::binary) {
    auto tmp = floatSwap(interpret_num(t));
    os.write(reinterpret_cast<const char*>(&tmp),sizeof(tmp));
  }
  template<class T> inline void write_formated_num(T&& t, io::format::ascii) {
    os << interpret_num(t) << " ";
  }
  template<class T> inline void write_formated_id(T&& t, io::format::binary) {
    auto tmp = intSwap(interpret_id(t));
    os.write(reinterpret_cast<const char*>(&tmp),sizeof(tmp));
  }
  template<class T> inline void write_formated_id(T&& t, io::format::ascii) {
    os << interpret_id(t) << " ";
  }

  inline void write_new_line() { write_new_line(Format()); }
  inline void write_new_line(io::format::ascii)  { os << "\n"; }
  inline void write_new_line(io::format::binary) { os << std::endl; }

  using Vertex = typename grid::CartesianHSP<nd>::CellVertices::Vertex;
  std::vector<StreamableVariable> streams_;
  const std::string fileName;
  std::stringstream header_;
  StreamableDomain<nd> domain;
  Ind precision;
  std::ofstream os;
  static const std::string header_format_string() { return header_format_string(Format()); }
  static const std::string header_format_string(io::format::ascii) { return "ASCII\n"; }
  static const std::string header_format_string(io::format::binary) { return "BINARY\n"; }
  void set_precision(io::precision::standard) { precision = 12; os.precision(precision); }
  void set_precision(io::precision::low) { precision = 6; os.precision(precision); }

  void initializeOutput() {
    os.open(fileName);
    if(!os) { TERMINATE("can't write to file"); }
    os << header_.str(); // writes header
  }

  bool isFileOpen, isHeaderWritten, isHeaderSet;

   inline int intSwap( int f ){ // Change endian of int
    union
    {
      int f;
      unsigned char b[4];
    } dat1, dat2;

    dat1.f = f;
    dat2.b[0] = dat1.b[3];
    dat2.b[1] = dat1.b[2];
    dat2.b[2] = dat1.b[1];
    dat2.b[3] = dat1.b[0];
    return dat2.f;

  };
  inline float floatSwap( float f ){ // Change endian of float
    union
    {
      float f;
      unsigned char b[4];
    } dat1, dat2;

    dat1.f = f;
    dat2.b[0] = dat1.b[3];
    dat2.b[1] = dat1.b[2];
    dat2.b[2] = dat1.b[1];
    dat2.b[3] = dat1.b[0];
    return dat2.f;

  };
};

////////////////////////////////////////////////////////////////////////////////
} } // hom3::io namespace
////////////////////////////////////////////////////////////////////////////////
#endif

#ifndef HOM3_MPI_HPP_
#define HOM3_MPI_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Includes:
#include <boost/mpi.hpp>
#include <type_traits>
////////////////////////////////////////////////////////////////////////////////
namespace boost { namespace mpi {
////////////////////////////////////////////////////////////////////////////////

namespace detail {

struct not_void{}; struct is_void{};
template<class T> struct is_void_type { using type = not_void; };
template<> struct is_void_type<void> { using type = is_void; };

template<class F, class... Args>
void root_do_(is_void, const communicator& comm, F&& f, Args&&... args) {
  if(comm.rank() == 0) {
    f(std::forward<Args>(args)...);
  }
  return;
}

template<class F, class... Args>
auto root_do_(not_void, const communicator& comm, F&& f, Args&&... args)
-> decltype(f(std::forward<Args>(args)...)) {
  decltype(f(std::forward<Args>(args)...)) result;
  if(comm.rank() == 0) {
    result = f(std::forward<Args>(args)...);
  }
  boost::mpi::broadcast(comm, result, 0);
  return result;
}

} // namespace detail

/// \brief Performs a computation _only_ at the root of the communicator
///
/// If the computation has a return type, this type is broadcasted to all
/// MPI Processes.
template<class F, class... Args>
auto root_do(const communicator& comm, F&& f, Args&&... args)
-> decltype(f(std::forward<Args>(args)...)) {
  using is_void = typename detail::is_void_type<
      decltype(f(std::forward<Args>(args)...))
    >::type;
  return detail::root_do_(is_void(),comm,
                          std::forward<F>(f),std::forward<Args>(args)...);
}

////////////////////////////////////////////////////////////////////////////////
}} // namespace boost::mpi
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
namespace hom3 {
////////////////////////////////////////////////////////////////////////////////

MPI_Info mpi_info() { return MPI_INFO_NULL; }

////////////////////////////////////////////////////////////////////////////////
} // hom3 namespace
////////////////////////////////////////////////////////////////////////////////
#endif

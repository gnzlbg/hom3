#ifndef HOM3_IO_HDF5_FILE_HPP_
#define HOM3_IO_HDF5_FILE_HPP_
////////////////////////////////////////////////////////////////////////////////
/// Include:
#include <hdf5.h>
#include <hdf5_hl.h>
#include "globals.hpp"
////////////////////////////////////////////////////////////////////////////////
/// File local macros (to be undefined at the end of the file):
#define assert_closed() ASSERT(!is_open(), "File is already open!")
#define assert_open()   ASSERT(is_open(),  "File is closed!")
#define herror(status)  handle_error((status), AT_)
#define assert_read_mode()                                              \
  static_assert(std::is_same<AccessMode, access_mode::read>::value,     \
                "This function can only be used in read-only mode!")
#define assert_write_mode()                                             \
  static_assert(std::is_same<AccessMode, access_mode::write>::value,    \
                "This function can only be used in write mode!")

////////////////////////////////////////////////////////////////////////////////
namespace hom3 { namespace io { namespace file {
////////////////////////////////////////////////////////////////////////////////

/// \brief File access modes
namespace access_mode {
struct read {};   ///< Read-only access tag
struct write {};  ///< Write access tag
}  // namespace access_mode

/// \brief HDF5 functionality
namespace hdf5 {

inline String extension() noexcept { return ".hdf5"; }

enum class array_order { row_major_order, col_major_order };

}  // namespace hdf5

/// \brief HDF5 file handler
template<class AccessMode> struct HDF5File {
  HDF5File(const String fileName,
           const boost::mpi::communicator& comm = boost::mpi::communicator())
  noexcept
    : fileName_(fileName)
    , comm_(comm)
    , isOpen_(false)
  { open(AccessMode()); }

  ~HDF5File() { close(); };  ///< \brief Releases file handler
  HDF5File() = delete;

  /// \brief Reads dataset \p n directly to the memory at \p target
  template<class T> inline void array
  (const String n, T *const target) { array_impl_(n, target); }

  /// \brief Reads dataset \p n into functor f
  template<class T> inline auto array
  (const String n, std::function<T(Ind, Ind)> f) { array_impl_(n, f); }

  /// \brief Writes dataset \p n from memory at \p data using
  /// the order \p no_rows, \p no_cols, \p order.
  template<class T> inline void array
  (const String n, const T *const data, const Ind no_rows,
  const Ind no_cols = 1,
  hdf5::array_order order = hdf5::array_order::col_major_order) noexcept
  { array_impl_(n, data, no_rows, no_cols, order); }

  /// \brief Writes \p data to file
  template<class R, class F> inline auto array
  (const String n, R&& range, F&& f, const Ind no_cols = 1,
  hdf5::array_order order = hdf5::array_order::col_major_order) noexcept
  -> decltype(typename std::remove_reference_t<R>::value_type{}, void())
  { array_impl_(n, range, std::forward<F>(f), no_cols, order); }

  /// \brief Read/Writes dataset's \p name attribute \p attrName to \p value
  template<class T> inline void attribute
  (const String name, const String attrName, T&& value) noexcept
  { attribute(AccessMode(), name, attrName, std::forward<T>(value)); }

 private:
  /// \todo Clean the implementation, API should remain stable
  const String fileName_;  ///< File name
  const boost::mpi::communicator comm_;  ///< Communicator of IO processes

  template<class T> inline void attribute
  (access_mode::read, const String name, const String attrName,
  T&& value) noexcept { get_attribute(name, attrName, std::forward<T>(value)); }

  template<class T> inline void attribute
  (access_mode::write, const String name, const String attrName,
  const T value) noexcept { set_attribute(name, attrName, value); }

  /// \brief Writes dataset's \p name attribute \p attrName to \p value
  template<class T> inline void set_attribute
  (const String name, const String attrName, const T value) noexcept
  { set_attribute_(name, attrName, value); }

  /// \brief Reads dataset's \p name attribute \p attrName to \p value
  template<class T> inline T get_attribute
  (const String name, const String attrName, T& value) noexcept {
    get_attribute_(name, attrName, value);
    return value;
  }

  template<class R, class F> inline auto array_impl_
  (const String n, R&& range, F&& f, const Ind no_cols = 1,
  hdf5::array_order order = hdf5::array_order::col_major_order) noexcept
  -> decltype(typename std::remove_reference_t<R>::value_type{}, void()) {
    array_impl_(AccessMode(), n, std::forward<R>(range),
                std::forward<F>(f), no_cols, order);
    set_attribute(n, "no_cols", no_cols);
  }


  template<class T>
  inline void array_impl_(const String n, T *const target) {
    assert_read_mode();
    herror(H5LTread_dataset(fileId_, n.c_str(), hdf5_t(T{}), target));
  }

  template<class T>
  inline auto array_impl_(const String n, std::function<T(Ind, Ind)> f) {
    assert_read_mode();
    using return_type = std::remove_reference_t<T>;
    Ind no_rows, no_cols;
    std::tie(no_rows, no_cols) = get_data_dimensions(n);
    auto b = buffer(no_rows * no_cols, return_type{});
    array_impl_(n, b.data());
    ordered_execute
        (Range<Ind>{Ind{0}, no_rows}, no_rows, no_cols, get_data_order(n),
         [&](const Ind offset, const Ind rowIdx, const Ind colIdx) {
          f(rowIdx, colIdx) = b(offset);
        });
  }

  /// \brief Writes \p data to file
  template<class T> inline void array_impl_
  (const String n, const T *const data,
  const Ind no_rows, const Ind no_cols = 1,
  hdf5::array_order order = hdf5::array_order::col_major_order) noexcept {
    array_impl_(AccessMode(), n, data, no_rows, no_cols);
    set_attribute(n, "no_rows", no_rows);
    set_attribute(n, "no_cols", no_cols);
    set_attribute(n, "order", array_order_to_value(order));
  }

  /// \brief Executes ternary F \p f in the specified \p order.
  ///
  /// R \p range is a range of row indices
  /// \p no_rows is the #of rows
  /// \p no_cols is the #of columns
  /// \p order is the data layout: row-major or col-major
  /// \p f is a ternary ficate taking a:
  ///   - one-dimensional index (a memoryOffset)
  ///   - and its corresponding row and column indices
  template<class R, class F> void ordered_execute
  (R range, const Ind no_rows, const Ind no_cols,
  const hdf5::array_order order, F f) const noexcept {
    if (order == hdf5::array_order::row_major_order) {
      for (auto&& i : range) {
        for (Ind d = 0; d < no_cols; ++d) {
          f(i * no_cols + d, i, d);
        }
      }
    } else if (order == hdf5::array_order::col_major_order) {
      for (Ind d = 0; d < no_cols; ++d) {
        for (auto&& i : range) {
          f(d * no_rows + i, i, d);
        }
      }
    } else {
      TERMINATE("Unknown order.");
    }
  }

  /// \name Read/Write array implementations
  ///@{

  /// \brief Writes 2D data
  template<class T> inline void array_impl_
  (access_mode::write, const String n, const T* const data,
  const Ind no_rows, const Ind no_cols = 1) noexcept {
    const std::array<hsize_t, 1> dim
      = {{ static_cast<hsize_t>(no_rows * no_cols) }};
    herror(H5LTmake_dataset(fileId_, n.c_str(), 1, dim.data(),
                            hdf5_t(T{}), data));
  }

  /// \brief Writes 2D data
  template<class R, class F> inline void array_impl_
  (access_mode::write, const String n, R&& range, F&& f,
  const Ind no_cols = 1,
  hdf5::array_order order = hdf5::array_order::col_major_order) noexcept {
    using T = typename std::remove_reference_t<R>::value_type;
    const Ind no_rows = boost::distance(range);
    auto b = buffer(no_rows * no_cols, T{});
    ordered_execute(range, no_rows, no_cols, order,
                    [&](const Ind lhsIdx, const Ind rowIdx, const Ind colIdx) {
                      b(lhsIdx) = f(rowIdx, colIdx);
                    });
    array(n, b.data(), no_rows, no_cols, order);
  }

  ///@}

  /// \name Type-map, hdf5_t: HOM3 Type -> HDF5 Type
  ///
  /// Note: The H5T_... are magic values of some type.
  ///@{

  auto hdf5_t(Ind)  RETURNS(H5T_NATIVE_LONG);
  auto hdf5_t(SInd) RETURNS(H5T_NATIVE_INT);
  auto hdf5_t(Num)  RETURNS(H5T_NATIVE_DOUBLE);

  ///@}

  /// \name Set attribute maps
  ///@{

  template<class T, class U> void assert_equal()
  { static_assert(std::is_same<T, U>::value, "type changed!"); }

  inline void set_attribute_
  (const String name, const String attrName, const Ind value) noexcept {
    assert_equal<Ind, long long>();
    herror(H5LTset_attribute_long_long
           (fileId_, name.c_str(), attrName.c_str(), &value, 1));
  }

  inline void set_attribute_
  (const String name, const String attrName, const SInd value) noexcept {
    assert_equal<SInd, int>();
    herror(H5LTset_attribute_int
           (fileId_, name.c_str(), attrName.c_str(), &value, 1));
  }

  inline void set_attribute_
  (const String name, const String attrName, const Num value) noexcept {
    assert_equal<Num, double>();
    herror(H5LTset_attribute_double
           (fileId_, name.c_str(), attrName.c_str(), &value, 1));
  }

  inline void get_attribute_
  (const String name, const String attrName, Ind& value) noexcept {
    assert_equal<Ind, long long>();
    herror(H5LTget_attribute_long_long
           (fileId_, name.c_str(), attrName.c_str(), &value));
  }

  inline void get_attribute_
  (const String name, const String attrName, SInd& value) noexcept {
    assert_equal<SInd, int>();
    herror(H5LTget_attribute_int
           (fileId_, name.c_str(), attrName.c_str(), &value));
  }

  inline void get_attribute_
  (const String name, const String attrName, Num& value) noexcept {
    assert_equal<Num, double>();
    herror(H5LTget_attribute_double
           (fileId_, name.c_str(), attrName.c_str(), &value));
  }

  std::tuple<Ind, Ind> get_data_dimensions(const String name) {
    Ind no_rows, no_cols;
    get_attribute(name, "no_rows", no_rows);
    get_attribute(name, "no_cols", no_cols);
    return std::make_tuple(no_rows, no_cols);
  }

  ///@}

  /// \brief Open/close files
  ///@{

  /// \brief Is a file open or closed?
  inline bool is_open() const noexcept { return isOpen_; }
  bool isOpen_;  ///< Is a file open or closed?

  /// \brief Opens file in read-only mode
  void open(access_mode::read) {
    assert_closed();
    if (!exists(fileName_, comm_)) {  // might be an user input error: no assert
      TERMINATE("file doesn't exist!");
    }
    fileId_ = H5Fopen(fileName_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    isOpen_ = true;
  }
  /// \brief Creates a file and opens it in read-write mode
  void open(access_mode::write) {
    assert_closed();
    if (!exists(fileName_, comm_)) {
      fileId_ = H5Fcreate(fileName_.c_str(), H5F_ACC_EXCL,
                          H5P_DEFAULT, H5P_DEFAULT);
      isOpen_ = true;
    } else {  // file exists: append mode
      fileId_ = H5Fopen(fileName_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      isOpen_ = true;
    }
    ASSERT(exists(fileName_, comm_), "file doesn't exist!");
  }

  /// \brief Closes file
  void close() {
    assert_open();
    herror(H5Fclose(fileId_));
    isOpen_ = false;
  }
  ///@}

  /// \brief Handles errors.
  void handle_error(herr_t status, const String position) noexcept {
    if (status < 0) {
      std::cerr << "There was an error in HDF5 at positon:\n"
                << position << "\n Printing the HDF5 Error stack:\n";
      H5Eprint(H5E_DEFAULT, stderr);
      TERMINATE("HDF5 error, see above!");
    }
  }

  Ind array_order_to_value(hdf5::array_order order) const noexcept {
    switch (order) {
      case hdf5::array_order::row_major_order: { return 0; }
      case hdf5::array_order::col_major_order: { return 1; }
      default: { TERMINATE("Unknown array order"); }
    }
  }

  hdf5::array_order get_data_order(const String name) {
    Ind order;
    get_attribute(name, "order", order);
    return array_value_to_order(order);
  }

  hdf5::array_order array_value_to_order(Ind order) const noexcept {
    switch (order) {
      case 0: { return hdf5::array_order::row_major_order; }
      case 1: { return hdf5::array_order::col_major_order; }
      default: { TERMINATE("Unknown array order"); }
    }
  }

  template<class T = Num>
  inline memory::buffer::Adapt<T> buffer(const Ind noTs, T) {
    return memory::buffer::Adapt<T>{buffer_, noTs};
  }

  hid_t fileId_;
  memory::Buffer buffer_;
};

using HDF5w = HDF5File<access_mode::write>;
using HDF5r = HDF5File<access_mode::read>;

////////////////////////////////////////////////////////////////////////////////
}  // namespace file
}  // namespace io
}  // namespace hom3
////////////////////////////////////////////////////////////////////////////////
#undef assert_open
#undef assert_closed
#undef herror
#undef assert_read_mode
#undef assert_write_mode
////////////////////////////////////////////////////////////////////////////////
#endif

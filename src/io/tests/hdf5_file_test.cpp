/// \file \brief Tests for HDF5 File handlers
#include <random>
#include "gtest/gtest.h"
#include "globals.hpp"
#include "io/hdf5_file.hpp"
////////////////////////////////////////////////////////////////////////////////
using namespace hom3; using namespace io; using namespace file;

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    mpi::environment env(argc, argv);
    return RUN_ALL_TESTS();
}

auto initialize(const String fName) {
  auto comm = boost::mpi::communicator();
  if (exists(fName, comm)) { remove(fName, comm); }
  return comm;
}

TEST(hdf5_file_test, filesystem_utilities) {
  auto p = boost::filesystem::current_path();  // This shouldn't produce garbage
  std::cout << p << std::endl;
}

/// \test Opening modes and fundamental filesystem functionality (exists/remove)
TEST(hdf5_file_test, opening_modes) {
  String fName{"test" + hdf5::extension()};
  auto comm = initialize(fName);

  ASSERT_FALSE(exists(fName, comm));
  // Opening a non-existent file in write mode creates a file:
  {
    auto file = HDF5w{fName, comm};
  }
  ASSERT_TRUE(exists(fName, comm));

  // Opening an existing file in read only mode:
  {
    auto file = HDF5r{fName, comm};
  }
  ASSERT_TRUE(exists(fName, comm));

  // Opening an existing file in write mode allows appending data to the file:
  {
    auto file = HDF5w{fName, comm};
  }
  ASSERT_TRUE(exists(fName, comm));

  remove(fName, comm);
  ASSERT_FALSE(exists(fName, comm));

  // Opening non-existing file in read only mode is not allowed!
  EXPECT_DEATH((HDF5r{fName, comm}), "[\\S\\s]+");
}

template<class T>
void check_equal_datasets(const std::vector<T>& v1, const std::vector<T>& v2) {
  for (auto i : zip(v1, v2)) { EXPECT_EQ(boost::get<0>(i), boost::get<1>(i)); }
}

void check_equal_datasets
(const std::vector<Num>& v1, const std::vector<Num>& v2) {
  for (auto i : zip(v1, v2)) {
    EXPECT_FLOAT_EQ(boost::get<0>(i), boost::get<1>(i));
  }
}

template<class T> void check_equal_attributes(const T& a, const T& b)
{EXPECT_EQ(a, b); }

void check_equal_attributes(const Num& a, const Num& b)
{ EXPECT_FLOAT_EQ(a, b); }

template<class T> void basic_read_write_memory() {
  String fName{"test" + hdf5::extension()};
  auto comm = initialize(fName);

  std::vector<T> original = {11, 22, 33, 44, 55, 66, 77, 88, 99};
  T originalRootAttr = T{33};
  T originalDataSetAttr = T{44};

  /// Write to file:
  {
    auto file = HDF5w{fName, comm};

    /// Write data to file:
    file.array("/dset", original.data(), original.size());

    /// Write attributes:
    file.attribute("/", "test_root_attr", T{originalRootAttr});
    file.attribute("/dset", "test_ds_attr", T{originalDataSetAttr});
    EXPECT_DEATH((file.attribute("/non_existent_data_set", "attr",
                                 T{originalDataSetAttr})), "[\\S\\s]+");
  }

  /// Read from file:
  std::vector<T> retrieved(9);
  {
    auto file = HDF5r{fName, comm};

    /// Read data from file:
    file.array("/dset", retrieved.data());

    /// Check data-set attributes
    Ind no_rows, no_cols;
    file.attribute("/dset", "no_rows", no_rows);
    file.attribute("/dset", "no_cols", no_cols);
    check_equal_attributes(9, no_rows);
    check_equal_attributes(1, no_cols);

    /// Read attributes:
    T rootAtt;
    file.attribute("/", "test_root_attr", rootAtt);
    check_equal_attributes(originalRootAttr, rootAtt);
    T dsAtt;
    file.attribute("/dset", "test_ds_attr", dsAtt);
    check_equal_attributes(originalDataSetAttr, dsAtt);
    EXPECT_DEATH((file.attribute("/dset", "non_existent_attr", dsAtt)),
                 "[\\S\\s]+");
    EXPECT_DEATH((file.attribute("/non_existent_dataset", "attr", dsAtt)),
                 "[\\S\\s]+");
  }

  check_equal_datasets(original, retrieved);
}

/// \test Writes directly from memory to file, reads the written file, and
/// compares it with initial memory and test attributes.
TEST(hdf5_file_test, basic_read_write_memory) {
  basic_read_write_memory<Ind> ();
  basic_read_write_memory<SInd>();
  basic_read_write_memory<Int> ();
  basic_read_write_memory<SInt>();
  basic_read_write_memory<Num> ();
}

/// \test Writes memory from complex data-structure using a functor + range,
/// reads the memory back to another complex data-structure, and compares both
/// of them
TEST(hdf5_file_test, complex_read_write_rt) {
  String fName{"test" + hdf5::extension()};
  auto comm = initialize(fName);

  /// Data-structure: Range of elements + two arrays + functor
  Range<Ind> row_range(0, 6), col_range(0, 2);
  std::array<Ind, 6> data0 = {{ 1, 3, 5, 7,  9, 11 }};
  std::array<Ind, 6> data1 = {{ 2, 4, 6, 8, 10, 12 }};
  auto f = [&](const Ind rowIdx, const Ind colIdx) {
    if (colIdx == 0) {
      return data0[rowIdx];
    } else if (colIdx == 1) {
      return data1[rowIdx];
    } else {
      TERMINATE("error!");
    }
  };

  // Write to file:
  {
    auto file = HDF5w{fName, comm};
    file.array("/dset", f, row_range, col_range, hdf5::order::row_major);
  }

  /// Read file to 1D data-structure:
  std::array<Ind, 12> target;
  {
    auto file = HDF5r{fName, comm};
    file.array("/dset", target.data());

    /// Check data-set attributes:
    Ind no_rows, no_cols;
    file.attribute("/dset", "no_rows", no_rows);
    file.attribute("/dset", "no_cols", no_cols);
    check_equal_attributes(6, no_rows);
    check_equal_attributes(2, no_cols);
  }

  for (Ind i = 0, j = 0, k = 0, e = target.size(); i < e; ++i) {
    if (i % 2 == 0) {
      EXPECT_EQ(target[i], data0[j]);
      ++j;
    } else {
      EXPECT_EQ(target[i], data1[k]);
      ++k;
    }
  }

  /// Read file back to another complex data-structure
  std::array<Ind, 6> data0r, data1r;
  auto fr = [&](Ind i, Ind d) mutable -> Ind& {
    if (d == 0) {
      return data0r[i];
    } else if (d == 1) {
      return data1r[i];
    } else {
      TERMINATE("error!");
    }
  };

  {
    auto file = HDF5r{fName, comm};
    file.array("/dset", fr);
  }

  for (Ind i = 0, e = 6; i < e; ++i) {
     EXPECT_EQ(data0[i], data0r[i]);
     EXPECT_EQ(data1[i], data1r[i]);
  }
}

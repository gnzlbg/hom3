/// \file \brief Implements tests for MPI related functionality
#include "misc/test.hpp"
#include "globals.hpp"
////////////////////////////////////////////////////////////////////////////////

using namespace hom3;

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    mpi::environment env(argc, argv);
    return RUN_ALL_TESTS();
}

TEST(mpi_test, initialization) {
  using namespace mpi;
  communicator world;

  if (is_root(world)) {
    String msg, out_msg = "Hello";
    requests<2> reqs = {{
      world.isend(1, 0, out_msg),
      world.irecv(1, 1, msg)
    }};
    wait_all(reqs);
    std::cout << msg << "!" << std::endl;
  } else {
    String msg, out_msg = "world";
    requests<2> reqs = {{
      world.isend(0, 1, out_msg),
      world.irecv(0, 0, msg)
    }};
    wait_all(reqs);
    std::cout << msg << ", ";
  }
}

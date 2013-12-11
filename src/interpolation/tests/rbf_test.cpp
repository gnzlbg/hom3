#include <iostream>
#include <fstream>
#include "globals.hpp"
#include "misc/test.hpp"
#include "interpolation/interpolation.hpp"
////////////////////////////////////////////////////////////////////////////////
using namespace hom3;

static const int nd = 2;

TEST(interpolation_rbf, first_test) {
  namespace rip = interpolation::rbf;
  std::vector<NumA<2>> xs = { NumA<2>::Constant(0.),
                              NumA<2>::Constant(0.25),
                              NumA<2>::Constant(0.5),
                              NumA<2>::Constant(0.75),
                              NumA<2>::Constant(1.) };
  std::vector<Num> vs = {0. , 0.25, 0.5, 0.75, 1. };

  auto weights = rip::build_weights(xs, vs);
  std::cerr << "rows: " << weights.rows() << " cols: " << weights.cols() << "\n";

  {
    std::ofstream orig;
    orig.open ("orig.dat");
    for(std::size_t i = 0; i < xs.size(); ++i) {
      for(SInd d = 0; d < nd; ++d) {
        orig << xs[i](d) << "    ";
      }
      orig << vs[i] << "\n";
    }
  }
  {
    std::ofstream interp;
    interp.open ("interp.dat");
    Num start = 0., end = 1.001, step = 0.05;
    while(start < end) {
      NumA<2> xs_i = NumA<2>::Constant(start);
      Num vs_i = rip::interpolate(xs_i, xs, weights);
      for(SInd d = 0; d < nd; ++d) {
        interp << xs_i(d) << "    ";
      }
      interp << vs_i << "\n";
      start += step;
    }
  }
}

TEST(interpolation_rbf, second_test) {
  namespace rip = interpolation::rbf;
  NumA<nd> p0 = NumA<nd>::Constant(0.0);
  NumA<nd> p1 = {1.0, 0.0};
  NumA<nd> p2 = {0.0, 1.0};
  NumA<nd> p3 = NumA<nd>::Constant(1.0);

  const std::vector<NumA<nd>> xs = {p0, p1, p2, p3};
  const std::vector<Num> vs = {-1.0, 1.0, 0.5, -0.5};

  auto weights = rip::build_weights(xs, vs);

  {
    std::ofstream orig;
    orig.open ("orig2D.dat");
    for(std::size_t i = 0; i < xs.size(); ++i) {
      for(SInd d = 0; d < nd; ++d) {
        orig << xs[i](d) << "    ";
      }
      orig << vs[i] << "\n";
    }
  }

  {
    std::ofstream interp;
    interp.open ("interp2D.dat");
    Num x = 0., x_end = 1.001, x_step = 0.05;
    Num y = 0., y_end = 1.001, y_step = 0.05;
    while(x < x_end) {
      while(y < y_end) {
        NumA<2> xs_i = {x, y};
        Num vs_i = rip::interpolate(xs_i, xs, weights);
        for(SInd d = 0; d < nd; ++d) {
          interp << xs_i(d) << "    ";
        }
        interp << vs_i << "\n";
        y += y_step;
      }
      interp << "\n";
      x += x_step;
      y = 0.;
    }
  }
}

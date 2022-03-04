#include "runtime/expression_tree.hpp"
#include "runtime/sollya_fix_function.hpp"
#include "runtime/sollya_handler.hpp"
#include <cmath>
#include <iostream>
#include <sollya.h>

using namespace archgenlib;

int main() {
  SollyaHandler sinus{SOLLYA_SIN(SOLLYA_X_)};
  FPDimRTRepr repr{5, -4, false};
  SollyaFunction sf{sinus, repr};
  auto reprvec = sf.faithful_at_weight(-4);

  constexpr float inc{0x1p-4};
  static_assert(inc * 16 == 1.);

  float cur_val = 0.;
  for (auto &val : reprvec) {
    double libm_res = std::sin(cur_val);
    double agl{val.get_d() / 16.};
    if (std::abs(libm_res - agl) >= inc) {
      std::cerr << "Err for " << cur_val << ": libm = " << libm_res
                << ", archgenlib = " << agl << "\n";
      return -1;
    }
    cur_val += inc; //Should be exact
  }
  std::cout << "All computed values seems to be faithful" << std::endl;
  return 0;
}

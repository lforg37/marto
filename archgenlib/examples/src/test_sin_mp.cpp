#include "runtime/expression_tree.hpp"
#include "runtime/fixfunction_multipartite.hpp"
#include "runtime/sollya_fix_function.hpp"
#include "runtime/sollya_handler.hpp"
#include <cmath>
#include <iostream>
#include <sollya.h>

using namespace archgenlib;

int main() {
  SollyaHandler sinus{SOLLYA_SIN(SOLLYA_X_)};
  FixedFormatRTRepr repr{0, -8, false};
  SollyaFunction sf{sinus, repr};


  MultipartiteFunction mpf{sf, -8};
  if (mpf.check_best_config(-16)) {
    std::cout << "All computed values seems to be faithful" << std::endl;
    return 0;
  } else {
    std::cerr << "Error !\n";
    return -1;
  }
}

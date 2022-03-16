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
  FPDimRTRepr repr{0, -8, false};
  SollyaFunction sf{sinus, repr};

  mpz_class a{-32};

  a >>= 3;
  assert(a == -4);

  MultipartiteFunction mpf{sf, -8};
  if (mpf.check_best_config(-37)) {
    std::cout << "All computed values seems to be faithful" << std::endl;
    return 0;
  } else {
    std::cerr << "Error !\n";
    return 0;
  }
}

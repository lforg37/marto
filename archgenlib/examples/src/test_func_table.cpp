#include <iostream> 
#include <sollya.h>
#include "runtime/expression_tree.hpp"
#include "runtime/sollya_fix_function.hpp"
#include "runtime/sollya_handler.hpp"

using namespace archgenlib;

int main() {
  SollyaHandler identity{sollya_lib_build_function_free_variable()};
  FPDimRTRepr repr{5, -4, true};
  SollyaFunction sf{identity, repr};
  auto reprvec = sf.faithful_at_weight(-4);
  for (auto & val : reprvec) {
    std::cout << val << "\n";
  }
  return 0;
}

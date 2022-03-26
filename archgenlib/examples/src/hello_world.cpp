#include <iostream>

#include <sollya.h>

#include "runtime/sollya_handler.hpp"
#include "runtime/output_formatter.hpp"
#include "fixedpoint/fixedpoint.hpp"

int main() {
  // No handler as the sollya_lib_build_* take ownership of the arguments.
  auto a = sollya_lib_constant_from_int(42); 
  auto b = archgenlib::SollyaHandler{sollya_lib_build_function_sqrt(a)};
  double res;
  sollya_lib_get_constant_as_double(&res, b);
  using constant_t = archgenlib::FixedConstant<archgenlib::FixedFormat<14, -23, unsigned>, 0x1111>;
  std::cout << archgenlib::detail::type_name<constant_t>() << std::endl;
  return 0;
}

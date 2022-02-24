#include "sollya_handler.hpp"
#include <sollya.h>

int main() {
  archgenlib::SollyaHandler a{sollya_lib_constant_from_int(42)};
  return 0;
}

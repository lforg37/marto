#include <cassert>

#include "bitint_tools/type_helpers.hpp"
#include "fixedpoint/fixedpoint.hpp"

using namespace archgenlib;

int main() {
  FixedFormat<4, -2, unsigned> format{};
  hint::detail::bitint_base_t<false, 7> val{0b01010};
  FixedNumber op1(format, {0b01010});
  FixedNumber<FixedFormat<3, -5, unsigned>>op2({0b11111});

  //      01010
  // +       11111
  // ---------------
  //     001101111

  FixedNumber<FixedFormat<5, -5, unsigned>> res({0b1101111});
  assert(res == (op1 + op2));
  return 0;
}

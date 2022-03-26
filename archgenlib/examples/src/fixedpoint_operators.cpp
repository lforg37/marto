#include "bitint_tools/type_helpers.hpp"
#include "fixedpoint/fixedpoint.hpp"

using namespace archgenlib;
int main() {
  constexpr FixedFormat<4, -2, unsigned> format{};
  constexpr FixedNumber op1(format, {0b01010});
  constexpr FixedNumber<FixedFormat<3, -5, unsigned>>op2({0b11111});

  //      01010
  // +       11111
  // ---------------
  //     001101111

  constexpr FixedNumber<FixedFormat<5, -5, unsigned>> res_add({0b1101111});
  static_assert(res_add == (op1 + op2));

  /*
        0101
    * 11111
    --------
        11111
    +  1111100
    -----------
    10011011
  */
  constexpr FixedNumber<FixedFormat<8, -7, unsigned>> res_mul{0b100110110};
  static_assert(res_mul == (op1 * op2));
  return 0;
}

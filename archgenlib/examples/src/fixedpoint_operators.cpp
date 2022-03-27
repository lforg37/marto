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


  constexpr FixedNumber<FixedFormat<5, -3, signed>> pattern {0b101101110};
  constexpr auto high_bits = pattern.template extract<5, 1>();
  static_assert(high_bits == FixedNumber<FixedFormat<5, 1, signed>>{0b10110});
  static_assert(pattern.template extract<0, -3>() == FixedNumber<FixedFormat<0, -3, unsigned>>{0b1110});

  // Conversions from int // Overlap
  using format_0 = FixedFormat<8, -2, signed>;
  using fpnum_0 = FixedNumber<format_0>;
  static_assert(fpnum_0::get_from_value(15) == fpnum_0{0b111100});
  static_assert(fpnum_0::get_from_value(256) == fpnum_0{0b1111111111});
  static_assert(fpnum_0::get_from_value(-1) ==  fpnum_0{0b11111111100});
  static_assert(fpnum_0::get_from_value(-256) ==  fpnum_0{0b10000000000});
  static_assert(fpnum_0::get_from_value(-272) ==  fpnum_0{0b10000000000});

  using format_1 = FixedFormat<-1, -4, unsigned>;
  using fpnum_1 = FixedNumber<format_1>;
  static_assert(fpnum_1::get_from_value(1) == fpnum_1{0b1111});
  static_assert(fpnum_1::get_from_value(256) == fpnum_1{0b1111});
  static_assert(fpnum_1::get_from_value(0) == fpnum_1{0b0});
  static_assert(fpnum_1::get_from_value(-37) == fpnum_1{0b0});
  static_assert(fpnum_1::get_from_value(-1) ==  fpnum_1{0b0});

  using format_2 = FixedFormat<-1, -4, signed>;
  using fpnum_2 = FixedNumber<format_2>;
  static_assert(fpnum_2::get_from_value(1) == fpnum_2{0b0111});
  static_assert(fpnum_2::get_from_value(256) == fpnum_2{0b0111});
  static_assert(fpnum_2::get_from_value(0) == fpnum_2{0b0});
  static_assert(fpnum_2::get_from_value(-37) == fpnum_2{0b1000});
  static_assert(fpnum_2::get_from_value(-1) ==  fpnum_2{0b1000});

  using format_3 = FixedFormat<8, 3, unsigned>;
  using fpnum_3 = FixedNumber<format_3>;
  static_assert(fpnum_3::get_from_value(0b1000) == fpnum_3{0b1});
  //Ties toward 0
  static_assert(fpnum_3::get_from_value(0b1100) == fpnum_3{0b1});
  static_assert(fpnum_3::get_from_value(0b1101) == fpnum_3{0b10});

  static_assert(fpnum_3::get_from_value(0b11111111111101) == fpnum_3{0b111111});

  return 0;
}

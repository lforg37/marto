#include "bitint_tools/type_helpers.hpp"
#include "fixedpoint/fixedpoint.hpp"
#include <cstdint>
#include <limits>

using namespace archgenlib;
int main() {
  constexpr FixedFormat<4, -2, unsigned> format{};
  constexpr FixedNumber op1(format, {0b01010});
  constexpr FixedNumber<FixedFormat<3, -5, unsigned>> op2({0b11111});

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

  constexpr FixedNumber<FixedFormat<5, -3, signed>> pattern{0b101101110};
  constexpr auto high_bits = pattern.template extract<5, 1>();
  static_assert(high_bits == FixedNumber<FixedFormat<5, 1, signed>>{0b10110});
  static_assert(pattern.template extract<0, -3>() ==
                FixedNumber<FixedFormat<0, -3, unsigned>>{0b1110});

  // Conversions from int // Overlap
  using format_0 = FixedFormat<7, -2, signed>;
  using fpnum_0 = FixedNumber<format_0>;
  static_assert(fpnum_0::get_from_value(15) == fpnum_0{0b111100});
  static_assert(fpnum_0::get_from_value(256) == fpnum_0{0b111111111});
  static_assert(fpnum_0::get_from_value(-1) == fpnum_0{0b1111111100});
  static_assert(fpnum_0::get_from_value(-256) == fpnum_0{0b1000000000});
  static_assert(fpnum_0::get_from_value(-272) == fpnum_0{0b1000000000});

  static_assert(fpnum_0{0b111100}.get_as<int>() == 15);
  static_assert(fpnum_0{0b111110}.get_as<int>() == 15);
  static_assert(fpnum_0{0b111111}.get_as<int>() == 16);

  static_assert(fpnum_0{0b11111111111}.get_as<int>() == 0);
  static_assert(fpnum_0{0b11111111110}.get_as<int>() == -1);
  static_assert(fpnum_0{0b11111111100}.get_as<int>() == -1);

  static_assert(fpnum_0{-16}.get_as<int>() == -4);
  static_assert(fpnum_0{-1 << 9}.get_as<std::uint8_t>() == 0);
  static_assert(fpnum_0{-1 << 9}.get_as<std::int8_t>() ==
                std::numeric_limits<std::int8_t>::min());

  using format_1 = FixedFormat<-1, -4, unsigned>;
  using fpnum_1 = FixedNumber<format_1>;
  static_assert(fpnum_1::get_from_value(1) == fpnum_1{0b1111});
  static_assert(fpnum_1::get_from_value(256) == fpnum_1{0b1111});
  static_assert(fpnum_1::get_from_value(0) == fpnum_1{0b0});
  static_assert(fpnum_1::get_from_value(-37) == fpnum_1{0b0});
  static_assert(fpnum_1::get_from_value(-1) == fpnum_1{0b0});
  static_assert(fpnum_1{0b1001}.get_as<int>() == 1);
  static_assert(fpnum_1{0b1000}.get_as<int>() == 0);
  static_assert(fpnum_1{0b0001}.get_as<int>() == 0);

  using format_2 = FixedFormat<-1, -4, signed>;
  using fpnum_2 = FixedNumber<format_2>;
  static_assert(fpnum_2::get_from_value(1) == fpnum_2{0b0111});
  static_assert(fpnum_2::get_from_value(256) == fpnum_2{0b0111});
  static_assert(fpnum_2::get_from_value(0) == fpnum_2{0b0});
  static_assert(fpnum_2::get_from_value(-37) == fpnum_2{0b1000});
  static_assert(fpnum_2::get_from_value(-1) == fpnum_2{0b1000});

  static_assert(fpnum_2{0b1000}.get_as<int>() == 0);

  using format_3 = FixedFormat<8, 3, unsigned>;
  using fpnum_3 = FixedNumber<format_3>;
  static_assert(fpnum_3::get_from_value(0b1000) == fpnum_3{0b1});
  // Ties toward 0
  static_assert(fpnum_3::get_from_value(0b1100) == fpnum_3{0b10});
  static_assert(fpnum_3::get_from_value(0b1101) == fpnum_3{0b10});

  static_assert(fpnum_3::get_from_value(0b11111111111101) == fpnum_3{0b111111});

  static_assert(fpnum_3{0b11}.get_as<int>() == 3 << 3);
  static_assert(fpnum_3{0b100000}.get_as<std::uint8_t>() ==
                std::numeric_limits<std::uint8_t>::max());

  // Floats
  auto from_val_0 = fpnum_0::get_from_value(2.);
  assert(from_val_0 == fpnum_0{0b1000});
  assert(from_val_0.get_as<double>() == 2.);
  auto from_val_1 = fpnum_0::get_from_value(.5);
  assert(from_val_1 == fpnum_0{0b10});
  assert(from_val_1.get_as<double>() == .5);

  auto from_val_2 = fpnum_0::get_from_value(3. / 16.);
  // round up
  assert(from_val_2 == fpnum_0{0b01});
  assert(from_val_2.get_as<double>() == .25);

  auto from_val_3 = fpnum_0::get_from_value(static_cast<double>(0b101101) / 16);
  assert(from_val_3 == fpnum_0{0b1011});

  auto from_val_4 = fpnum_0::get_from_value(static_cast<double>(0b101110) / 16);
  assert(from_val_4 == fpnum_0{0b1100});

  auto from_val_5 = fpnum_0::get_from_value(static_cast<double>(0b101111) / 16);
  assert(from_val_5 == fpnum_0{0b1100});

  auto from_val_6 = fpnum_0::get_from_value(static_cast<double>(0b1) / 64);
  assert(from_val_6 == fpnum_0{0b0});

  auto from_val_7 =
      fpnum_0::get_from_value(static_cast<double>(0b101011010) / 4);
  assert(from_val_7 == fpnum_0{0b101011010});

  auto from_val_8 = fpnum_0::get_from_value(static_cast<double>(0b111010110101) /
                                            2); // Overflow
  assert(from_val_8 == fpnum_0{0b111111111});

  // TODO check inf and negative numbers... Not an emergency



  return 0;
}

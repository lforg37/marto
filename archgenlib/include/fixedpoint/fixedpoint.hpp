#ifndef FIXEDPOINT_FIXEDPOINT_HPP
#define FIXEDPOINT_FIXEDPOINT_HPP

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <string_view>
#include <type_traits>

#ifndef BITINT_BACKEND
#define BITINT_BACKEND
#endif
#include "bitint_tools/bitint_constant.hpp"
#include "bitint_tools/type_helpers.hpp"
#include "hint.hpp"
#include "tools/static_math.hpp"

namespace archgenlib {
using bitweight_t = std::int32_t;
using vecwidth_t = std::uint32_t;

template<bool>
struct is_signed_to_type {
  using type = unsigned;
};

template<>
struct is_signed_to_type<true> {
  using type = signed;
};

template<bool b>
using is_signed_to_type_t = typename is_signed_to_type<b>::type;

template<bool b>
using is_unsigned_to_type_t = typename is_signed_to_type<!b>::type;

template<typename T>
using same_sign_t = is_signed_to_type_t<std::is_signed_v<T>>;

template <bitweight_t MSBWeight, bitweight_t LSBWeight, typename sign_t>
struct FixedFormat {
  static constexpr bool IsUnsigned = std::is_same<unsigned, sign_t>::value;
  static constexpr bool IsSigned = std::is_same<signed, sign_t>::value;
  static_assert(IsUnsigned || IsSigned,
                "invalid sign type expected signed or unsigned");
  static_assert(MSBWeight >= LSBWeight,
                "MSB weight cannot be smaller than LSBs");

  /**
   * @brief msb_weight Weight of the most significant bit
   */
  static constexpr bitweight_t msb_weight = MSBWeight;

  /**
   * @brief lsb_weight
   */
  static constexpr bitweight_t lsb_weight = LSBWeight;

  /**
   * @brief WFF Full significand width
   */
  static constexpr vecwidth_t width = MSBWeight - LSBWeight + 1;

  /**
   * @brief is_signed is the number signed
   */
  static constexpr bool is_signed = IsSigned;

  /**
   * @brief max_positive_bitweight maximum bitweight for which associated bit is
   * positive.
   *
   */
  static constexpr bitweight_t max_positive_bitweight =
      is_signed ? msb_weight - 1 : msb_weight;

  using bitint_type = hint::detail::bitint_base_t<is_signed, width>;

  constexpr auto get_bit_int(auto const & val) const {
    return static_cast<bitint_type>(val);
  }
};

namespace detail {
template <typename T> constexpr bool is_fixed_format = false;

template <bitweight_t MSBWeight, bitweight_t LSBWeight, typename T>
constexpr bool is_fixed_format<FixedFormat<MSBWeight, LSBWeight, T>> = true;
} // namespace detail

template <typename T>
concept FixedFormatType = detail::is_fixed_format<T>;

namespace detail {
template <FixedFormatType T1, FixedFormatType T2>
constexpr bool can_extend_to(T1 source = {}, T2 dest = {}) {
  if (T1::msb_weight > T2::msb_weight || T1::lsb_weight < T2::lsb_weight ||
      T1::is_signed && !T2::is_signed) {
    return false;
  }
  // We can extend an unsigned value to signed value, but only if the signed
  // value msb weight is strictly greater
  if (!T1::is_signed && T2::is_signed && T1::msb_weight == T2::msb_weight) {
    return false;
  }
  return true;
}
}; // namespace detail

template <FixedFormatType FF1, FixedFormatType FF2>
constexpr auto operator+(FF1 format1, FF2 format2) {
  constexpr auto lsb_out = std::min(format1.lsb_weight, format2.lsb_weight);
  constexpr auto max_msb = std::max(format1.msb_weight, format2.msb_weight);
  constexpr auto max_pos_msb =
      std::max(format1.max_positive_bitweight, format2.max_positive_bitweight);
  constexpr auto one_signed = format1.is_signed || format2.is_signed;
  constexpr auto msb =
      1 + (((max_msb == max_pos_msb) && one_signed) ? max_msb + 1 : max_msb);
  return FixedFormat<msb, lsb_out, is_signed_to_type_t<one_signed>>{};
}

template <FixedFormatType FF1, FixedFormatType FF2>
constexpr auto operator*(FF1 format1, FF2 format2) {
  using arith_prop = hint::Arithmetic_Prop<FF1::width, FF2::width, FF1::is_signed, FF2::is_signed>;
  constexpr auto lsb_out = FF1::lsb_weight + FF2::lsb_weight;
  constexpr auto msb_out = lsb_out + arith_prop::_prodSize - 1;
  return FixedFormat<msb_out, lsb_out, is_signed_to_type_t<arith_prop::_prodSigned>>{};
}

template <std::integral IT>
using fixedformat_from_integral =
    FixedFormat<std::numeric_limits<IT>::digits + std::is_signed_v<IT>, 0,
          same_sign_t<IT>>;

template <FixedFormatType Format> class FixedNumber {
public:
  using storage_t =
      hint::detail::bitint_base_t<Format::is_signed, Format::width>;
  storage_t const value_;

  using format_t = Format;
  static constexpr auto width = Format::width;
  constexpr FixedNumber(storage_t const &val) : value_{val} {}

  template <FixedFormatType FT>
  constexpr FixedNumber(FT, storage_t const &val) : value_{val} {}

  explicit constexpr FixedNumber(hint::BitIntWrapper<Format::width, Format::is_signed> const &val)
      : value_{val.unravel()} {}
  constexpr storage_t value() const { return value_; }
  hint::BitIntWrapper<Format::width, Format::is_signed> as_hint() const {
    return {value_};
  }

  static constexpr format_t format{};

  template <FixedFormatType NewFormat>
  constexpr FixedNumber<NewFormat> extend_to(NewFormat new_format = {}) const {
    static_assert(detail::can_extend_to(format, new_format),
                  "Trying to extend value to a format that is not a superset "
                  "of current format");
    auto new_value =
        static_cast<typename FixedNumber<NewFormat>::storage_t>(value_);
    return {new_value << (format.lsb_weight - new_format.lsb_weight)};
  }
};

template <FixedFormatType FT> FixedNumber(FT, auto) -> FixedNumber<FT>;

template <FixedFormatType FT, typename T>
FixedNumber(FT, std::initializer_list<T>) -> FixedNumber<FT>;

namespace detail {
template <typename T> constexpr bool is_fixed_num = false;

template <FixedFormatType Dim>
constexpr bool is_fixed_num<FixedNumber<Dim>> = true;
} // namespace detail

template <typename T>
concept FixedNumberType = detail::is_fixed_num<T>;

// FixedPointConstant
template <typename T, typename Dim>
concept FixedBitIntCst =
    hint::BitIntConstT<T> && FixedFormatType<Dim> && T::width ==
Dim::width &&T::isSigned == Dim::is_signed;

template <FixedFormatType Dim,
          hint::detail::bitint_base_t<Dim::is_signed, Dim::width> Val>
struct FixedConstant {
  using dimension_t = Dim;
  using value_type = hint::detail::bitint_base_t<Dim::is_signed, Dim::width>;
  static constexpr value_type value = Val;
};

namespace detail {
template <typename T> constexpr bool _is_fixed_constant = false;
template <FixedFormatType Dim,
          hint::detail::bitint_base_t<Dim::is_signed, Dim::width> Val>
constexpr bool _is_fixed_constant<FixedConstant<Dim, Val>> = true;
} // namespace detail

template <typename T>
concept FixedConstantType = detail::_is_fixed_constant<T>;

template <FixedNumberType T1, FixedNumberType T2>
constexpr auto operator+(T1 const &op1, T2 const &op2) {
  constexpr auto res_format = T1::format + T2::format;
  auto resized_1 = op1.extend_to(res_format);
  auto resized_2 = op2.extend_to(res_format);
  auto res_val = resized_1.value() + resized_2.value();
  return FixedNumber(res_format, res_val);
}

template <FixedNumberType T1, FixedNumberType T2>
constexpr auto operator-(T1 const &op1, T2 const &op2) {
  constexpr auto res_format = T1::format + T2::format;
  auto resized_1 = op1.extend_to(res_format);
  auto resized_2 = op2.extend_to(res_format);
  auto res_val = resized_1.value() - resized_2.value();
  return FixedNumber(res_format, res_val);
}

template <FixedNumberType T1, FixedNumberType T2>
constexpr auto operator*(T1 const &op1, T2 const &op2) {
  constexpr auto res_format = T1::format * T2::format;
  auto val1 =  res_format.get_bit_int(op1.value());
  auto val2 = res_format.get_bit_int(op2.value());
  return FixedNumber(res_format, val1 * val2);
}

template <FixedNumberType FT> constexpr bool operator==(FT const &op1, FT const &op2) {
  return op1.value() == op2.value();
}

} // namespace archgenlib

#endif

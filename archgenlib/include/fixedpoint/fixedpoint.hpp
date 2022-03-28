#ifndef FIXEDPOINT_FIXEDPOINT_HPP
#define FIXEDPOINT_FIXEDPOINT_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
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

template <bool> struct is_signed_to_type {
  using type = unsigned;
};

template <> struct is_signed_to_type<true> {
  using type = signed;
};

template <bool b>
using is_signed_to_type_t = typename is_signed_to_type<b>::type;

template <bool b>
using is_unsigned_to_type_t = typename is_signed_to_type<!b>::type;

template <typename T>
using same_sign_t = is_signed_to_type_t<std::is_signed_v<T>>;

template <bitweight_t MSBWeight, bitweight_t LSBWeight, typename SignT>
struct FixedFormat {
  static constexpr bool IsUnsigned = std::is_same<unsigned, SignT>::value;
  static constexpr bool IsSigned = std::is_same<signed, SignT>::value;
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

  static_assert(width != 277, "");

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

  constexpr auto get_bit_int(auto const &val) const {
    return static_cast<bitint_type>(val);
  }

  using sign_t = SignT;
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
  if (T1::msb_weight > T2::msb_weight || T1::lsb_weight < T2::lsb_weight) {
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
  using arith_prop = hint::Arithmetic_Prop<FF1::width, FF2::width,
                                           FF1::is_signed, FF2::is_signed>;
  constexpr auto lsb_out = FF1::lsb_weight + FF2::lsb_weight;
  constexpr auto msb_out = lsb_out + arith_prop::_prodSize - 1;
  return FixedFormat<msb_out, lsb_out,
                     is_signed_to_type_t<arith_prop::_prodSigned>>{};
}

template <std::integral IT>
using fixedformat_from_integral =
    FixedFormat<std::numeric_limits<IT>::digits + std::is_signed_v<IT>, 0,
                same_sign_t<IT>>;

template <FixedFormatType Format> class FixedNumber {
public:
  template <FixedFormatType>
  friend class FixedNumber;

  using storage_t = typename Format::bitint_type;
private:
  storage_t value_;
public:

  using format_t = Format;
  static constexpr auto width = Format::width;
  constexpr FixedNumber(storage_t const &val = {0}) : value_{val} {}

  template <FixedFormatType FT>
  constexpr FixedNumber(FT, storage_t const &val = {0}) : value_{val} {}


  template <FixedFormatType FT>
  constexpr FixedNumber modular_add(FixedNumber<FT> const & val) {
    static_assert(FT::msb_weight <= format.msb_weight);
    static_assert(FT::lsb_weight >= format.lsb_weight);
    auto aligned = val.extend_to(format);
    storage_t res{value_ + aligned.value_};
    return res;
  }

  /**
   * @brief Performs accumulation operation
   * 
   * @tparam FT 
   * @param val 
   * @return constexpr FixedNumber& 
   */
  template <FixedFormatType FT>
  FixedNumber & operator+=(FixedNumber<FT> const & val) {
    value_ = modular_add(val)._value;
    return *this;
  }


  template <FixedFormatType FT>
  constexpr FixedNumber modular_sub(FixedNumber<FT> const & val) {
    static_assert(FT::msb_weight <= format.msb_weight);
    static_assert(FT::lsb_weight >= format.lsb_weight);
    auto aligned = val.extend_to(format);
    storage_t res{value_ - aligned.value_};
    return res;
  }

  template <FixedFormatType FT>
  constexpr FixedNumber & modular_mult(FixedNumber<FT> const & op) {
    constexpr auto prod_fmt = FT{} * format;
    static_assert(prod_fmt.msb_weight >= format.lsb_weight);
    static_assert(prod_fmt.lsb_weight <= format.msb_weight);
    auto prod = *this * op;
    auto res_val = prod.template extract<std::min(format.msb_weight, prod_fmt.msb_weight), std::max(format.lsb_weight, prod_fmt.lsb_weight)>();
    auto extended = res_val.extend_to(format);
    return value_ + extended.value_; 
  }

  template <FixedFormatType FT>
  FixedNumber & operator*=(FixedNumber<FT> const & op) {
    // Avoid cases were result is always zero
    value_ = modular_mult(op).value_;
    return *this;
  }

    /**
   * @brief Performs accumulation operation
   * 
   * @tparam FT 
   * @param val 
   * @return constexpr FixedNumber& 
   */
  template <FixedFormatType FT>
  FixedNumber & operator-=(FixedNumber<FT> const & val) {
    value_ = modular_sub(val)._value;
    return *this;
  }

  template <std::floating_point FT>
  /**
   * @brief Get the fixed number with the value closest to input
   *
   * @param in_val value to convert to FixedNumber
   * @return FixedNumber read above
   */
  static FixedNumber get_from_value(FT const in_val) {
    using limits = std::numeric_limits<FT>;
    static_assert(limits::is_iec559);
    static_assert(limits::radix == 2,
                  "The library only handles radix 2 floats.");
    assert(!std::isnan(in_val)); // No good value for NaN;
    if (in_val == 0) {
      return {0};
    } else if (in_val == limits::max()) {
      return max_val();
    } else if (in_val == limits::min()) {
      return min_val();
    }
    constexpr auto fp_width = limits::digits;
    int exp_get;
    auto normalized = std::frexp(in_val, &exp_get);
    int exp_msb = exp_get - 1; // Get the weight of the LSB
    int exp_lsb = exp_get - fp_width; // Get the weight of the LSB
    if (exp_msb > format.msb_weight ||
        (format.is_signed && exp_msb == format.msb_weight)) {
      if (in_val < 0) {
        return min_val();
      } else {
        return max_val();
      }
    }

    if (exp_msb < format.lsb_weight - 1) {
      return 0;
    }

    if (exp_msb == format.lsb_weight - 1) {
      if (normalized > 0) {
        return min_pos();
      } else {
        return min_neg();
      }
    }

    if (!format.is_signed && in_val < 0) {
      return {0};
    }

    // Now there should be at least some overlap;
    constexpr auto scaler_width = fp_width + format.width;
    using fp_storage_t =
        hint::detail::bitint_base_t<format.is_signed,
                                    fp_width + (format.is_signed)>;
    auto scaled = std::ldexp(normalized, fp_width);
    auto fp_bitint_val = static_cast<fp_storage_t>(scaled);
    constexpr auto extended_format_lsb_weight = format.lsb_weight - fp_width;
    using extended_fixed_dim =
        FixedFormat<format.msb_weight, extended_format_lsb_weight,
                    typename Format::sign_t>;
    using extended_storage = typename extended_fixed_dim::bitint_type;
    FixedNumber<extended_fixed_dim> extended{
        static_cast<extended_storage>(fp_bitint_val)
        << (exp_lsb - extended_format_lsb_weight)};

    return extended.template round_to<format.lsb_weight>();
  }

  template <std::integral IT>
  /**
   * @brief Get the fixed number with the value closest to input
   *
   * @param in_val value to convert to FixedNumber
   * @return FixedNumber read above
   */
  static constexpr FixedNumber get_from_value(IT const in_val) {
    using limits = std::numeric_limits<IT>;
    if constexpr (Format::msb_weight < 0) {
      // Only fractional bits
      if (in_val >= 1) {
        return max_val();
      } else if (in_val < 0) {
        return min_val();
      } else {
        return 0;
      }
    } else if constexpr (limits::max() < min_pos().value_) {
      // We have a lsb which is above the maximum representable value
      constexpr auto half_val =
          hint::detail::bitint_base_t<true, Format::lsb_weight - 1>{1}
          << (Format::lsb_weight - 2);
      if constexpr (Format::isSigned) {
        if (in_val < 0) {
          return (-in_val > half_val) ? -1 : 0;
          return;
        }
      }
      return (in_val > half_val) ? 1 : 0;
    } else {
      using full_int_format =
          FixedFormat<Format::msb_weight, 0, typename Format::sign_t>;

      constexpr auto get_shift = std::max(Format::lsb_weight, 0);
      constexpr auto min_val_as_int =
          min_val().template extract<Format::msb_weight, get_shift>().extend_to(
              full_int_format{});
      constexpr auto max_val_as_int =
          max_val().template extract<Format::msb_weight, get_shift>().extend_to(
              full_int_format{});
      if (in_val < min_val_as_int.value()) {
        return min_val();
      } else if (in_val > max_val_as_int.value()) {
        return max_val();
      } else {
        if constexpr (get_shift == 0) {
          auto val = static_cast<hint::detail::bitint_base_t<
              Format::is_signed, Format::msb_weight + 1>>(in_val);
          return FixedNumber<full_int_format>(val).extend_to(format).value();
        } else if constexpr (get_shift >= 1) {
          constexpr auto round_bit_mask = IT{1} << (get_shift - 1);
          using fui_storage_t = typename full_int_format::bitint_type;
          auto val = static_cast<fui_storage_t>((in_val >> get_shift));
          bool round_up = (in_val & round_bit_mask) && (val != ~fui_storage_t{0});
          return val + round_up;
        }
      }
    }
  }

  template<bitweight_t RoundingPos>
  constexpr auto round_to() const {
    static_assert(RoundingPos < format.msb_weight);
    if constexpr (RoundingPos <= format.lsb_weight) {
      return *this;
    } else {
      auto round = extract<RoundingPos - 1, RoundingPos - 1>();
      bool round_up = round.value_ != 0;
      auto up = extract<format.msb_weight, RoundingPos>();
      using res_t = decltype(up);
      return res_t{up.value_ + (round_up && up != res_t::max_val())};
    }
  }

  template<bitweight_t ShiftVal>
  constexpr auto shift() const {
    using out_format = FixedFormat<Format::msb_weight + ShiftVal, Format::lsb_weight + ShiftVal, typename Format::sign_t>;
    return out_format{value_};
  }

  /**
   * @brief Get the FT value nearest to represented value.
   *        Ties are rounded toward +infinity (with saturation).          
   * 
   * @tparam FT type of the converted result
   * @return constexpr FT 
   */
  template <std::floating_point FT> constexpr FT get_as() const {
    using limits = std::numeric_limits<FT>;
    static_assert(limits::max_exponent > format.msb_weight, "Overflow to infinity are not handled yet");
    static_assert(limits::min_exponent < format.lsb_weight);
    static_assert(limits::radix == 2);
    constexpr bool extra_bits = std::max(static_cast<int>(format.width - format.is_signed) - static_cast<int>(limits::digits), 0);
    if constexpr (extra_bits == 0) {  
      auto unscaled = static_cast<FT>(value_);
      auto scaled = std::ldexp(unscaled, format.lsb_weight);
      return scaled;
    } else {
      auto round = extract<format.lsb_weight + extra_bits - 1, format.lsb_weight + extra_bits - 1>();
      auto up = extract<format.msb_weight, format.lsb_weight + extra_bits>();
      using extended_format = FixedFormat<format.msb_weight + 1, format.lsb_weight, typename Format::sign_t>;
      auto val = up.as_hint();
      auto rounded = val.add_with_carry({0}, {round.value_});
      auto unscaled = static_cast<FT>(rounded.unravel());
      auto scaled = std::ldexp(unscaled, format.lsb_weight);
      return scaled;
    }
  }

  /**
   * @brief Get the OT nearest to represented value.
   *        Ties are rounded toward +infinity.
   * 
   * @tparam IT type of the converted result
   * @return constexpr IT 
   */
  template <std::integral IT> constexpr IT get_as() const {
    using limits = std::numeric_limits<IT>;
    static_assert(limits::radix == 2);

    if constexpr (format.msb_weight < 0) {
      // No integer part
      if constexpr (format.msb_weight < -1 || format.is_signed) {
        // If format is signed : values are in [-.5; .5) so ties toward zero
        // always return 0
        return 0;
      }
      constexpr auto val_mask = storage_t{1} << (format.width) - 1;
      return (value_ > val_mask) ? 1 : 0;
    } else {

      if constexpr (!limits::is_signed && format.is_signed) {
        if (value_ < 0)
          return 0;
      }
      // Max positive value: 1 << pos_bits - 1
      constexpr auto pos_bits_msb = limits::digits - 1;
      constexpr auto output_format_msb =
          pos_bits_msb + (limits::is_signed ? 1 : 0);
      constexpr auto format_pos_msb = format.msb_weight - format.is_signed;
      constexpr bool can_overflow = format_pos_msb > pos_bits_msb;
      if constexpr (can_overflow) {
        auto top_bits = extract<format.msb_weight, pos_bits_msb + 1>().value();
        if (top_bits < 0) {
          return limits::min();
        } else if (top_bits > 0) {
          return limits::max();
        }
      }
      constexpr bool needs_rounding_management = format.lsb_weight < 0;
      constexpr auto slice_msb = std::min(output_format_msb, format.msb_weight);
      constexpr auto slice_lsb = std::max(0, format.lsb_weight);
      auto int_slice = extract<slice_msb, slice_lsb>();
      using slice_sign_t =
          std::conditional_t<slice_msb == format.msb_weight,
                             typename Format::sign_t, unsigned>;
      if constexpr (needs_rounding_management) {
        // No need of extension
        auto frac_slice = extract<-1, format.lsb_weight>();
        constexpr auto cut_val = typename decltype(frac_slice)::storage_t{1}
                                 << -(frac_slice.format.lsb_weight + 1);
        bool round_up =
            frac_slice.value() > cut_val && int_slice != int_slice.max_val();
        return static_cast<IT>(int_slice.value_ + (round_up ? 1 : 0));
      } else {
        using out_format = FixedFormat<slice_msb, 0, slice_sign_t>;
        return static_cast<IT>(int_slice.extend_to(out_format{}).value());
      }
    }
  }

  static constexpr FixedNumber min_val() {
    return (Format::is_signed) ? storage_t{1} << (Format::width - 1)
                               : storage_t{0};
  }
  static constexpr FixedNumber max_val() {
    if constexpr (Format::is_signed) {
      return -(min_val().value_ + storage_t{1});
    } else {
      return ~storage_t{0};
    }
  }
  static constexpr FixedNumber min_pos() { return 1; }
  static constexpr FixedNumber min_neg() {
    if constexpr (format.is_signed) {
      return ~storage_t{0};
    } else {
      return 0;
    }
  }

  explicit constexpr FixedNumber(
      hint::BitIntWrapper<Format::width, Format::is_signed> const &val)
      : value_{val.unravel()} {}
  constexpr storage_t value() const { return value_; }
  constexpr hint::BitIntWrapper<Format::width, Format::is_signed>
  as_hint() const {
    return {value_};
  }

  template <bitweight_t high, bitweight_t low> constexpr auto extract() const {
    static_assert(high <= Format::msb_weight);
    static_assert(low >= Format::lsb_weight);
    constexpr bool signed_ret = Format::is_signed && high == Format::msb_weight;
    using ret_dim = FixedFormat<high, low, is_signed_to_type_t<signed_ret>>;
    auto val_hint = this->as_hint();
    auto ret_val = val_hint.template slice<high - Format::lsb_weight,
                                           low - Format::lsb_weight>();
    if constexpr (signed_ret) {
      return FixedNumber<ret_dim>(ret_val.as_signed().unravel());
    } else {
      return FixedNumber<ret_dim>(ret_val.unravel());
    }
  }

  static constexpr format_t format{};

  template <FixedFormatType NewFormat>
  constexpr FixedNumber<NewFormat> extend_to(NewFormat new_format = {}) const {
    static_assert(detail::can_extend_to(format, new_format),
                  "Trying to extend value to a format that is not a superset "
                  "of current format");
    using interm_format = hint::detail::bitint_base_t<format.is_signed, NewFormat::width>;
    auto interm_val = interm_format{value_};
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
  auto val1 = res_format.get_bit_int(op1.value());
  auto val2 = res_format.get_bit_int(op2.value());
  return FixedNumber(res_format, val1 * val2);
}

template <FixedNumberType FT>
constexpr bool operator==(FT const &op1, FT const &op2) {
  return op1.value() == op2.value();
}

template <FixedNumberType FT1, FixedNumberType FT2>
constexpr auto operator<=>(FT1 const & op1, FT2 const & op2) {
  return op1.value() <=> op2.value();
}

} // namespace archgenlib

#endif

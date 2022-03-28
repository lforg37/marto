#ifndef FIXEDPOINT_LITERAL_HPP
#define FIXEDPOINT_LITERAL_HPP

#include "fixedpoint.hpp"
#include "expression_types.hpp"
#include "bitint_tools/bitint_constant.hpp"

namespace archgenlib {

namespace detail {

template <typename... Ts> struct pack {};

template <bool fail, typename Types, auto... vals> struct print_assert {
  static_assert(fail, "print");
};

template <auto... vals>
struct print : print_assert<false, pack<decltype(vals)...>, vals...> {};

template <typename T, T base> struct strconv {
  static constexpr T get_value(char c) {
    if (c >= '0' && c <= '9')
      return c - '0';
    if (base == 16 && c >= 'a' && c <= 'f')
      return 10 + (c - 'a');
    if (base == 16 && c >= 'A' && c <= 'F')
      return 10 + (c - 'A');
    return base;
  }

  static constexpr T stoi_impl(std::string_view str, int idx, T value) {
    if (idx == str.size())
      return value;
    if (str[idx] == '-' && idx == 0)
      return -stoi_impl(str, idx + 1, 0);
    if (get_value(str[idx]) == base)
      throw "failed to parse digit";
    return stoi_impl(str, idx + 1, get_value(str[idx]) + value * base);
  }

  static constexpr T stoi(std::string_view str) { return stoi_impl(str, 0, 0); }
};

template <char... Bits> struct UDLHexFloatHelper {
  static constexpr unsigned count_hex_char(std::string_view str) {
    unsigned c = 0;
    while (c < str.size() && strconv<int, 16>::get_value(str[c]) != 16)
      c++;
    return c;
  }

  static constexpr char tab[]{Bits..., '\0'};
  static constexpr std::string_view signed_prefixed_num{tab};
  static constexpr bool is_signed = signed_prefixed_num[0] == '-';
  static constexpr std::string_view prefixed_num =
      signed_prefixed_num.substr(is_signed || signed_prefixed_num[0] == '+');
  static_assert(prefixed_num.starts_with("0x") ||
                    prefixed_num.starts_with("0X"),
                "Only hexadecimal literals are accepted");
  static constexpr std::string_view num = prefixed_num.substr(2);

  static constexpr std::string_view raw_num = num.substr(is_signed);
  static constexpr unsigned num_size = count_hex_char(raw_num);
  static constexpr std::string_view raw_num_only =
      num.substr(is_signed, num_size);
  template<auto size>
  static constexpr auto bit_int_size = size * 4 + 5; //to make the base always fit inside
  using tmp_int_storage = hint::detail::bitint_base_t<false, bit_int_size<num_size>>;
  static constexpr auto raw_num_val = strconv<tmp_int_storage, 16>::stoi(raw_num_only);
  static constexpr int highest_bit_set = bit_int_size<num_size> - hint::detail::clz<raw_num_val>();

  static constexpr std::string_view after_int_part = num.substr(num_size);
  static constexpr bool has_dot = after_int_part.starts_with(".");
  static constexpr std::string_view after_dot = after_int_part.substr(has_dot);
  static constexpr unsigned decimal_part_size = count_hex_char(after_dot);
  static constexpr std::string_view decimal_part = after_dot.substr(0, decimal_part_size);
  using tmp_dec_storage = hint::detail::bitint_base_t<false, bit_int_size<decimal_part_size>>;
  static constexpr bool has_no_dec = (!has_dot || decimal_part_size == 0);
  static constexpr auto decimal_value =
      has_no_dec ? static_cast<tmp_dec_storage>(0)
                 : strconv<tmp_dec_storage, 16>::stoi(decimal_part);
  static constexpr int lowest_bit_set = std::min<int>(
      0, hint::detail::ctz<decimal_value>() - (decimal_part_size * 4));
  static constexpr auto corrected_decimal_value =
      decimal_value >>
      std::min<int>((decimal_part_size * 4), hint::detail::ctz<decimal_value>());
  // print<decimal_part_size, lowest_bit_set, hint::detail::ctz<decimal_value>(),
  //       decimal_part_size * 4>
  //     p0;

  static constexpr std::string_view postfix = after_dot.substr(decimal_part_size);
  static constexpr bool has_p = postfix.starts_with("p");
  static constexpr std::string_view low_prec_str = postfix.substr(has_p);
  static constexpr auto prec_shift = has_p ? strconv<int, 10>::stoi(low_prec_str) : 0;

  // print<has_p, prec_shift> p;
  static constexpr auto high_prec = highest_bit_set + prec_shift + is_signed;
  static constexpr auto low_prec = lowest_bit_set + prec_shift;
  static constexpr auto bit_count = highest_bit_set - lowest_bit_set + is_signed;
  using final_storage = hint::detail::bitint_base_t<is_signed, std::max(bit_count, 1)>;
  // print<is_signed, high_prec, low_prec, bit_count> p2;
  static constexpr final_storage final_value = {
      static_cast<final_storage>(raw_num_val)
          << static_cast<final_storage>(-lowest_bit_set) |
      static_cast<final_storage>(corrected_decimal_value)};
  // print<high_prec, low_prec, raw_num_val, corrected_decimal_value, final_value> p2;
  static_assert((high_prec - low_prec) == bit_count, "");

  using dim = archgenlib::FixedFormat<high_prec, low_prec,
                        archgenlib::is_signed_to_type_t<is_signed>>;

  static constexpr final_storage value{final_value};
};
} // namespace detail

template <char... Bits> auto operator"" _cst() {
  using parser = detail::UDLHexFloatHelper<Bits...>;
  return archgenlib::Constant<
      archgenlib::FixedConstant<typename parser::dim, parser::value>>{};
}

template <char... Bits> auto operator"" _fixed() {
  using parser = detail::UDLHexFloatHelper<Bits...>;
  return FixedNumber<typename parser::dim>{
      detail::UDLHexFloatHelper<Bits...>::value};
}

}

using archgenlib::operator"" _cst;
using archgenlib::operator"" _fixed;

#endif

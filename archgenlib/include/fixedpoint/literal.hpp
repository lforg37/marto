#ifndef FIXEDPOINT_LITERAL_HPP
#define FIXEDPOINT_LITERAL_HPP

#include "fixedpoint.hpp"
#include "expression_types.hpp"

namespace archgenlib {

namespace detail {

template <typename T, int base> struct strconv {
  static constexpr unsigned get_value(char c) {
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
    while (strconv<int, 16>::get_value(str[c]) != 16)
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
  static constexpr auto raw_num = num.substr(is_signed);
  static constexpr auto num_size = count_hex_char(raw_num);
  static constexpr std::string_view raw_num_only =
      num.substr(is_signed, num_size);
  static constexpr auto raw_num_val = strconv<uint64_t, 16>::stoi(raw_num_only);
  static constexpr std::string_view postfix = num.substr(num_size);
  static_assert(postfix.starts_with(".p"), "invalid format");
  static constexpr std::string_view low_prec_str = postfix.substr(2);
  static constexpr auto low_prec = strconv<int, 10>::stoi(low_prec_str);
  static constexpr auto highest_bit_set = 64 - __builtin_clzll(raw_num_val) - 1;
  static constexpr auto high_prec = highest_bit_set + low_prec + is_signed;
  static constexpr auto bit_count = high_prec - low_prec + 1;
  // print_many<is_signed, raw_num_val, low_prec, high_prec> p;

  using dim = archgenlib::FixedFormat<high_prec, low_prec,
                        archgenlib::is_signed_to_type_t<is_signed>>;

  static constexpr hint::detail::bitint_base_t<is_signed, bit_count> value{raw_num_val};
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

#include <cmath>
#include <iostream>
#include <math.h>
#include <type_traits>
#include <cassert>

#include "bitint_tools/type_helpers.hpp"
#include "hint.hpp"

#include "fixedpoint/literal.hpp"
#include "fixedpoint/operators.hpp"
#include "fixedpoint/evaluate.hpp"
#include "fixedpoint/fixedpoint.hpp"

template <typename T> static constexpr bool ok = false;

template <archgenlib::ExpressionType T> static constexpr bool ok<T> = true;

using archgenlib::bitweight_t;

auto compute_ref(auto val) {
  auto val_d = val.template get_as<double>();
  auto sin = std::sin(val_d * M_PI);
  return sin;
}

bool check_against_ref(auto val, auto res) {
  static int count = 0;
  count++;

  auto x = val.template get_as<double>();
  auto outprec = decltype(res)::format_t::lsb_weight;
  auto resval = res.template get_as<double>();
  auto ref = compute_ref(val);
  auto diffabs = std::abs(ref - resval);
  static const double err_budget = ldexp(double{1}, outprec);
  if (diffabs < err_budget)
    return true;
  static int errors = 0;
  errors++;
  std::cerr << std::hex << "f(0x" << (unsigned)val.get_representation() << ")=0x" << (unsigned)res.get_representation()<< " ";
  std::cerr << std::dec << "error(" << errors << "/" << count << ") for f(" << x << ") (" << ref << " - " << resval << ") = " << diffabs << " > " << err_budget << std::endl;
  return false; 
} 

using fpdim_t = archgenlib::FixedFormat<-1, -3, unsigned>;
using fpnum_t = archgenlib::FixedNumber<fpdim_t>;

using storage_t =
    hint::detail::bitint_base_t<fpdim_t::is_signed, fpdim_t::width>;

constexpr bitweight_t outprec = -10;

inline auto convert_to_double(storage_t val) {
  auto val_d = static_cast<double>(val);
  val_d = ldexp(val_d, fpdim_t::lsb_weight);
  return val_d;
}

auto compute_ref(storage_t val) {
  auto val_d = convert_to_double(val);
  auto sin = std::sin(val_d / 64 * std::numbers::pi);
  return sin;
}

bool compare_ref(storage_t val, auto res) {
  auto res_int = res.value();
  if (res_int > (1 << (res.width - 1))) {
    res_int -= (1 << res.width);
  }
  auto resval = static_cast<double>(res_int);
  resval = ldexp(resval, outprec);
  auto ref = compute_ref(val);
  auto diffabs = std::abs(ref - resval);
  static const double err_budget = ldexp(double{1}, outprec);
  return diffabs < err_budget;
}

int main() {
  using storage_t = unsigned _BitInt(fpdim_t::width);
  for (unsigned int i = 0; i < (1 << fpdim_t::width); ++i) {
    auto val = static_cast<storage_t>(i);
    fpnum_t val_fixed{val};
    auto a = archgenlib::FreeVariable(val_fixed);
    auto s = a * 0x1p-7_cst * archgenlib::pi;
    auto c = archgenlib::sin(s);
    auto res = archgenlib::evaluate<archgenlib::FixedFormat<1, -10, signed>>(c);
    if constexpr (archgenlib::has_specialization_header) {
      assert(compare_ref(val, res));
    }
  }
  if constexpr (archgenlib::has_implementation) { 
    std::cout << "All results seems correct !\n";
  }
  return 0;
}

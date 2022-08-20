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

/// debug function provided by archgen MLIR for now
extern "C" double polynom_double_0(double);

template <typename T> static constexpr bool ok = false;

template <archgenlib::ExpressionType T> static constexpr bool ok<T> = true;

using archgenlib::bitweight_t;

auto compute_ref(auto val) {
  auto val_d = val.template get_as<double>();
  auto sin = std::sin(val_d /2 * M_PI);
  return sin;
}

bool check_against_ref(auto val, auto res) {
  static int count = 0;
  count++;

  auto x = val.template get_as<double>();
  auto outprec = decltype(res)::format_t::lsb_weight;
  // auto resval = res.template get_as<double>();
  auto ref = compute_ref(val);
  auto resval = polynom_double_0(x);
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
 
int main() {
  using storage_t = unsigned _BitInt(fpdim_t::width);
  for (unsigned int i = 0; i < (1 << fpdim_t::width); ++i) {
    auto val = static_cast<storage_t>(i);
    fpnum_t val_fixed{val};
    auto a = archgenlib::FreeVariable(val_fixed);
    auto c = archgenlib::sin(a * archgenlib::pi / 0x2p0_cst);
    auto res = archgenlib::evaluate<archgenlib::FixedFormat<-1, -3, unsigned>>(c);
    if constexpr (archgenlib::has_implementation) {
      check_against_ref(val_fixed, res); 
    }
  }
  if constexpr (archgenlib::has_implementation) { 
    std::cout << "All results seems correct !\n";
  }
  return 0;
}

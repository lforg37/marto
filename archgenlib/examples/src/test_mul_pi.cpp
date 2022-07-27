#include <cstdint>
#include <string_view>
#include <type_traits>

#include "bitint_tools/type_helpers.hpp"
#include "hint.hpp"

#include "fixedpoint/constants.hpp"
#include "fixedpoint/expression.hpp"
#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"
#include "fixedpoint/literal.hpp"

__VITIS_KERNEL auto
mul_pi(archgenlib::FixedNumber<archgenlib::FixedFormat<2, -1, unsigned>> val) {
  archgenlib::Variable x{val};
  auto expr = archgenlib::pi * x;
  return archgenlib::evaluate<-2>(expr);
}

#include <iostream>
#include <type_traits>

#include "fixedpoint/expression.hpp"
#include "fixedpoint/expression_types.hpp"

#ifdef INCLUDE_GENERATED_HEADER
#include "specialization_header.hpp"
#endif

template <typename T> static constexpr bool ok = false;

template <archgenlib::ExpressionType T> static constexpr bool ok<T> = true;

int main() {
  archgenlib::Variable<int> a{17};
  archgenlib::Constant<std::integral_constant<int, 42>> b{};
  auto c = a + b;
  auto res = archgenlib::evaluate<-14>(c);
  std::cout << res << std::endl;
  return 0;
}
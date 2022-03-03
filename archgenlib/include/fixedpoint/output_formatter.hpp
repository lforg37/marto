#ifndef FIXEDPOINT_OUTPUT_FORMATTER_HPP
#define FIXEDPOINT_OUTPUT_FORMATTER_HPP

#include <array>
#include <fstream>
#include <string>
#include <string_view>
#include <utility>

namespace archgenlib {
namespace detail {

// Utility to get the typename as a string,
// Taken almost as is from
// https://bitwizeshift.github.io/posts/2021/03/09/getting-an-unmangled-type-name-at-compile-time/

template <typename T> constexpr auto type_name_str() {
#if defined(__clang__)
  constexpr auto prefix = std::string_view{"[T = "};
  constexpr auto suffix = std::string_view{"]"};
  constexpr auto function = std::string_view{__PRETTY_FUNCTION__};
#elif defined(__GNUC__)
  constexpr auto prefix = std::string_view{"with T = "};
  constexpr auto suffix = std::string_view{"]"};
  constexpr auto function = std::string_view{__PRETTY_FUNCTION__};
#elif defined(_MSC_VER)
  constexpr auto prefix = std::string_view{"type_name_str<"};
  constexpr auto suffix = std::string_view{">(void)"};
  constexpr auto function = std::string_view{__FUNCSIG__};
#else
#error Unsupported compiler
#endif

  constexpr auto start = function.find(prefix) + prefix.size();
  constexpr auto end = function.rfind(suffix);

  static_assert(start < end);

  constexpr auto name = function.substr(start, (end - start));
  return name;
}

template <typename T> struct type_name_holder {
  static inline constexpr auto value = type_name_str<T>();
};

template <typename T> constexpr auto type_name() {
  constexpr auto &value = type_name_holder<T>::value;
  return std::string_view{value};
}

template <typename T> constexpr auto type_name_from_val(T) {
  return type_name<T>();
}

/// End of utility

std::string doReturn(std::string_view expr);

class SpecialisationFormatter {
public:
  std::ofstream output;
  SpecialisationFormatter();
  ~SpecialisationFormatter();
};
SpecialisationFormatter &getFormatter();
} // namespace detail

} // namespace archgenlib

#endif

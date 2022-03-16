#include "runtime/expression_tree.hpp"
#include "runtime/output_formatter.hpp"
#include <ostream>
#include <string_view>

namespace archgenlib {
namespace detail {
std::ostream &print_path(std::ostream &os,
                         ExpressionRTRepr::path_t const &path) {
  for (auto iter = path.rbegin(); iter != path.rend(); iter++) {
    switch (*iter) {
    case ExpressionRTRepr::DirectionInTree::LEFT:
      os << ".left";
      break;
    case ExpressionRTRepr::DirectionInTree::RIGHT:
      os << ".right";
      break;
    case ExpressionRTRepr::DirectionInTree::DIRECT:
      os << ".op";
      break;
    }
  }
  os << ".value";
  return os;
}
} // namespace detail

}; // namespace archgenlib

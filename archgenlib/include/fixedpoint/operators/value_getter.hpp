#ifndef FIXEDPOINT_OPERATORS_VALUE_GETTER_HPP
#define FIXEDPOINT_OPERATORS_VALUE_GETTER_HPP

#include "runtime/expression_tree.hpp"
#include "runtime/output_formatter.hpp"
#include <sstream>
#include <string_view>

namespace archgenlib {
struct ExprLeafGetter {
  static auto getOperator(ExpressionRTRepr::path_t const & path) {
    return [path](std::string_view expr) {
        std::stringstream ss;
        ss << expr;
        for (auto iter = path.rbegin() ; iter != path.rend() ; iter++) 
        {
            switch(*iter) {
                case ExpressionRTRepr::DirectionInTree::LEFT:
                ss << ".left";
                break;
                case ExpressionRTRepr::DirectionInTree::RIGHT:
                ss << ".right";
                break;
                case ExpressionRTRepr::DirectionInTree::DIRECT:
                ss << ".op";
                break;
            }
        }
        ss << ".value";
        return ss.str();
    };
  } 
};
} // namespace archgenlib

#endif

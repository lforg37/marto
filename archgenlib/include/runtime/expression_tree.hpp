#ifndef FIXEDPOINT_EXPRESSION_RT_TREE
#define FIXEDPOINT_EXPRESSION_RT_TREE

#include <cassert>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "fixedpoint/expression_types.hpp"
#include "fixedpoint/fixedpoint.hpp"
#include "fixedpoint/operations.hpp"
#include "output_formatter.hpp"

namespace archgenlib {
struct FixedFormatRTRepr {
  bitweight_t const msb_weight;
  bitweight_t const lsb_weight;
  vecwidth_t width;
  bool const is_signed;
  FixedFormatRTRepr(bitweight_t msb, bitweight_t lsb, bool signedness)
      : msb_weight{msb}, lsb_weight{lsb}, is_signed{signedness} {
    assert((msb >= lsb) && "MSB should have a greater weight than LSB");
    width = static_cast<vecwidth_t>(msb - lsb + 1);
  }
  std::string toFixedFormatName() const;
  std::string toFPNumName() const;
};

template <ExpressionType ET> struct ExprTypeHolder {};

struct BinaryOpRTNodeContent;
struct UnaryOpRTNodeContent;
struct NullaryOpLeafRTNodeContent;
struct ConstantLeafRTNodeContent;
struct VariableLeafRTNodeContent;

template <typename ContentType> class OpNodePrototype {
  using contents_t = std::unique_ptr<ContentType>;
  contents_t contents;
  using ul_t = ContentType;

public:
  std::uint32_t const id;
  ContentType const &operator*() const { return *contents.get(); }
  ContentType const *operator->() const { return contents.get(); }
  OpNodePrototype(contents_t &&c, std::uint32_t id)
      : contents{std::move(c)}, id{id} {}
};

using BinaryOpRTNode = OpNodePrototype<BinaryOpRTNodeContent>;
using UnaryOpRTNode = OpNodePrototype<UnaryOpRTNodeContent>;
using NullaryOpLeafRTNode = OpNodePrototype<NullaryOpLeafRTNodeContent>;
using ConstantLeafRTNode = OpNodePrototype<ConstantLeafRTNodeContent>;
using VariableLeafRTNode = OpNodePrototype<VariableLeafRTNodeContent>;

using RTNode = std::variant<BinaryOpRTNode, UnaryOpRTNode, NullaryOpLeafRTNode,
                            ConstantLeafRTNode, VariableLeafRTNode>;

struct BinaryOpRTNodeContent {
  RTNode left;
  RTNode right;
  OperationKind const operation_kind;
  BinaryOpRTNodeContent(RTNode &left, RTNode &right, OperationKind op)
      : left{std::move(left)}, right{std::move(right)}, operation_kind{op} {}
};

struct UnaryOpRTNodeContent {
  RTNode operand;
  OperationKind const operation_kind;
  UnaryOpRTNodeContent(RTNode &operand, OperationKind op_type)
      : operand{std::move(operand)}, operation_kind{op_type} {}
};

struct NullaryOpLeafRTNodeContent {
  OperationKind const operation_kind;
  NullaryOpLeafRTNodeContent(OperationKind op_type) : operation_kind{op_type} {}
};

struct ConstantLeafRTNodeContent {
  FixedFormatRTRepr const dim;
  std::string const constant_representation;
  ConstantLeafRTNodeContent(FixedFormatRTRepr const &fpdim,
                            std::string_view repres)
      : dim{fpdim}, constant_representation{repres} {}
};

struct VariableLeafRTNodeContent {
  FixedFormatRTRepr const dim;
  VariableLeafRTNodeContent(FixedFormatRTRepr const &fpdim) : dim{fpdim} {}
};

struct ExpressionRTRepr {
  enum struct DirectionInTree { LEFT, RIGHT, DIRECT };

  // Stored reversed
  using path_t = std::vector<DirectionInTree>;

  struct SymbolRecord {
    VariableLeafRTNodeContent const &description;
    path_t path_from_root;
  };

  using sym_table_t = std::map<std::uint32_t, SymbolRecord>;

private:
  template <FixedFormatType T> auto descriptor_from_fpdim(T const &dim) {
    return FixedFormatRTRepr{dim.msb_weight, dim.lsb_weight, dim.is_signed};
  }

  ///
  ///@brief Extract the constant value string representation from the constant
  /// typename
  ///
  ///@param constant_typename
  ///@return std::string
  static std::string extract_constant_repr(std::string_view constant_typename);

  template <BinaryExprType ET> auto extract_subtree() {
    auto [left_child_op_node, left_child_sym_table] =
        extract_subtree<typename ET::LeftType>();
    auto [right_child_op_node, right_child_sym_table] =
        extract_subtree<typename ET::RightType>();
    for (auto &[addr, symbol] : left_child_sym_table) {
      symbol.path_from_root.emplace_back(DirectionInTree::LEFT);
    }
    for (auto &[addr, symbol] : right_child_sym_table) {
      symbol.path_from_root.emplace_back(DirectionInTree::RIGHT);
    }
    left_child_sym_table.merge(right_child_sym_table);
    RTNode op_node{BinaryOpRTNode{std::make_unique<BinaryOpRTNodeContent>(
                                      left_child_op_node, right_child_op_node,
                                      ET::operation_t::operation_kind),
                                  first_available_node_id++}};
    return std::make_tuple(std::move(op_node), left_child_sym_table);
  }

  template <UnaryExprType ET> auto extract_subtree() {
    auto [child_op_node, child_sym_table] =
        extract_subtree<typename ET::ChildType>();
    for (auto &[addr, symbol] : child_sym_table) {
      symbol.path_from_root.emplace_back(DirectionInTree::DIRECT);
    }
    RTNode op_node{
        UnaryOpRTNode{std::make_unique<UnaryOpRTNodeContent>(
                          child_op_node, ET::operation_t::operation_kind),
                      first_available_node_id++}};
    return std::make_tuple(std::move(op_node), child_sym_table);
  }

  template <NullaryOpExprType ET> auto extract_subtree() {
    RTNode res{NullaryOpLeafRTNode{std::make_unique<NullaryOpLeafRTNodeContent>(
                                       ET::operation_t::operation_kind),
                                   first_available_node_id++}};
    return std::make_tuple(std::move(res), sym_table_t{});
  }

  template <VariableExprType ET> auto extract_subtree() {
    auto id = first_available_node_id++;
    auto fpdim = descriptor_from_fpdim(typename ET::dimension_t{});
    RTNode res{VariableLeafRTNode{
        std::make_unique<VariableLeafRTNodeContent>(fpdim), id}};

    sym_table_t sym_table{};
    sym_table.insert(std::make_pair(
        id, SymbolRecord{*std::get<VariableLeafRTNode>(res), {}}));
    return std::make_tuple(std::move(res), sym_table);
  }

  template <ConstantExprType ET> auto extract_subtree() {
    auto dim = descriptor_from_fpdim(typename ET::dimension_t{});
    auto descr = extract_constant_repr(detail::type_name<typename ET::type>());
    RTNode res{ConstantLeafRTNode{
        std::make_unique<ConstantLeafRTNodeContent>(dim, descr),
        first_available_node_id++}};
    return std::make_tuple(std::move(res), sym_table_t{});
  }

public:
  template <ExpressionType ET>
  ExpressionRTRepr(ExprTypeHolder<ET>) : first_available_node_id{0} {
    auto [root_node, sym_table] = extract_subtree<ET>();
    symbol_table = std::move(sym_table);
    root.emplace(std::move(root_node));
  }

  /***
   * Return a map that associates to each node the number of dependant variables
   */
  std::map<std::uint32_t, std::size_t> get_count_for_nodes() const;

  std::vector<RTNode const *> get_singlevar_dominants() const;

  sym_table_t symbol_table;
  std::optional<RTNode> root;
  std::uint32_t first_available_node_id;
};
} // namespace archgenlib

#endif

#ifndef FIXEDPOINT_EXPRESSION_RT_TREE
#define FIXEDPOINT_EXPRESSION_RT_TREE

#include <atomic>
#include <cstdint>
#include <list>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "expression_types.hpp"
#include "fixedpoint.hpp"
#include "operations.hpp"

namespace archgenlib {
struct FPDimRTRepr {
  bitweight_t msb_weight;
  bitweight_t lsb_weight;
  bool is_signed;
};

struct ExpressionRTRepr {
  enum struct DirectionInTree { LEFT, RIGHT, DIRECT };

  enum struct LeafType { VARIABLE, CONSTANT };

  // Stored reversed
  using path_t = std::list<DirectionInTree>;

  struct LeafDescriptor {
    LeafType leaf_type;
    FPDimRTRepr dim;
  };

  struct LeafId {
    LeafDescriptor description;
    path_t path_from_root;
  };

  struct OperationNode {
    std::optional<OperationKind> op_kind;
    std::optional<LeafDescriptor> leaf_descr;
    std::list<OperationNode> children;
    bool is_constant;
    void DFS(auto &lambda) const {
      for (auto const &child : children)
        child.DFS(lambda);
      lambda(*this);
    }

    /*Predicate should return whether the exploration should end on this node */
    void BFS(auto &lambda) const {
      if (!lambda(*this)) {
        for (auto &child : children) {
          child.BFS(lambda);
        }
      }
    }
  };

  using sym_table_t = std::map<std::uint32_t, LeafId>;

private:
  template <FPDimType T> auto descriptor_from_fpdim(T const &dim) {
    return FPDimRTRepr{dim.msb_weight, dim.lsb_weight, dim.is_signed};
  }

  template <BinaryExprType ET> auto extract_subtree() {
    auto [left_child_op_node, left_child_sym_table] =
        extract_subtree<typename ET::LeftType>();
    auto [right_child_op_node, right_child_sym_table] =
        extract_subtree<typename ET::RightType>();
    OperationNode op_node{};
    op_node.op_kind.emplace(ET::operation_t::operation_kind);
    op_node.children.emplace_back(std::move(left_child_op_node));
    op_node.children.emplace_back(std::move(right_child_op_node));
    for (auto &[addr, symbol] : left_child_sym_table) {
      symbol.path_from_root.emplace_back(DirectionInTree::LEFT);
    }
    for (auto &[addr, symbol] : right_child_sym_table) {
      symbol.path_from_root.emplace_back(DirectionInTree::RIGHT);
    }
    left_child_sym_table.merge(right_child_sym_table);
    op_node.is_constant =
        left_child_op_node.is_constant && right_child_op_node.is_constant;
    return std::make_tuple(op_node, left_child_sym_table);
  }

  template <UnaryExprType ET> auto extract_subtree() {
    auto [child_op_node, child_sym_table] =
        extract_subtree<typename ET::ChildType>();
    OperationNode op_node{};
    op_node.op_kind.emplace(ET::operation_t::operation_kind);
    op_node.children.emplace_back(std::move(child_op_node));
    for (auto &[addr, symbol] : child_sym_table) {
      symbol.path_from_root.emplace_back(DirectionInTree::DIRECT);
    }
    op_node.is_constant = child_op_node.is_constant;
    return std::make_tuple(op_node, child_sym_table);
  }

  template <VariableExprType ET> auto extract_subtree() {
    OperationNode op_node{};
    LeafDescriptor ldescr{LeafType::VARIABLE,
                          descriptor_from_fpdim(typename ET::dimension_t{})};

    sym_table_t sym_table{};
    sym_table.insert(
        std::make_pair(first_available_leaf_id++, LeafId{ldescr, {}}));
    op_node.leaf_descr = ldescr;
    op_node.is_constant = false;
    return std::make_tuple(op_node, sym_table);
  }

  template <ConstantExprType ET> auto extract_subtree() {
    OperationNode op_node{};
    LeafDescriptor leaf_descr{
        LeafType::CONSTANT, descriptor_from_fpdim(typename ET::dimension_t{})};

    sym_table_t sym_table{};
    sym_table.insert(
        std::make_pair(first_available_leaf_id++, LeafId{leaf_descr, {}}));
    op_node.leaf_descr = leaf_descr;
    op_node.is_constant = true;
    return std::make_tuple(op_node, sym_table);
  }

public:
  ExpressionRTRepr() : first_available_leaf_id{0} {}
  template <ExpressionType ET> void do_init() {
    auto [root_node, sym_table] = extract_subtree<ET>();
    symbol_table.emplace(std::move(sym_table));
    root.emplace(std::move(root_node));
  }

  /***
   * Return a map that associates to each node the number of dependant variables
   */
  std::map<OperationNode const *, std::size_t>
  get_count_for_nodes() const;

  std::list<OperationNode const *> get_singlevar_dominants() const;

  std::uint32_t first_available_leaf_id;
  std::optional<sym_table_t> symbol_table;
  std::optional<OperationNode> root;
};
} // namespace archgenlib

#endif

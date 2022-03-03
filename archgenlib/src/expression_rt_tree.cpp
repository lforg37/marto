#include "fixedpoint/expression_rt_tree.hpp"
#include <utility>

namespace archgenlib {

std::map<ExpressionRTRepr::OperationNode const *, std::size_t> ExpressionRTRepr::get_count_for_nodes() const {
  std::map<ExpressionRTRepr::OperationNode const *, std::size_t> var_count{};
  auto check_node = [&var_count](OperationNode const & op_node) {
    size_t current_count;
    if (op_node.children.empty()) {
      current_count = (op_node.is_constant) ? 0 : 1;
    } else {
      current_count = 0;
      for (auto& child : op_node.children) {
        current_count += var_count[&child];
      }
    }
    var_count.insert(std::make_pair(&op_node, current_count));
  };
  root->DFS(check_node);
  return var_count;
}

std::list<ExpressionRTRepr::OperationNode const *> ExpressionRTRepr::get_singlevar_dominants() const{
  std::list<ExpressionRTRepr::OperationNode const *> ret_list{};
  auto var_count = get_count_for_nodes();
  auto pred = [&ret_list, &var_count](ExpressionRTRepr::OperationNode const & op_node) {
      if (var_count[&op_node] <= 1) {
        // Node is a maximum
        ret_list.push_back(&op_node);
        return true;
      }
      return false;
  };
  root->BFS(pred);
  return ret_list;
}

}

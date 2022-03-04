#include <cstddef>
#include <cstdint>
#include <utility>
#include <variant>

#include "fixedpoint/fixedpoint.hpp"
#include "runtime/expression_tree.hpp"
#include "runtime/output_formatter.hpp"

namespace {
std::uint32_t id(archgenlib::RTNode const &val) {
  return std::visit([]<typename T>(T const &v) { return v.id; }, val);
}

template <typename FuncType> struct HierarchicalVisitor {

  void operator()(auto const &val) {}

  void operator()(archgenlib::UnaryOpRTNode const &val) {
    (*this)(val->operand);
  }

  void operator()(archgenlib::BinaryOpRTNode const &val) {
    (*this)(val->left);
    (*this)(val->right);
  }

  void operator()(archgenlib::RTNode const &node) {
    if (!functor(node)) {
      std::visit(*this, node);
    }
  }

  FuncType &functor;

  HierarchicalVisitor(FuncType &functor) : functor{functor} {}
};

struct NodeFreeVarCounterVisitor {
private:
  std::map<std::uint32_t, std::size_t> account_holder;

public:
  void visit(archgenlib::ConstantLeafRTNode const &val) {
    account_holder.emplace(val.id, 0);
  }

  void visit(archgenlib::VariableLeafRTNode const &val) {
    account_holder.emplace(val.id, 1);
  }

  void visit(archgenlib::UnaryOpRTNode const &val) {
    visit(val->operand);
    auto childval = account_holder.at(id(val->operand));
    account_holder.insert(std::make_pair(val.id, childval));
  }

  void visit(archgenlib::BinaryOpRTNode const &val) {
    visit(val->left);
    visit(val->right);
    auto leftval = account_holder.at(id(val->left));
    auto rightval = account_holder.at(id(val->right));
    account_holder.insert(std::make_pair(val.id, leftval + rightval));
  }

  void visit(archgenlib::RTNode const &node) {
    std::visit([this]<typename T>(T const &n) { visit(n); }, node);
  }

public:
  auto operator()(auto const &val) {
    visit(val);
    return account_holder;
  }
};

} // namespace

namespace archgenlib {

std::map<std::uint32_t, std::size_t>
ExpressionRTRepr::get_count_for_nodes() const {
  return NodeFreeVarCounterVisitor{}(*root);
}

std::string
ExpressionRTRepr::extract_constant_repr(std::string_view const_tname) {
  using fpdim_t = FPDim<37, -12, false>;
  using indicator_t =
      FixedConstant<fpdim_t,
                    hint::detail::bitint_base_t<false, 37 + 12 + 1>{94}>;
  constexpr auto fpdim_name = detail::type_name<fpdim_t>();
  constexpr auto indicator_name = detail::type_name<indicator_t>();
  constexpr auto start_idx = indicator_name.find(fpdim_name);
  // End of the FPDim type representation
  auto fdim_endpos = const_tname.find('>', start_idx);
  auto number_start_pos = const_tname.find_first_of("0123456789-", fdim_endpos);
  auto number_end_pos = const_tname.find_first_not_of("0123456789xabcdefABCDEF",
                                                      number_start_pos + 1);
  return std::string{
      const_tname.substr(number_start_pos, number_end_pos - number_start_pos)};
}

std::vector<RTNode const *> ExpressionRTRepr::get_singlevar_dominants() const {
  std::vector<RTNode const *> ret_list{};
  auto var_count = get_count_for_nodes();
  auto pred = [&ret_list, &var_count](RTNode const &op_node) {
    if (var_count.at(id(op_node)) <= 1) {
      // Node is a maximum
      ret_list.emplace_back(&op_node);
      return true;
    }
    return false;
  };
  HierarchicalVisitor{pred}(*root);
  return ret_list;
}

} // namespace archgenlib

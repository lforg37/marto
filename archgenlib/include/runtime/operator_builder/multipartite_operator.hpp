#ifndef RUNTIME_OPERATOR_BUILDER_MULTIPARTITE_OPERATOR
#define RUNTIME_OPERATOR_BUILDER_MULTIPARTITE_OPERATOR

#include "runtime/fixfunction_multipartite.hpp"
#include "runtime/operator_builder/operator.hpp"
namespace archgenlib {

struct MultipartiteOperator {
  MultipartiteFunction & mpf;
  CPPOperator get_operator();
};

}

#endif
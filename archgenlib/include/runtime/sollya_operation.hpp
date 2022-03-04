#ifndef FIXEDPOINT_SOLLYA_OPERATION_HPP
#define FIXEDPOINT_SOLLYA_OPERATION_HPP

#include "sollya_handler.hpp"
#include "expression_tree.hpp"
#include <cstddef>
#include <sollya.h>

namespace archgenlib {

///
///@brief Build a sollya function object from a node
///
///@param op_node a node that should contain at most one variable 
///       leaf and any number of constants 
///@return a sollya_obj_t containing a constant expression or unary function equivalent 
///        to the operation tree described in the node. 
sollya_obj_t sollyaFunctionFromNode(RTNode const & op_node);
}
#endif

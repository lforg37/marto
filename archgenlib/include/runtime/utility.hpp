#ifndef RUNTIME_UTILITY_HPP
#define RUNTIME_UTILITY_HPP

#include <sollya.h>

#include "runtime/expression_tree.hpp"

namespace archgenlib {
sollya_obj_t sollya_interval_from_rtfpdim(FPDimRTRepr dim);
}

#endif

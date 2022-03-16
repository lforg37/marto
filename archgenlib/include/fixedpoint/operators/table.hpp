#ifndef FIXEDPOINT_OPERATORS_TABLE_HPP
#define FIXEDPOINT_OPERATORS_TABLE_HPP

#include <array>

#include "fixedpoint/fixedpoint.hpp"

namespace archgenlib {

template<vecwidth_t inputwidth, FPDimType OutputDimT>
using Table = std::array<FPNumber<OutputDimT>, 1 << inputwidth>;

}
#endif
#ifndef FIXEDPOINT_OPERATORS_TABLE_HPP
#define FIXEDPOINT_OPERATORS_TABLE_HPP

#include <array>

#include "fixedpoint/fixedpoint.hpp"

namespace archgenlib {

template<vecwidth_t inputwidth, FixedFormatType OutputDimT>
using Table = std::array<FixedNumber<OutputDimT>, 1 << inputwidth>;

}
#endif
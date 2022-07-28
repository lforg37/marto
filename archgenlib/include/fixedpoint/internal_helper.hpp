/// This is intended to be included from the specialization header only

#ifndef FIXEDPOINT_OPERATORS_TABLE_HPP
#define FIXEDPOINT_OPERATORS_TABLE_HPP

#include <array>

#include "fixedpoint/fixedpoint.hpp"
#include "bitint_tools/type_helpers.hpp"

using hint::operator""_sbi;
using hint::operator""_ubi;

namespace archgenlib {

template<vecwidth_t inputwidth, FixedFormatType OutputDimT>
using Table = std::array<FixedNumber<OutputDimT>, 1 << inputwidth>;

}
#endif

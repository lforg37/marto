#ifndef IEEE_ROUNDING_HPP
#define IEEE_ROUNDING_HPP

#include "ieeetype.hpp"

template<template<unsigned int, bool> class Wrapper>
Wrapper<1, false> ieee_getRoundBit(
		Wrapper<1, false> sign,
		Wrapper<1, false> lastFracBit,
		Wrapper<1, false> roundBit,
		Wrapper<1, false> sticky,
		IEEERoundingMode roundingMode
		)
{
	Wrapper<3, false> roundingCode{static_cast<uint8_t>(roundingMode)};
	auto b0 = roundingCode.template get<0>();
	auto b1 = roundingCode.template get<1>();
	auto b2 = roundingCode.template get<2>();

	auto any = b2 & b0 & b1 ^ sign.invert();

	auto always = b0 & (b1 ^ sign);

	auto tieBreak = Wrapper<1, false>::mux(b1, sign, lastFracBit);

	auto result = b2 & ((sticky & always) | (roundBit & (sticky | tieBreak | always)));
	return result;
}

#endif // IEEE_ROUNDING_HPP

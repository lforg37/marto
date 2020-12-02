#ifndef IEEE_ROUNDING_HPP
#define IEEE_ROUNDING_HPP

#include "ieeetype.hpp"

//#define DEBUG_IEEE_ROUNDING
#ifdef DEBUG_IEEE_ROUNDING
#include <iostream>
#include "tools/printing.hpp"

using std::cout;
using std::endl;

using hint::to_string;
#endif


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


	auto any = b2 & (b0.invert() | (b1 ^ sign));

	auto always = b0 & (b1 ^ sign);

	auto tieBreak = b1 | lastFracBit;

	auto result = any & ((sticky & always) | (roundBit & (sticky | tieBreak | always)));

#ifdef DEBUG_IEEE_ROUNDING
	cout << "Sign :" << to_string(sign) << endl;
	cout << "LastFracBit:" << to_string(lastFracBit) << endl;
	cout << "RoundBit :" << to_string(roundBit) << endl;
	cout << "sticky :" << to_string(sticky) << endl;
	cout << static_cast<int>(roundingMode) << endl;
	cout << "Rounding code :" <<  to_string(roundingCode) << endl;
	cout << "Any : " << to_string(any) << endl;
	cout << "Always : " << to_string(always) << endl;
	cout << "tieBreak" << to_string(tieBreak) << endl;
	cout << "Result : " << to_string(result) << endl;
#endif
	return result;
}

#endif // IEEE_ROUNDING_HPP

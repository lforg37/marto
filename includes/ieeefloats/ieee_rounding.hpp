#ifndef IEEE_ROUNDING_HPP
#define IEEE_ROUNDING_HPP

#include "ieeetype.hpp"
/*
#include <iostream>
#include "tools/printing.hpp"

using std::cout;
using std::endl;

using hint::to_string;
*/

template<template<unsigned int, bool> class Wrapper>
Wrapper<1, false> ieee_getRoundBit(
		Wrapper<1, false> sign,
		Wrapper<1, false> lastFracBit,
		Wrapper<1, false> roundBit,
		Wrapper<1, false> sticky,
		IEEERoundingMode roundingMode
		)
{
	/*
	cout << "Sign :" << to_string(sign) << endl;
	cout << "LastFracBit:" << to_string(lastFracBit) << endl;
	cout << "RoundBit :" << to_string(roundBit) << endl;
	cout << "sticky :" << to_string(sticky) << endl;
	*/
	Wrapper<3, false> roundingCode{static_cast<uint8_t>(roundingMode)};
	//cout << static_cast<int>(roundingMode) << endl;
	auto b0 = roundingCode.template get<0>();
	auto b1 = roundingCode.template get<1>();
	auto b2 = roundingCode.template get<2>();

	//cout << "Rounding code :" <<  to_string(roundingCode) << endl;

	auto any = b2 & (b0.invert() | (b1 ^ sign));
	//cout << "Any : " << to_string(any) << endl;

	auto always = b0 & (b1 ^ sign);
	//cout << "Always : " << to_string(always) << endl;

	auto tieBreak = Wrapper<1, false>::mux(b1, sign, lastFracBit);
	//cout << "tieBreak" << to_string(tieBreak) << endl;

	auto result = any & ((sticky & always) | (roundBit & (sticky | tieBreak | always)));
	//cout << "Result : " << to_string(result) << endl;
	return result;
}

#endif // IEEE_ROUNDING_HPP

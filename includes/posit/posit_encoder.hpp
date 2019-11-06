#pragma once

#include "posit_dim.hpp"
#include "primitives/shifter_sticky.hpp"

using namespace hint;
using namespace std;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
PositEncoding<N, WES, Wrapper> posit_encoder(PositIntermediateFormat<N, WES, Wrapper> positValue)
{
	constexpr auto S_WF = PositDim<N, WES>::WF;
	constexpr auto S_WE = PositDim<N, WES>::WE;
	constexpr auto S_WES = WES;
	constexpr auto K_SIZE = S_WE - S_WES;

	auto expWoBias = positValue.getExp().modularSub(Wrapper<S_WE, false>{PositDim<N, WES>::EXP_BIAS});
	//cerr << "expwobias : " << to_string(expWoBias) << endl;

	auto sign = positValue.getSignBit();
	auto sign_sequence_wes = Wrapper<S_WES, false>::generateSequence(sign);

	auto es_wo_xor = expWoBias.template slice<S_WES-1, 0>();
	// TODO: why this does not work when hint<S_WES> es = es_wo_xor xor sign; ???
	// Compiler throws error: error: conversion from ‘ap_int_base<1, false>::RType<1, false>::logic’ {aka ‘ap_uint<1>’} to non-scalar type ‘hint<1>’ {aka ‘hint_base<1, false>’} requested
	auto es = es_wo_xor.bitwise_xor(sign_sequence_wes);

	//cerr << "ES : " << to_string(es) << endl;

	//K_SIZE
	auto k = expWoBias.template slice<S_WE-1, S_WES>();

	//cerr << "range : " << to_string(k) << endl;

	// N - (3 + WES)
	auto significand = positValue.getFraction();

	//cerr << "significand : " << to_string(significand) << endl;

	// N-3
	auto esAndSignificand = es.concatenate(significand);

	Wrapper<2, false> zero_one{1};
	Wrapper<2, false> one_zero{2};

	auto isNegative = k.template get<K_SIZE-1>().bitwise_xor(sign);
	auto leading = Wrapper<2, false>::mux(isNegative, zero_one, one_zero);

	//cerr << "Leading : " << to_string(leading) << endl;

	//N-1
	auto reverseBitAndEsAndSignificand = leading.concatenate(esAndSignificand);

	// K_SIZE - 1
	auto low_k = k.template slice<K_SIZE-2, 0>();
	auto absK = Wrapper<K_SIZE-1, false>::mux(
				k.template get<K_SIZE-1>(),
				low_k.invert(),
				low_k
			);
	//cerr << "AbsK : " << to_string(absK) << endl;

	//N
	auto reverseBitAndEsAndSignificandAndGuardBit = reverseBitAndEsAndSignificand.concatenate(positValue.getGuardBit());

	//cerr << "Before shift : " << to_string(reverseBitAndEsAndSignificandAndGuardBit) << endl;
	// readyToShift.print();

	//N+1
	auto shifted = shifter_sticky(
				reverseBitAndEsAndSignificandAndGuardBit,
				absK,
				reverseBitAndEsAndSignificandAndGuardBit.template get<N-1>()
		); //TODO rajouter le fillbit

	//cerr << "Shifted : " << to_string(shifted) << endl;

	//N-1
	auto unroundedResult =  shifted.template slice<N, 2>();

	//cerr << "Unrounded : " << to_string(unroundedResult) << endl;


	auto guard = shifted.template get<1>();
	auto sticky = shifted.template get<0>().bitwise_or(positValue.getStickyBit());

	//cerr << "guard : " << to_string(guard) << endl << "sticky : " << to_string(sticky) << endl;

	auto roundingBit = guard.bitwise_and(sticky.bitwise_or(unroundedResult.template get<0>()));

	//cerr << "rounding : " << to_string(roundingBit) << endl;

	auto roundedResult = unroundedResult.modularAdd(roundingBit.template leftpad<N-1>());

	auto normalOutput = sign.concatenate(roundedResult);
	auto zero = Wrapper<N-1, false>::generateSequence({0});
	auto isNaRBit = positValue.getIsNaR();
	auto specialCasesValue = isNaRBit.concatenate(zero);

	auto isSpecial = positValue.getSignBit().invert().bitwise_and(positValue.getImplicitBit().invert()).bitwise_or(isNaRBit);
	return Wrapper<N, false>::mux(isSpecial, specialCasesValue, normalOutput);
}

#pragma once

#include <type_traits>
#include "posit_dim.hpp"
#include "primitives/shifter_sticky.hpp"
#include "primitives/shifter.hpp"

#ifdef POSIT_ENCODER_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cout;
#endif

using namespace hint;
using namespace std;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<PositDim<N, WES>::WF + WES, false> buildEsSignifSequence(
		Wrapper<1, false> sign,
		Wrapper<PositDim<N, WES>::WF, false> significand,
		Wrapper<PositDim<N, WES>::WE, false> exponent,
		typename enable_if<PositDim<N, WES>::HAS_ES>::type* = 0
	)
{
	auto sign_sequence_wes = Wrapper<WES, false>::generateSequence(sign);

	auto es_wo_xor = exponent.template slice<WES-1, 0>();
	auto es = es_wo_xor.bitwise_xor(sign_sequence_wes);
	auto ret = es.concatenate(significand);
#ifdef POSIT_ENCODER_DEBUG
	cerr << "=== buildEsSignif (WES!=0) ===" << endl;
	cerr << "sign_sequence_wes: " << to_string(sign_sequence_wes) << endl;
	cerr << "es_wo_xor: " << to_string(es_wo_xor) << endl;
	cerr << "es: " << to_string(es) << endl;
	cerr << "ret: " << to_string(ret) << endl;
	cerr << "==============================" << endl;
#endif
	return ret;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<PositDim<N, WES>::WF, false> buildEsSignifSequence(
		Wrapper<1, false>,
		Wrapper<PositDim<N, WES>::WF, false> significand,
		Wrapper<PositDim<N, WES>::WE, false>,
		typename enable_if<not PositDim<N, WES>::HAS_ES>::type* = 0
	)
{
#ifdef POSIT_ENCODER_DEBUG
	cerr << "=== buildEsSignif (WES==0) ===" << endl;
	cerr << "ret: " << to_string(significand) << endl;
	cerr << "==============================" << endl;
#endif
	return significand;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<1, false> round_would_overflow(
			Wrapper<1, false>,
			Wrapper<PositDim<N, WES>::WE, false>,
			typename enable_if<PositDim<N, WES>::HAS_ES>::type* = 0
		)
{
	#ifdef POSIT_ENCODER_DEBUG
		cerr << "=== round_would_overflow (WES != 0) ===" << endl;
		cerr << "ret: 0" << endl;
		cerr << "=======================================" << endl;
	#endif
		return {0};
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<1, false> round_would_overflow(
			Wrapper<1, false> sign,
			Wrapper<PositDim<N, WES>::WE, false> exp,
			typename enable_if<not PositDim<N, WES>::HAS_ES>::type* = 0
		)
{
	constexpr auto WE = PositDim<N, WES>::WE;
	Wrapper<WE, false> emax{PositDim<N, WES>::EMax}, emin{PositDim<N, WES>::EMin_repr};
	auto isEmin = (exp.as_signed() <= emin.as_signed());
	auto isEmax = (exp.as_signed() >= emax.as_signed());
	auto res = (isEmin & sign) | (isEmax & sign.invert());
	#ifdef POSIT_ENCODER_DEBUG
		cerr << "=== round_would_overflow (WES != 0) ===" << endl;
		cerr << "emax: " << to_string(emax) << endl;
		cerr << "emin: " << to_string(emin) << endl;
		cerr << "isEmin: " << to_string(isEmin) << endl;
		cerr << "isEmax: " << to_string(isEmax) << endl;
		cerr << "res: " << to_string(res) << endl;
		cerr << "=======================================" << endl;
	#endif
	return res;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositEncoding<N, WES, Wrapper> posit_encoder(PositIntermediateFormat<N, WES, Wrapper, false> positValue)
{
	constexpr auto S_WF = PositDim<N, WES>::WF;
	constexpr auto S_WE = PositDim<N, WES>::WE;
	constexpr auto S_WES = WES;
	constexpr auto K_SIZE = S_WE - S_WES;

	auto exp = positValue.getExp();
	auto sign = positValue.getSignBit();
	auto significand = positValue.getFraction();

	//K_SIZE
	auto k = exp.template slice<S_WE-1, S_WES>();

	// N-3
	auto esAndSignificand = buildEsSignifSequence<N, WES, Wrapper>(sign, significand, exp);

	Wrapper<2, false> zero_one{1};
	Wrapper<2, false> one_zero{2};

	auto isNegative = k.template get<K_SIZE-1>().bitwise_xor(sign);
	auto leading = Wrapper<2, false>::mux(isNegative, zero_one, one_zero);

	//N-1
	auto reverseBitAndEsAndSignificand = leading.concatenate(esAndSignificand);

	// K_SIZE - 1
	auto low_k = k.template slice<K_SIZE-2, 0>();
	auto absK = Wrapper<K_SIZE-1, false>::mux(
				k.template get<K_SIZE-1>(),
				low_k.invert(),
				low_k
			);
	//N
	auto reverseBitAndEsAndSignificandAndGuardBit = reverseBitAndEsAndSignificand.concatenate(positValue.getGuardBit());

	// readyToShift.print();

	//N+1
	auto shifted = shifter_sticky(
				reverseBitAndEsAndSignificandAndGuardBit,
				absK,
				reverseBitAndEsAndSignificandAndGuardBit.template get<N-1>()
		); //TODO rajouter le fillbit

	//N-1
	auto unroundedResult =  shifted.template slice<N, 2>();

	auto guard = shifted.template get<1>();
	auto sticky = shifted.template get<0>().bitwise_or(positValue.getStickyBit());

	auto roundOverflow = round_would_overflow<N, WES, Wrapper>(sign, exp);
	auto roundingBit = guard & (sticky | unroundedResult.template get<0>()) & roundOverflow.invert();

	auto roundedResult = unroundedResult.modularAdd(roundingBit.template leftpad<N-1>());

	auto normalOutput = sign.concatenate(roundedResult);
	auto zero = Wrapper<N-1, false>::generateSequence({0});
	auto isNaRBit = positValue.getIsNaR();
	auto specialCasesValue = isNaRBit.concatenate(zero);

	auto isSpecial = (positValue.getSignBit().invert() & positValue.getImplicitBit().invert()) | isNaRBit;

	auto ret = Wrapper<N, false>::mux(isSpecial, specialCasesValue, normalOutput);
#ifdef POSIT_ENCODER_DEBUG
	cerr << "=== POSIT_ENCODER_INEXACT ===" << endl;
	cerr << "exp: " << to_string(exp) << endl;
	cerr << "sign: " << to_string(sign) << endl;
	cerr << "k: " << to_string(k) << endl;
	cerr << "significand: " << to_string(significand) << endl;
	cerr << "esAndSignificand: " << to_string(esAndSignificand) << endl;
	cerr << "isNegative: " << to_string(isNegative) << endl;
	cerr << "leading: " << to_string(leading) << endl;
	cerr << "reverseBitAndEsAndSignificand: " << to_string(reverseBitAndEsAndSignificand) << endl;
	cerr << "low_k: " << to_string(low_k) << endl;
	cerr << "absK: " << to_string(absK) << endl;
	cerr << "reverseBitAndEsAndSignificandAndGuardBit: " << to_string(reverseBitAndEsAndSignificandAndGuardBit) << endl;
	cerr << "shifted: " << to_string(shifted) << endl;
	cerr << "unroundedResult: " << to_string(unroundedResult) << endl;
	cerr << "guard: " << to_string(guard) << endl;
	cerr << "sticky: " << to_string(sticky) << endl;
	cerr << "roundOverflow: " << to_string(roundOverflow) << endl;
	cerr << "roundingBit: " << to_string(roundingBit) << endl;
	cerr << "roundedResult: " << to_string(roundedResult) << endl;
	cerr << "normalOutput: " << to_string(normalOutput) << endl;
	cerr << "zero: " << to_string(zero) << endl;
	cerr << "isNaRBit: " << to_string(isNaRBit) << endl;
	cerr << "specialCasesValue: " << to_string(specialCasesValue) << endl;
	cerr << "isSpecial: " << to_string(isSpecial) << endl;
	cerr << "ret: " << to_string(ret) << endl;
	cerr << "=============================" << endl;
#endif

	return ret;
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositEncoding<N, WES, Wrapper> posit_encoder(PositIntermediateFormat<N, WES, Wrapper, true> positValue)
{
	constexpr auto S_WF = PositDim<N, WES>::WF;
	constexpr auto S_WE = PositDim<N, WES>::WE;
	constexpr auto S_WES = WES;
	constexpr auto K_SIZE = S_WE - S_WES;

	auto exp = positValue.getExp();
	auto significand = positValue.getFraction();
	auto sign = positValue.getSignBit();

	//K_SIZE
	auto k = exp.template slice<S_WE-1, S_WES>();

	// N - (3 + WES)


	// N-3
	auto esAndSignificand = buildEsSignifSequence<N, WES, Wrapper>(sign, significand, exp);

	Wrapper<2, false> zero_one{1};
	Wrapper<2, false> one_zero{2};

	auto isNegative = k.template get<K_SIZE-1>() ^ sign;
	auto leading = Wrapper<2, false>::mux(isNegative, zero_one, one_zero);

	//N-1
	auto reverseBitAndEsAndSignificand = leading.concatenate(esAndSignificand);

	// K_SIZE - 1
	auto low_k = k.template slice<K_SIZE-2, 0>();
	auto absK = Wrapper<K_SIZE-1, false>::mux(
				k.template get<K_SIZE-1>(),
				low_k.invert(),
				low_k
			);

	//N-1
	auto shifted = shifter<true>(
				reverseBitAndEsAndSignificand,
				absK,
				reverseBitAndEsAndSignificand.template get<N-1-1>()
		);

	auto normalOutput = sign.concatenate(shifted);
	auto zero = Wrapper<N-1, false>::generateSequence({0});
	auto isNaRBit = positValue.getIsNaR();
	auto specialCasesValue = isNaRBit.concatenate(zero);

	auto isSpecial = positValue.getSignBit().invert().bitwise_and(positValue.getImplicitBit().invert()).bitwise_or(isNaRBit);
	auto ret = Wrapper<N, false>::mux(isSpecial, specialCasesValue, normalOutput);
#ifdef POSIT_ENCODER_DEBUG
	cerr << "=== POSIT_ENCODER_EXACT ===" << endl;
	cerr << "expWoBias: " << to_string(exp) << endl;
	cerr << "sign: " << to_string(sign) << endl;
	cerr << "k: " << to_string(k) << endl;
	cerr << "significand: " << to_string(significand) << endl;
	cerr << "esAndSignificand: " << to_string(esAndSignificand) << endl;
	cerr << "isNegative: " << to_string(isNegative) << endl;
	cerr << "leading: " << to_string(leading) << endl;
	cerr << "reverseBitAndEsAndSignificand: " << to_string(reverseBitAndEsAndSignificand) << endl;
	cerr << "low_k: " << to_string(low_k) << endl;
	cerr << "absK: " << to_string(absK) << endl;
	cerr << "shifted: " << to_string(shifted) << endl;
	cerr << "normalOutput: " << to_string(normalOutput) << endl;
	cerr << "zero: " << to_string(zero) << endl;
	cerr << "isNaRBit: " << to_string(isNaRBit) << endl;
	cerr << "specialCasesValue: " << to_string(specialCasesValue) << endl;
	cerr << "isSpecial: " << to_string(isSpecial) << endl;
	cerr << "ret: " << to_string(ret) << endl;
	cerr << "==========================" << endl;
#endif
	return ret;
}

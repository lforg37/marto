#pragma once
#include "posit_dim.hpp"
#include "primitives/lzoc_shifter.hpp"

#ifdef POSIT_DECODER_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<PositDim<N, WES>::WE, true> getExponent(
		Wrapper<hint::Static_Val<N-2>::_storage + 1, false> range_count,
		Wrapper<N-3, false> shifted_fraction,
		Wrapper<1, false> sign,
		typename enable_if<PositDim<N, WES>::HAS_ES>::type* = 0
	)
{
	auto es = shifted_fraction.template slice<N - 4, N-3-WES>();
	auto ext_sign = Wrapper<WES, false>::generateSequence(sign);

	auto decoded_es = es.bitwise_xor(ext_sign);

	return range_count.concatenate(decoded_es).as_signed();
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline Wrapper<PositDim<N, WES>::WE, true> getExponent(
		Wrapper<hint::Static_Val<N-2>::_storage + 1, false> range_count,
		Wrapper<N-3, false>,
		Wrapper<1, false>,
		typename enable_if<not PositDim<N, WES>::HAS_ES>::type* = 0
	)
{
	return range_count.as_signed();
}

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, true> posit_decoder(PositEncoding<N, WES, Wrapper> positN)
{
	//Sign bit
	auto s = positN.template get<N-1>();
	//First regime bit
	auto count_type = positN.template get<N-2>();
	//Remainder
	auto input_shift = positN.template slice<N-3, 0>();

	auto zero_NAR = input_shift.or_reduction().bitwise_or(count_type).invert();
	auto is_NAR = zero_NAR.bitwise_and(s);
	auto is_zero = zero_NAR.bitwise_and(s.invert());

	auto implicit_bit = s.invert().bitwise_and(zero_NAR.invert());

	/*constexpr int logN = Static_Val<N>::_log2;
	auto lzoc_shifted = lzoc_shifter<logN>(input_shift, count_type);
	ap_uint<logN + 1> RangeCount = lzoc_shifted.range(N+logN-1, N);
	ap_uint<N-3> usefulBits = lzoc_shifted.range(N-2, 2);
	ap_uint<N-3-WES> fraction = usefulBits(N-4-WES, 0);

	ap_int<1> neg_count = not(s^count_type);
	ap_int<logN + 1> signed_extended_neg_count = neg_count;
	ap_uint<logN + 1> extended_neg_count = signed_extended_neg_count;
	ap_uint<logN + 1> comp2_range_count = RangeCount^extended_neg_count;

	ap_uint<PositIntermediateFormat<N, WES>::ExpSize> exponent = getExponent<N, WES>(
			comp2_range_count,
			usefulBits,
			s
		);*/

	constexpr int rangeCountSize = hint::Static_Val<N-2>::_storage;
	auto lzoc_shifted = hint::LZOC_shift<N-2, N-2>(input_shift, count_type);

	auto rangeCount = lzoc_shifted.lzoc.template leftpad<rangeCountSize+1>();
	auto usefulBits = lzoc_shifted.shifted.template slice<N-4, 0>();
	auto fraction = usefulBits.template slice<N-4-WES, 0>();

	auto neg_count = s.bitwise_xor(count_type).invert();
	auto extended_neg_count = Wrapper<rangeCountSize + 1, false>::generateSequence(neg_count);
	auto comp2_range_count = rangeCount.bitwise_xor(extended_neg_count);

	auto exponent = getExponent<N, WES>(
				comp2_range_count,
				usefulBits,
				s
			);

	//auto biased_exp = exponent.modularAdd({PositDim<N, WES>::EXP_BIAS});

	auto is_not_zero = is_zero.invert();
	//auto extended_is_not_zero = Wrapper<PositDim<N, WES>::WE, false>::generateSequence(is_not_zero);
	//auto final_biased_exp = biased_exp.bitwise_and(extended_is_not_zero);
	auto final_biased_exp = exponent.as_unsigned();//.bitwise_and(extended_is_not_zero);

#ifdef POSIT_DECODER_DEBUG
	cerr << "=== POSIT_DECODER ===" << endl;
	cerr << "positN: " << to_string(positN) << endl;
	cerr << "s: " << to_string(s) << endl;
	cerr << "count_type: " << to_string(count_type) << endl;
	cerr << "input_shift: " << to_string(input_shift) << endl;
	cerr << "zero_NAR: " << to_string(zero_NAR) << endl;
	cerr << "is_NAR: " << to_string(is_NAR) << endl;
	cerr << "is_zero: " << to_string(is_zero) << endl;
	cerr << "implicit_bit: " << to_string(implicit_bit) << endl;
	cerr << "lzoc: " << to_string(lzoc_shifted.lzoc) << endl;
	cerr << "shifted: " << to_string(lzoc_shifted.shifted) << endl;
	cerr << "rangeCount: " << to_string(rangeCount) << endl;
	cerr << "usefulBits: " << to_string(usefulBits) << endl;
	cerr << "fraction: " << to_string(fraction) << endl;
	cerr << "neg_count: " << to_string(neg_count) << endl;
	cerr << "extended_neg_count: " << to_string(extended_neg_count) << endl;
	cerr << "comp2_range_count: " << to_string(comp2_range_count) << endl;
	cerr << "exponent: " << to_string(exponent) << endl;
	cerr << "is_not_zero: " << to_string(is_not_zero) << endl;
	cerr << "extended_is_not_zero: " << to_string(extended_is_not_zero) << endl;
	cerr << "final_biased_exp: " << to_string(final_biased_exp) << endl;
	cerr << "=====================" << endl;
#endif

	return PositIntermediateFormat<N, WES, Wrapper, true>(
				is_NAR,
				final_biased_exp,
				s,
				implicit_bit,
				fraction
			);
}

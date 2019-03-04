#ifndef LZOC_SHIFTER_HPP
#define LZOC_SHIFTER_HPP
#include <iostream>
#include <type_traits>

#include "ap_int.h"

template<int S>
struct LZOCStageInfo
{
	static constexpr bool NeedsRecursion = (S>0);
	static constexpr bool IsFinalStage = (S==0);
};

template<int N, int S>
//N : Power of 2 of the size of the whole LZOC, 
//S power of two of the size of the stage
inline ap_uint<S + 1 + (1 << N)> lzoc_shifter_stage(
		ap_uint<1<<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<LZOCStageInfo<S>::NeedsRecursion>::type* dummy = 0
	)
{
	#pragma HLS INLINE
	ap_uint<1<<S> zeros = 0;
	ap_uint<1<<S> ones = -1;
	ap_int<1<<S> padding_s = (ap_int<1>) fill_bit;
	ap_uint<1<<S> padding = padding_s;

	ap_uint<1 << S> high = input.range((1 << N) - 1, (1 << N) - (1 << S));
	ap_uint<(1 << N) - (1 << S)> low = input.range((1 << N) - (1 << S) - 1, 0); 

	ap_uint<1> leader;
	ap_uint<1<<N> next_stage_input;

	if ((leading && (high == ones)) || (!leading && high == zeros) ) {
		next_stage_input = low.concat(padding);
		leader = 1;
	} else {
		next_stage_input = input;
		leader = 0;
	}
	auto lower_stage = lzoc_shifter_stage<N, S-1>(next_stage_input, leading, fill_bit);
	return leader.concat(lower_stage);
}

template<int N, int S>
inline ap_uint<S + 1 + (1 << N)> lzoc_shifter_stage(
		ap_uint<1<<N> input,
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<LZOCStageInfo<S>::IsFinalStage>::type* dummy = 0
	)
{
	#pragma HLS INLINE
	if (input[(1<<N) - 1] == leading) {
		ap_uint<(1<<N) - 1> low = input.range((1<<N) - 2, 0);
		ap_uint<(1<<N)> res = low.concat(fill_bit);
		return ap_uint<1>(1).concat(res);
	} else {
		return ap_uint<1>(0).concat(input);
	}
}

template<int N>
ap_uint<N + (1<<N)> lzoc_shifter(
		ap_uint<1<<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0)
{
	#pragma HLS INLINE
	return lzoc_shifter_stage<N, N-1>(input, leading, fill_bit);
}
#endif

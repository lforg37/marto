#ifndef LZOC_SHIFTER_HPP
#define LZOC_SHIFTER_HPP
#include <iostream>

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
		typename enable_if<LZOCStageInfo<S>::NeedsRecursion>::type* dummy = 0
	)
{
	ap_uint<1<<S> zeros = 0;
	ap_uint<1<<S> ones = -1;
	ap_uint<1<<S> padding = 0;

	ap_uint<1 << S> high = input.range((1 << N) - 1, (1 << N) - (1 << S));
	ap_uint<(1 << N) - (1 << S)> low = input.range((1 << N) - (1 << S) - 1, 0); 


	if ((leading && (high == ones)) || (!leading && high == zeros) ) {
		ap_uint<1<<N> next_stage_input = low.concat(padding);
		auto lower_stage = lzoc_shifter_stage<N, S-1>(next_stage_input, leading);
		return ap_uint<1>(1).concat(lower_stage);
	} else {
		auto lower_stage = lzoc_shifter_stage<N, S-1>(input, leading);
		return ap_uint<1>(0).concat(lower_stage);
	}
}

template<int N, int S>
inline ap_uint<S + 1 + (1 << N)> lzoc_shifter_stage(
		ap_uint<1<<N> input,
		ap_uint<1> leading,
		typename enable_if<LZOCStageInfo<S>::IsFinalStage>::type* dummy = 0
	)
{
	if (input[(1<<N) - 1] == leading) {
		ap_uint<(1<<N) - 1> low = input.range((1<<N) - 2, 0);
		ap_uint<(1<<N)> res = low.concat(ap_uint<1>(0));
		return ap_uint<1>(1).concat(res);
	} else {
		return ap_uint<1>(0).concat(input);
	}
}

template<int N>
ap_uint<N + (1<<N)> lzoc_shifter(ap_uint<1<<N> input, ap_uint<1> leading)
{
	return lzoc_shifter_stage<N, N-1>(input, leading);
}
#endif

#ifndef SHIFTER_HPP
#define SHIFTER_HPP
#include <iostream>
#include <type_traits>

#include "ap_int.h"


template<int S>
struct ShifterStageInfo
{
	static constexpr bool NeedsRecursion = (S>1);
	static constexpr bool IsFinalStage = (S==1);
};

template<int N, int S>
//N : Power of 2 of the size of the whole LZOC, 
//S : size of the shift to consider
inline ap_uint<(1 << N)> shifter_stage(
		ap_uint<1<<N> input, 
		ap_uint<S> count,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<ShifterStageInfo<S>::NeedsRecursion>::type* dummy = 0
	)
{
	ap_int<1<<(S-1)> padding_s = (ap_int<1>) fill_bit;
	ap_uint<1<<(S-1)> padding = padding_s;
	ap_uint<1> stageNeedsShift = count[S-1];
	ap_uint<S-1> countnext = count.range(S-2, 0);

	ap_uint<(1 << N) - (1 << (S-1))> low = input.range((1 << N) - (1 << (S-1)) - 1, 0); 
	ap_uint<1<<N> next_stage_input;

	if (stageNeedsShift) {
		next_stage_input = low.concat(padding);
	} else {
		next_stage_input = input;
	}
	return shifter_stage<N, S-1>(next_stage_input, countnext, fill_bit);
}

template<int N, int S>
inline ap_uint<(1 << N)> shifter_stage(
		ap_uint<1<<N> input,
		ap_uint<S> count,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<ShifterStageInfo<S>::IsFinalStage>::type* dummy = 0
	)
{
	ap_uint<1<<N> result;
	if (count[0] == 1) {
		ap_uint<(1<<N) - 1> low = input.range((1<<N) - 2, 0);
		result = low.concat(fill_bit);
	} else {
		result = input;
	}
	return result;
}

template<int N>
ap_uint<(1<<N)> shifter(
		ap_uint<1<<N> input, 
		ap_uint<N-1> count,
		ap_uint<1> fill_bit = 0)
{
	return shifter_stage<N, N-1>(input, count, fill_bit);
}
#endif

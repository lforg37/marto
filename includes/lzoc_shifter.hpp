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

	ap_uint<1<<S> padding;
	if(fill_bit){
		padding = -1;
	}
	else{
		padding = 0;
	}

	ap_uint<(1 << N) - (1 << S)> low = input.range((1 << N) - (1 << S) - 1, 0); 


	ap_uint<1> cmp = 1;
	for(int i = (1 << N) - 1; i>=((1 << N) - (1 << S)); i--){
		#pragma HLS UNROLL
		cmp &= (input[i]==leading);
	}

	ap_uint<1<<N> next_stage_input = (cmp) ? low.concat(padding) : input;
	ap_uint<1> leader = (cmp) ? 1 : 0;

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







template<int N>
struct GenericLZOCStageInfo
{
	static constexpr bool is_odd = ((N%2) == 1);
	static constexpr bool is_pair = ((N%2) == 0);
};




template<int N, int S>
ap_uint<Static_Val<N>::_log2 + N> generic_lzoc_shifter_stage(
		ap_uint<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<(S==1)>::type* dummy = 0)
{
	#pragma HLS INLINE
	if (input[N - 1] == leading) {
		ap_uint<N - 1> low = input.range(N - 2, 0);
		ap_uint<N> res = low.concat(fill_bit);
		return ap_uint<1>(1).concat(res);
	} else {
		return ap_uint<1>(0).concat(input);
	}

}



template<int N, int S>
ap_uint<Static_Val<N>::_log2 + N> generic_lzoc_shifter_stage(
		ap_uint<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<(S>1)>::type* dummy = 0)
{
	#pragma HLS INLINE
	static constexpr int half_size = Static_Ceil_Div<S,2>::val;

	ap_uint<S> padding;
	if(fill_bit){
		padding = -1;
	}
	else{
		padding = 0;
	}

	ap_uint<N-S> low = input.range(N-S - 1, 0); 


	ap_uint<1> cmp = 1;
	for(int i = N - 1; i>=(N-S); i--){
		#pragma HLS UNROLL
		cmp &= (input[i]==leading);
	}

	ap_uint<N> next_stage_input = (cmp) ? low.concat(padding) : input;
	ap_uint<Static_Val<S>::_log2> stage_lzoc = (cmp) ? S : 0;

	ap_uint<Static_Val<N>::_log2 + N> lower_stage = generic_lzoc_shifter_stage<N, half_size>(next_stage_input, leading, fill_bit);
	ap_uint<Static_Val<N>::_log2> lower_stage_lzoc = lower_stage.range(Static_Val<N>::_log2 + N-1, N);
	ap_uint<N> lower_stage_shift = lower_stage.range(N-1, 0);

	ap_int<Static_Val<N>::_log2> lzoc = stage_lzoc + lower_stage_lzoc;
	return lzoc.concat(lower_stage_shift);
}


template<int N>
ap_uint<Static_Val<N>::_log2 + N> generic_lzoc_shifter(
		ap_uint<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<GenericLZOCStageInfo<N>::is_pair>::type* dummy = 0)
{
	#pragma HLS INLINE
	static constexpr int half_size = Static_Ceil_Div<N,2>::val;

	return generic_lzoc_shifter_stage<N, half_size>(input, leading, fill_bit);
}


template<int N>
ap_uint<Static_Val<N>::_log2 + N> generic_lzoc_shifter(
		ap_uint<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<GenericLZOCStageInfo<N>::is_odd>::type* dummy = 0)
{
	#pragma HLS INLINE
	static constexpr int pair_size = N-1;

	ap_uint<Static_Val<N>::_log2> lzoc;
	ap_uint<N> shift;

	ap_uint<Static_Val<N-1>::_log2+ N-1> pair_lzoc_shift;
	if(input[N-1] == leading){

		ap_uint<pair_size> shrinked_input = input.range(N-2,0);
		pair_lzoc_shift = generic_lzoc_shifter<pair_size>(shrinked_input, leading, fill_bit);
		lzoc = ((ap_uint<Static_Val<N>::_log2>)pair_lzoc_shift.range(Static_Val<N-1>::_log2+N-1-1, N-1)) + 1;
		ap_uint<N-1> pair_shift = pair_lzoc_shift.range(N-1-1, 0);
		shift = pair_shift.concat(fill_bit);

	}
	else{
		shift = input;
		lzoc = 0;
	}
	return  lzoc.concat(shift);
}



#endif

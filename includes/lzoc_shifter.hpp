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
	static constexpr bool is_a_power_of_2 = ((Static_Val<N>::_2pow) == N);
	static constexpr bool is_one = (N==1);
};

template<int N, int S>
ap_uint<Static_Val<S>::_rlog2 + N> generic_lzoc_shifter_stage(
		ap_uint<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<GenericLZOCStageInfo<S>::is_a_power_of_2 and GenericLZOCStageInfo<S>::is_one>::type* dummy = 0)
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
ap_uint<Static_Val<S>::_rlog2 + N> generic_lzoc_shifter_stage(
		ap_uint<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<not(GenericLZOCStageInfo<S>::is_a_power_of_2)>::type* dummy = 0)
{
	#pragma HLS INLINE
	static constexpr int log2S = Static_Val<S>::_rlog2-1;
	static constexpr int power_of_2_in_S = 1<<log2S;
	static constexpr int rest_of_S = S-power_of_2_in_S;

// fprintf(stderr, "N: %d, S: %d, log2s: %d, powerlo2s: %d, rest: %d \n", N, S, log2S, power_of_2_in_S, rest_of_S);


	ap_uint<Static_Val<S>::_rlog2 + N> lzoc_shift = generic_lzoc_shifter_stage<N, power_of_2_in_S>( input, leading, fill_bit);

	ap_uint<Static_Val<S>::_rlog2> lzoc = lzoc_shift.range(Static_Val<S>::_rlog2 + N-1, N);
	ap_uint<N> shift = lzoc_shift.range(N-1, 0);

	ap_uint<Static_Val<rest_of_S>::_rlog2 + N> lzoc_shift_rest = generic_lzoc_shifter_stage<N, rest_of_S>( shift, leading, fill_bit);


	ap_uint<Static_Val<rest_of_S>::_rlog2> lzoc_rest = lzoc_shift_rest.range(Static_Val<rest_of_S>::_rlog2 + N-1, N);
	ap_uint<N> shift_rest = lzoc_shift_rest.range(N-1, 0);


	ap_uint<N> final_shift = (lzoc==-1) ? shift_rest : shift;
	ap_uint<Static_Val<S>::_rlog2> lzoc_if_rest = S+lzoc_rest;

	ap_uint<Static_Val<S>::_rlog2> final_lzoc = (lzoc==-1) ? lzoc_if_rest : lzoc;

	return final_lzoc.concat(final_shift);

/*	ap_uint<power_of_2_in_S> padding;
	if(fill_bit){
		padding = -1;
	}
	else{
		padding = 0;
	}

	ap_uint<N-power_of_2_in_S> low = input.range(N - power_of_2_in_S - 1, 0); 


	ap_uint<1> cmp = 1;
	for(int i = N - 1; i>=(N - power_of_2_in_S); i--){
		#pragma HLS UNROLL
		cmp &= (input[i]==leading);
	}

	ap_uint<1<<N> next_stage_input = (cmp) ? low.concat(padding) : input;
	ap_uint<1> leader = (cmp) ? 1 : 0;

	auto lower_stage = lzoc_shifter_stage<N, (S>>1)>(next_stage_input, leading, fill_bit);

	return lower_stage;
*/
}	

template<int N, int S>
ap_uint<Static_Val<S>::_rlog2 + N> generic_lzoc_shifter_stage(
		ap_uint<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<GenericLZOCStageInfo<S>::is_a_power_of_2  and not(GenericLZOCStageInfo<S>::is_one)>::type* dummy = 0)
{
	#pragma HLS INLINE
	ap_uint<S> padding;
	if(fill_bit){
		padding = -1;
	}
	else{
		padding = 0;
	}

	ap_uint<N-S> low = input.range(N - S - 1, 0); 


	ap_uint<1> cmp = 1;
	for(int i = N - 1; i>=(N - S); i--){
		#pragma HLS UNROLL
		cmp &= (input[i]==leading);
	}

	ap_uint<N> next_stage_input = (cmp) ? low.concat(padding) : input;
	ap_uint<1> leader = (cmp) ? 1 : 0;

	auto lower_stage = generic_lzoc_shifter_stage<N, (S>>1) >(next_stage_input, leading, fill_bit);
	return leader.concat(lower_stage);
}	




template<int N>
ap_uint<Static_Val<N>::_rlog2 + N> generic_lzoc_shifter(
		ap_uint<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<not(GenericLZOCStageInfo<N>::is_a_power_of_2)>::type* dummy = 0
		)
{
	#pragma HLS INLINE

	ap_uint<(Static_Val<N>::_rlog2 + N)> lzoc_shift =  generic_lzoc_shifter_stage<N, N>(input, leading, fill_bit);
	return lzoc_shift;
}	

template<int N>
ap_uint<Static_Val<N>::_rlog2 + N> generic_lzoc_shifter(
		ap_uint<N> input, 
		ap_uint<1> leading,
		ap_uint<1> fill_bit = 0,
		typename std::enable_if<GenericLZOCStageInfo<N>::is_a_power_of_2>::type* dummy = 0
		)
{
	#pragma HLS INLINE
	static constexpr int log2N = Static_Val<N>::_log2;

	ap_uint<(log2N + (1<<log2N))> lzoc_shift =  lzoc_shifter<log2N>(input, leading, fill_bit);
	return lzoc_shift;
}	


#endif

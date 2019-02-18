#pragma once
#include "ap_int.h"
#include "posit_dim.hpp"


template<int N>
Quire<N> add_sub_quire(
		Quire<N> quire, 
		PositProd<N> input,
	   	ap_uint<1> isSub
){
	
	ap_int<PositDim<N>::ProdSignificandSize> inputSignificand = input.getSignificand();
	ap_int<PositDim<N>::ProdSignificandSize> complementedInputIfIsSub;
	
	#pragma HLS UNROLL
	for (int i=0; i<PositDim<N>::ProdSignificandSize; i++){
		complementedInputIfIsSub[i] = input[i] ^ isSub;
	}
	
	// fprintf(stderr, "=== original ===\n");
	// printApInt(inputSignificand);
	// fprintf(stderr, "=== complemented ===\n");
	// printApInt(complementedInputIfIsSub);


	ap_uint<PositDim<N>::ProdExpSize> shiftValue = input.getExp();
	ap_int<PositDim<N>::ExtQuireSize-1> shiftedInput = (ap_int<PositDim<N>::ExtQuireSize-1>)complementedInputIfIsSub<<(shiftValue);
	fprintf(stderr, "=== shiftdValue ===\n");
	// printApUint(shiftValue);

	// fprintf(stderr, "=== shiftedInput ===\n");
	// printApInt(shiftedInput);

	ap_uint<PositDim<N>::ExtQuireSize-1> quireWithoutSignAndNARBit = quire.getQuireWithoutNaR();

	ap_uint<PositDim<N>::ExtQuireSize-1> sumResult = shiftedInput + quireWithoutSignAndNARBit + isSub;

	ap_uint<1> resultIsNaR = quire.getIsNaR() || input.getIsNaR();

	return Quire<N>(resultIsNaR.concat(sumResult));
}

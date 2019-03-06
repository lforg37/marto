#pragma once
#include <cstdio>

#include "ap_int.h"
#include "posit_dim.hpp"
#include "lzoc_shifter.hpp"
#include "utils.hpp"


#define S_WF PositDim<N>::WF
#define S_WE PositDim<N>::WE
#define S_WES PositDim<N>::WES
#define K_SIZE (S_WE-S_WES)

template<int N>
PositValue<N> posit_add(
		PositValue<N> in1, 
		PositValue<N> in2, 
		ap_uint<1> isSub
){
	#pragma HLS INLINE
	static constexpr int EXT_SUM_SIZE = Static_Val<S_WF+2 + S_WF +1>::_2pow;
	static constexpr int LOG2_EXT_SUM_SIZE = Static_Val<EXT_SUM_SIZE>::_log2;
	
	bool in1IsGreater = in1.getExp() > in2.getExp();

	ap_uint<S_WE> subExpOp1, subExpOp2;
	ap_uint<S_WE+1> shiftValue;
	ap_int<S_WF+2> mostSignificantSignificand, lessSignificantSignificand;

	ap_int<S_WF+2> input2Significand = in2.getSignedSignificand();
	ap_int<S_WF+2> complementedInputIfIsSub;

	for (int i=0; i<S_WF+2; i++){
	#pragma HLS UNROLL
		complementedInputIfIsSub[i] = input2Significand[i] ^ isSub;
	}

	if(in1IsGreater){
		subExpOp1 = in1.getExp();
		subExpOp2 = in2.getExp();
		mostSignificantSignificand = in1.getSignedSignificand();
		lessSignificantSignificand = complementedInputIfIsSub;
	}
	else{
		subExpOp1 = in2.getExp();
		subExpOp2 = in1.getExp();
		mostSignificantSignificand = complementedInputIfIsSub;
		lessSignificantSignificand = in1.getSignedSignificand();	
	}

	ap_uint<S_WF> WFZeros=0;
	ap_uint<S_WF> WFOnes=-1;
	ap_uint<S_WF> toConcatMost = (isSub and not(in1IsGreater))? WFOnes: WFZeros;
	ap_uint<S_WF> toConcatLess = (isSub and in1IsGreater)? WFOnes: WFZeros;

	shiftValue = subExpOp1 - subExpOp2;
	
	ap_int<S_WF+2 + S_WF> shiftedSignificand = ((ap_int<S_WF+2 + S_WF>)lessSignificantSignificand.concat(toConcatLess)) >> shiftValue;
	ap_int<S_WF+2 + S_WF> unShiftedSignificand = mostSignificantSignificand.concat(toConcatMost);



	ap_int<S_WF+2 + S_WF +1 +1> sum = shiftedSignificand + unShiftedSignificand + isSub;

	// printApInt(shiftedSignificand);
	// printApInt(unShiftedSignificand);
	// printApInt(sum);


	ap_uint<(1<<LOG2_EXT_SUM_SIZE)> extSum = ((ap_uint<(1<<LOG2_EXT_SUM_SIZE)>) sum) << ((1<<LOG2_EXT_SUM_SIZE) - (S_WF+2 + S_WF +1)-1);



	ap_uint<1> extsumSign = extSum[EXT_SUM_SIZE -1];
	ap_uint<(LOG2_EXT_SUM_SIZE + (1<<LOG2_EXT_SUM_SIZE))> lzocShifter = lzoc_shifter<LOG2_EXT_SUM_SIZE>(extSum, extsumSign);

	ap_uint<LOG2_EXT_SUM_SIZE> lzoc = lzocShifter.range(LOG2_EXT_SUM_SIZE + EXT_SUM_SIZE-1,EXT_SUM_SIZE);
	ap_uint<EXT_SUM_SIZE> shiftedSum = lzocShifter.range(EXT_SUM_SIZE-1,0);



	// We add two bits to check for both overflows and negative exponents
	ap_uint<S_WE +1 +1> computedExp = subExpOp1 +1 - (lzoc-2);
	ap_uint<1> expIsNegative = computedExp[S_WE +1 +1 -1];
	ap_uint<1> expOverflowed = computedExp[S_WE +1 +1 -1 -1] == 1;


	ap_uint<S_WF+1> resultSignificand = shiftedSum.range(EXT_SUM_SIZE-1,EXT_SUM_SIZE-1 -(S_WF+1)+1);
	ap_uint<EXT_SUM_SIZE -(S_WF+1)> resultRest = shiftedSum.range(EXT_SUM_SIZE-1 -(S_WF+1),0);



	ap_uint<1> guardBit = resultRest[EXT_SUM_SIZE -(S_WF+1)-1];
	ap_uint<1> stickyBit = !(resultRest.range(EXT_SUM_SIZE -(S_WF+1) -1-1, 0) == 0);

	ap_uint<1> resultIsNaR = in1.getIsNaR() || in2.getIsNaR();
	ap_uint<1> isZero = ((resultSignificand == 0) && ((guardBit == 0) || ((guardBit == 1) && (stickyBit == 0)))) && !(extsumSign);
	ap_uint<1> resultS =  (isZero) ? 0 : (!resultSignificand[S_WF+1 -1]); 

	ap_uint<S_WE> resultExp = (isZero) ? 0 : computedExp.range(S_WE-1,0);

	PositValue<N> result = PositValue<N>(
				guardBit,
				stickyBit,
				resultIsNaR,
				resultExp,
				resultS,
				resultSignificand[S_WF+1 -1],
				resultSignificand.range(S_WF+1 -1 -1, 0));


	return result;


}

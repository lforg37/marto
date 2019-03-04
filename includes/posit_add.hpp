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
		PositValue<N> in2
){
	#pragma HLS INLINE
	static constexpr int EXT_SUM_SIZE = Static_Val<S_WF+2 + S_WF +1>::_2pow;
	static constexpr int LOG2_EXT_SUM_SIZE = Static_Val<EXT_SUM_SIZE>::_log2;
	
	bool in1IsGreater = in1.getExp() > in2.getExp();

	ap_uint<S_WE> subExpOp1, subExpOp2;
	ap_uint<S_WE+1> shiftValue;
	ap_int<S_WF+2> mostSignificantSignificand, lessSignificantSignificand;


	if(in1IsGreater){
		subExpOp1 = in1.getExp();
		subExpOp2 = in2.getExp();
		mostSignificantSignificand = in1.getSignedSignificand();
		lessSignificantSignificand = in2.getSignedSignificand();
	}
	else{
		subExpOp1 = in2.getExp();
		subExpOp2 = in1.getExp();
		mostSignificantSignificand = in2.getSignedSignificand();
		lessSignificantSignificand = in1.getSignedSignificand();		
	}

	shiftValue = subExpOp1 - subExpOp2;
	
	ap_uint<S_WF> WFZeros;
	ap_int<S_WF+2 + S_WF> shiftedSignificand = ((ap_int<S_WF+2 + S_WF>)lessSignificantSignificand.concat(WFZeros)) >> shiftValue;
	ap_int<S_WF+2 + S_WF> unShiftedSignificand = mostSignificantSignificand.concat(WFZeros);



	ap_int<S_WF+2 + S_WF +1 +1> sum = shiftedSignificand + unShiftedSignificand;



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

	ap_uint<1> resultIsNaR = in1.getIsNaR() || in1.getIsNaR();
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

extern template PositValue<8> posit_add(PositValue<8> in1, PositValue<8> in2);
extern template PositValue<16> posit_add(PositValue<16> in1, PositValue<16> in2);
extern template PositValue<32> posit_add(PositValue<32> in1, PositValue<32> in2);
extern template PositValue<64> posit_add(PositValue<64> in1, PositValue<64> in2);

#pragma once
#include "ap_int.h"
#include "posit_dim.hpp"
#include "lzoc_shifter.hpp"

#define EXT_SUM_SIZE ceil2Power(PositDim<N>::WF+1 + PositDim<N>::WF +1)
#define LOG2_EXT_SUM_SIZE ceilLog2(EXT_SUM_SIZE)

template<int N>
PositValue<N> posit_add(
		PositValue<N> in1, 
		PositValue<N> in2
){

	bool in1IsGreater = in1.getExp() > in2.getExp();

	ap_uint<PositDim<N>::WE> subExpOp1, subExpOp2;
	ap_uint<ceilLog2(PositDim<N>::WF)> shiftValue;
	ap_uint<PositDim<N>::WF+1> mostSignificantSignificand, lessSignificantSignificand;

	if(in1IsGreater){
		subExpOp1 = in1.getExp();
		subExpOp2 = in2.getExp();
		mostSignificantSignificand = in1.getSignificand();
		lessSignificantSignificand = in2.getSignificand();
	}
	else{
		subExpOp1 = in2.getExp();
		subExpOp2 = in1.getExp();
		mostSignificantSignificand = in2.getSignificand();
		lessSignificantSignificand = in1.getSignificand();		
	}

	shiftValue = subExpOp1 - subExpOp2;
	
	ap_uint<PositDim<N>::WF> WfZeros;
	ap_int<PositDim<N>::WF+1 + PositDim<N>::WF> shiftedSignificand = ((ap_int<PositDim<N>::WF+1 + PositDim<N>::WF>)lessSignificantSignificand.concat(WfZeros)) >> shiftValue;
	ap_int<PositDim<N>::WF+1 + PositDim<N>::WF> unShiftedSignificand = mostSignificantSignificand.concat(WfZeros);

	ap_int<PositDim<N>::WF+1 + PositDim<N>::WF +1> sum = shiftedSignificand + unShiftedSignificand;

	ap_uint<(1<<LOG2_EXT_SUM_SIZE)> extSum = sum << (EXT_SUM_SIZE - PositDim<N>::WF+1 + PositDim<N>::WF +1);
	ap_uint<(LOG2_EXT_SUM_SIZE + (1<<LOG2_EXT_SUM_SIZE))> lzocShifter = lzoc_shifter<LOG2_EXT_SUM_SIZE>(extSum, extSum[EXT_SUM_SIZE -1]);

	ap_uint<LOG2_EXT_SUM_SIZE> lzoc = lzocShifter.range(LOG2_EXT_SUM_SIZE + EXT_SUM_SIZE-1,EXT_SUM_SIZE);
	ap_uint<EXT_SUM_SIZE> shiftedSum = lzocShifter.range(EXT_SUM_SIZE-1,0);

	ap_uint<PositDim<N>::WF+1> resultSignificand = shiftedSum.range(EXT_SUM_SIZE-1,EXT_SUM_SIZE -PositDim<N>::WF+1);
	ap_uint<EXT_SUM_SIZE -PositDim<N>::WF+1> resultRest = shiftedSum.range(EXT_SUM_SIZE -PositDim<N>::WF+1 -1,0);

	ap_uint<1> firstRestBit = resultRest[EXT_SUM_SIZE -PositDim<N>::WF+1 -1];
	ap_uint<1> remainingRestBitsAreZeros = resultRest.range(EXT_SUM_SIZE -PositDim<N>::WF+1 -1 -1, 0) == 0;

	ap_uint<1> roudingBit = (firstRestBit && !remainingRestBitsAreZeros) || (firstRestBit && remainingRestBitsAreZeros && (resultSignificand[0]==1));
	// What if this sum overflows?
	ap_uint<PositDim<N>::WF+1> roundedResultSignificand = resultSignificand + roudingBit;

	ap_uint<PositDim<N>::WE +1 +1> computedExp = subExpOp1 +1 - lzoc;
	ap_uint<1> expIsNegative = computedExp[PositDim<N>::WE +1 +1 -1];
	ap_uint<1> expOverflowed = computedExp[PositDim<N>::WE +1 +1 -1 -1] == 1;

	ap_uint<1> resultIsNaR = in1.getIsNaR() || in1.getIsNaR();
	ap_uint<1> resultS =  (roundedResultSignificand == 0) ? 1 : (!roundedResultSignificand[PositDim<N>::WF+1 -1]); 

	if (expIsNegative){
		// return minpos;
	}
	else if (expOverflowed){
		// return maxpos;
	}
	else{
		return PositValue<N>(	resultIsNaR,	
								(roundedResultSignificand == 0) ? 1 : computedExp.range(PositDim<N>::WE-1,0),
								resultS,
								roundedResultSignificand[PositDim<N>::WF+1 -1],
								roundedResultSignificand.range(PositDim<N>::WF+1 -1 -1,0));
		
	}

}

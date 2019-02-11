#pragma once
#include "ap_int.h"
#include "posit_dim.hpp"
#include "lzoc_shifter.hpp"
#include "utils.hpp"


#define S_WF PositDim<N>::WF
#define S_WE PositDim<N>::WE
#define S_WES PositDim<N>::WES
#define K_SIZE (S_WE-S_WES)
#define EXT_SUM_SIZE ceil2Power(S_WF+1 + S_WF +1)
#define LOG2_EXT_SUM_SIZE ceilLog2(EXT_SUM_SIZE)


#define DEBUG_ADDER


template<int N>
ap_uint<N> posit_add(
		PositValue<N> in1, 
		PositValue<N> in2
){
	bool in1IsGreater = in1.getExp() > in2.getExp();

	ap_uint<S_WE> subExpOp1, subExpOp2;
	ap_uint<S_WE+1> shiftValue;
	ap_uint<S_WF+1> mostSignificantSignificand, lessSignificantSignificand;

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
	
	ap_uint<S_WF> WFZeros;
	ap_uint<S_WF+1 + S_WF> shiftedSignificand = ((ap_int<S_WF+1 + S_WF>)lessSignificantSignificand.concat(WFZeros)) >> shiftValue;
	ap_uint<S_WF+1 + S_WF> unShiftedSignificand = mostSignificantSignificand.concat(WFZeros);

#ifdef DEBUG_ADDER
	fprintf(stderr, "sum operand 1: \t  ");
	printApUint(shiftedSignificand);

	fprintf(stderr, "sum operand 2: \t  ");
	printApUint(unShiftedSignificand);
#endif

	ap_uint<S_WF+1 + S_WF +1 +1> sum = shiftedSignificand + unShiftedSignificand;

#ifdef DEBUG_ADDER
	fprintf(stderr, "sum result: \t");
	printApUint(sum);
#endif

	ap_uint<(1<<LOG2_EXT_SUM_SIZE)> extSum = ((ap_uint<(1<<LOG2_EXT_SUM_SIZE)>) sum) << ((1<<LOG2_EXT_SUM_SIZE) - (S_WF+1 + S_WF +1)-1);

#ifdef DEBUG_ADDER
	fprintf(stderr, "extsum: \t");
	printApUint(extSum);
#endif

	ap_uint<(LOG2_EXT_SUM_SIZE + (1<<LOG2_EXT_SUM_SIZE))> lzocShifter = lzoc_shifter<LOG2_EXT_SUM_SIZE>(extSum, extSum[EXT_SUM_SIZE -1]);

	ap_uint<LOG2_EXT_SUM_SIZE> lzoc = lzocShifter.range(LOG2_EXT_SUM_SIZE + EXT_SUM_SIZE-1,EXT_SUM_SIZE);
	ap_uint<EXT_SUM_SIZE> shiftedSum = lzocShifter.range(EXT_SUM_SIZE-1,0);

#ifdef DEBUG_ADDER
	fprintf(stderr, "lzoc: \t %d\n", (int) lzoc);
	fprintf(stderr, "shifted sum: \t");
	printApUint(shiftedSum);

	fprintf(stderr, "greatest exp: \t");
	printApUint(subExpOp1);

#endif

#ifdef DEBUG_ADDER
	fprintf(stderr, "lzoc: \t %d\n", (int) lzoc);
	fprintf(stderr, "shifted sum: \t");
	printApUint(shiftedSum);

	fprintf(stderr, "greatest exp: \t");
	printApUint(subExpOp1);

#endif

	// We add two bits to check for both overflows and negative exponents
	ap_uint<S_WE +1 +1> computedExp = subExpOp1 +1 - (lzoc-1);
	ap_uint<1> expIsNegative = computedExp[S_WE +1 +1 -1];
	ap_uint<1> expOverflowed = computedExp[S_WE +1 +1 -1 -1] == 1;
	ap_uint<S_WE> finalExp =computedExp.range(S_WE-1, 0) - PositDim<N>::EXP_BIAS;

#ifdef DEBUG_ADDER
	fprintf(stderr, "final exp: \t");
	printApUint(finalExp);
#endif

	ap_uint<S_WES> es = finalExp.range(S_WES-1,0);
	ap_int<K_SIZE> k = finalExp.range(S_WE-1,S_WES);

	ap_uint<EXT_SUM_SIZE-1> shiftedSumWoImplicitBit = shiftedSum.range(EXT_SUM_SIZE-1 -1, 0);
	ap_uint<S_WES+EXT_SUM_SIZE-1> esAndExactSignificand = es.concat(shiftedSumWoImplicitBit);

	ap_uint<2> zero_one = 0b01;
	ap_uint<2> one_zero = 0b10;

	ap_int<2+S_WES+EXT_SUM_SIZE-1> reverseBitAndEsAndExactSignificand;

	if(k[K_SIZE-1] == 1){
		reverseBitAndEsAndExactSignificand = zero_one.concat(esAndExactSignificand);
	}
	else{
		reverseBitAndEsAndExactSignificand = one_zero.concat(esAndExactSignificand);
	}

#ifdef DEBUG_ADDER
	fprintf(stderr, "exact posit before shift : \n");
	printApInt(reverseBitAndEsAndExactSignificand);
#endif

	ap_uint<K_SIZE-1> abs_k;
	
	if(k[K_SIZE-1] == 1){
		abs_k = ~k;
	}
	else{
		abs_k = k;
	}

#ifdef DEBUG_ADDER
	fprintf(stderr, "shift value when reconstructing the posit: \t");
	printApUint(abs_k);
#endif

	ap_int<2+S_WES+EXT_SUM_SIZE+K_SIZE-1> shiftedReverseBitAndEsAndExactSignificand = (((ap_int<2+S_WES+EXT_SUM_SIZE+K_SIZE-1>) reverseBitAndEsAndExactSignificand) << K_SIZE-1) >> abs_k;

#ifdef DEBUG_ADDER
	fprintf(stderr, "result as posit with infinite accuracy: \n");
	printApInt(shiftedReverseBitAndEsAndExactSignificand);
#endif

	ap_uint<N-1> unroundedPositWoSign = shiftedReverseBitAndEsAndExactSignificand.range(2+S_WES+EXT_SUM_SIZE+K_SIZE-1-1, 2+S_WES+EXT_SUM_SIZE+K_SIZE-1-1 -(N-1) +1);
	ap_uint<2+S_WES+EXT_SUM_SIZE+K_SIZE-1 - (N-1)> remainingBits = shiftedReverseBitAndEsAndExactSignificand.range(2+S_WES+EXT_SUM_SIZE+K_SIZE-1-1 -(N-1) +1 +1, 0);

#ifdef DEBUG_ADDER
	fprintf(stderr, "result as posit N without sign: \n");
	printApUint(unroundedPositWoSign);
	fprintf(stderr, "remaining bits: \n");
	printApUint(remainingBits);
#endif

	ap_uint<1> resultSign = 0;

	ap_uint<N> result = resultSign.concat(unroundedPositWoSign);

	ap_uint<S_WF+1> resultSignificand = shiftedSum.range(EXT_SUM_SIZE-1,EXT_SUM_SIZE -S_WF+1);
	ap_uint<EXT_SUM_SIZE -S_WF+1> resultRest = shiftedSum.range(EXT_SUM_SIZE -S_WF+1 -1,0);

	ap_uint<1> firstRestBit = resultRest[EXT_SUM_SIZE -S_WF+1 -1];
	ap_uint<1> remainingRestBitsAreZeros = resultRest.range(EXT_SUM_SIZE -S_WF+1 -1 -1, 0) == 0;

	ap_uint<1> roudingBit = (firstRestBit && !remainingRestBitsAreZeros) || (firstRestBit && remainingRestBitsAreZeros && (resultSignificand[0]==1));
	// What if this sum overflows?
	ap_uint<S_WF+1> roundedResultSignificand = resultSignificand + roudingBit;


	ap_uint<1> resultIsNaR = in1.getIsNaR() || in1.getIsNaR();
	ap_uint<1> resultS =  (roundedResultSignificand == 0) ? 1 : (!roundedResultSignificand[S_WF+1 -1]); 

	ap_uint<S_WE> resultExp = (roundedResultSignificand == 0) ? 1 : computedExp.range(S_WE-1,0);

	if (expIsNegative){
		// return minpos;
	}
	else if (expOverflowed){
		// return maxpos;
	}
	else{
		// return result;		
	}

	return result;


}

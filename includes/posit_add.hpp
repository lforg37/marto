#pragma once
#include <cstdio>

#include "ap_int.h"
#include "posit_dim.hpp"
#include "lzoc_shifter.hpp"
#include "utils.hpp"


#define S_WF PositIntermediateFormat<N, WES>::FractionSize
#define S_WE PositIntermediateFormat<N, WES>::ExpSize
#define S_WES WES
#define K_SIZE (S_WE-S_WES)

/*
	Uppon testing, the cheapest way to perform the add_sub component is to 
	negate the input before the addition instead of merging the negation in
	the operator.
*/

template<int N, int WES>
PositIntermediateFormat<N, WES> posit_add(
        PositIntermediateFormat<N, WES> in1,
        PositIntermediateFormat<N, WES> in2
){
	#pragma HLS INLINE
	static constexpr int EXT_SUM_SIZE = Static_Val<S_WF+2 + S_WF +1>::_2pow;
	static constexpr int LOG2_EXT_SUM_SIZE = Static_Val<EXT_SUM_SIZE>::_log2;
	
	bool in1IsGreater = in1.getExp() > in2.getExp();

	ap_uint<S_WE> subExpOp1, subExpOp2;
	ap_uint<S_WE+1> shiftValue;
	ap_int<S_WF+2> mostSignificantSignificand, lessSignificantSignificand;

	ap_int<S_WF+2> input2Significand = in2.getSignedSignificand();
	ap_int<S_WF+2> input1Significand = in1.getSignedSignificand();

#ifdef DEBUG_ADDER
	fprintf(stderr, "Inputs 1,2\n");
	printApInt(input1Significand);
	printApInt(input2Significand);
#endif

	if(in1IsGreater){
		subExpOp1 = in1.getExp();
		subExpOp2 = in2.getExp();
		mostSignificantSignificand = input1Significand;
		lessSignificantSignificand = input2Significand;
	}
	else{
		subExpOp1 = in2.getExp();
		subExpOp2 = in1.getExp();
		mostSignificantSignificand = input2Significand;
		lessSignificantSignificand = input1Significand;	
	}

#ifdef DEBUG_ADDER
	fprintf(stderr, "Most, less\n");
	printApInt(mostSignificantSignificand);
	printApInt(lessSignificantSignificand);
#endif

	shiftValue = subExpOp1 - subExpOp2;

#ifdef DEBUG_ADDER
	fprintf(stderr, "Shift value\n");
	printApUint(shiftValue);
#endif

	ap_uint<S_WF+1> toConcatLess = ap_uint<S_WF+1>(0);
	ap_int<S_WF+2 + S_WF+1> shiftedSignificand = ((ap_int<S_WF+2 + S_WF+1>)lessSignificantSignificand.concat(toConcatLess)) >> shiftValue;

#ifdef DEBUG_ADDER
	fprintf(stderr, "Shifted\n");
	printApInt(shiftedSignificand);
#endif

	ap_uint<S_WF+3> shifted_top = shiftedSignificand.range(S_WF+2 + S_WF+1 -1, S_WF+1-1);
	ap_uint<S_WF> sticky_bits = shiftedSignificand.range(S_WF-1, 0);

#ifdef DEBUG_ADDER
	fprintf(stderr, "Shifted top, sticky bits\n");
	printApUint(shifted_top);
	printApUint(sticky_bits);
#endif

	ap_uint<1> shifted_guard = sticky_bits[S_WF-1];
	ap_uint<S_WF-1> rest_sticky_shift = sticky_bits.range(S_WF-2,0);

	ap_uint<1> sticky = rest_sticky_shift.or_reduce();
#ifdef DEBUG_ADDER
	fprintf(stderr, "sticky_full\n");
	printApUint(sticky_full);
#endif

	ap_uint<2> guard_sticky_shifted = shifted_guard.concat(ap_uint<1>(sticky));

#ifdef DEBUG_ADDER
	fprintf(stderr, "guard_sticky_shifted\n");
	printApUint(guard_sticky_shifted);
#endif

	ap_int<S_WF+4> sum_op_1 = (ap_int<S_WF+3>)mostSignificantSignificand.concat(ap_uint<1>(0));
	ap_int<S_WF+4> sum_op_2 = (ap_int<S_WF+3>)shifted_top;

#ifdef DEBUG_ADDER
	fprintf(stderr, "Sum op1, op2, lastbit\n");
	printApInt(sum_op_1);
	printApInt(sum_op_2);
	printApUint(carry_bit);
#endif


	ap_uint<S_WF+4> sum = sum_op_1 + sum_op_2 ;
	ap_uint<1> sum_sign = sum[S_WF+4-1];

#ifdef DEBUG_ADDER
	fprintf(stderr, "Sum\n");
	printApUint(sum);
#endif	

	ap_uint<2> to_append = guard_sticky_shifted;

	ap_uint<S_WF+6> sum_ext = sum.concat(to_append);
	ap_uint<Static_Val<S_WF+6>::_rlog2 + S_WF+6> lzocShifter = generic_lzoc_shifter<S_WF+6>(sum_ext, sum_sign);

	ap_uint<Static_Val<S_WF+6>::_rlog2> lzoc = lzocShifter.range(Static_Val<S_WF+6>::_rlog2 + S_WF+6-1,S_WF+6);
	ap_uint<S_WF+6> shiftedSum = lzocShifter.range(S_WF+6-1,0);

#ifdef DEBUG_ADDER
	fprintf(stderr, "Lzoc, shift\n");
	printApUint(lzoc);
	printApUint(shiftedSum);
#endif

	ap_uint<S_WE +1 +1> computedExp = subExpOp1+1 - (lzoc-1);
	// ap_uint<1> expIsNegative = computedExp[S_WE +1 +1 -1];
	// ap_uint<1> expOverflowed = computedExp[S_WE +1 +1 -1 -1] == 1;


	ap_uint<S_WF+1> resultSignificand = shiftedSum.range(S_WF+6-1,S_WF+6-1 -(S_WF+1)+1);
	ap_uint<S_WF+6 -(S_WF+1)> resultRest = shiftedSum.range(S_WF+6 -(S_WF+1)-1,0);
#ifdef DEBUG_ADDER
	fprintf(stderr, "Significand, rest\n");
	printApUint(resultSignificand);
	printApUint(resultRest);
#endif


	ap_uint<1> guardBit = resultRest[S_WF+6 -(S_WF+1)-1];
	ap_uint<S_WF+6 -(S_WF+1)-1> rest_sticky = resultRest.range(S_WF+6 -(S_WF+1) -1-1, 0);
	ap_uint<1> stickyBit = rest_sticky.or_reduce();

#ifdef DEBUG_ADDER
	fprintf(stderr, "Guard, sticky\n");
	printApUint(guardBit);
	printApUint(stickyBit);
#endif

	ap_uint<1> resultIsNaR = in1.getIsNaR() || in2.getIsNaR();
	ap_uint<1> isZero = ((resultSignificand == 0) && ((guardBit == 0) || ((guardBit == 1) && (stickyBit == 0)))) && !(sum_sign);
	ap_uint<1> resultS =  (isZero) ? 0 : (!resultSignificand[S_WF+1 -1]); 

	ap_uint<S_WE> resultExp = (isZero) ? 0 : computedExp.range(S_WE-1,0);

    PositIntermediateFormat<N, WES> result {
				guardBit,
				stickyBit,
				resultIsNaR,
				resultExp,
				resultS,
				resultSignificand[S_WF+1 -1],
                resultSignificand.range(S_WF+1 -1 -1, 0)};


	return result;
}

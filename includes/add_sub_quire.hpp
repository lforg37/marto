#pragma once
#include "ap_int.h"
#include "shifter.hpp"
#include "posit_dim.hpp"


template<int N>
Quire<N> add_sub_quire(
		Quire<N> quire, 
		PositProd<N> input,
		ap_uint<1> isSub
){
	static constexpr int LOG2_EXT_SUM_SIZE = Static_Val<PositDim<N>::ExtQuireSize-1 +PositDim<N>::ProdSignificandSize>::_log2;

	ap_int<PositDim<N>::ProdSignificandSize+1> inputSignificand = input.getSignedSignificand();
	ap_int<PositDim<N>::ProdSignificandSize+1> complementedInputIfIsSub;
	ap_uint<1> sign = inputSignificand[PositDim<N>::ProdSignificandSize];
//	printApUint(sign);
	#pragma HLS UNROLL
	for (int i=0; i<PositDim<N>::ProdSignificandSize+1; i++){
		complementedInputIfIsSub[i] = input[i] ^ isSub;
	}
	
/*	fprintf(stderr, "=== original ===\n");
	printApInt(inputSignificand);
	fprintf(stderr, "=== complemented ===\n");
	printApInt(complementedInputIfIsSub);
*/


	ap_uint<PositDim<N>::ProdExpSize> shiftValue = input.getExp();

	// fprintf(stderr, "=== shiftValue ===\n");
	// printApUint(shiftValue);

	ap_int<1<<LOG2_EXT_SUM_SIZE> ext = complementedInputIfIsSub;
	//ap_int<PositDim<N>::ExtQuireSize-1 + PositDim<N>::ProdSignificandSize> shiftedInput = 
	//(ap_int<PositDim<N>::ExtQuireSize-1 +PositDim<N>::ProdSignificandSize>)complementedInputIfIsSub<<(shiftValue);
	ap_int<(1<<LOG2_EXT_SUM_SIZE)> shiftedInput = shifter<LOG2_EXT_SUM_SIZE>(ext,shiftValue,isSub);
	ap_int<PositDim<N>::ExtQuireSize-1> shiftedInputShrinked = shiftedInput.range(PositDim<N>::ExtQuireSize-1 + PositDim<N>::ProdSignificandSize-1, PositDim<N>::ProdSignificandSize);
/*	fprintf(stderr, "=== shiftdValue ===\n");
	printApUint(shiftValue);

	fprintf(stderr, "=== shiftedInput ===\n");
	printApInt(shiftedInput);
*/
	ap_uint<PositDim<N>::ExtQuireSize-1> quireWithoutSignAndNARBit = quire.getQuireWithoutNaR();

	ap_uint<PositDim<N>::ExtQuireSize-1> sumResult = shiftedInputShrinked + quireWithoutSignAndNARBit + isSub;

	ap_uint<1> resultIsNaR = quire.getIsNaR() || input.getIsNaR();

	return Quire<N>(resultIsNaR.concat(sumResult));
}


template<int bankSize> 
constexpr int getIndex(int index, bool isUpper)
{
	return bankSize*index - isUpper;
}



template <int N, int bankSize, int spread>
ap_uint<bankSize> getToAddRec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<PositDim<N>::ProdExpSize - getShiftSize<bankSize>()> stageSelect,
    ap_uint<1> inputSign,
    ap_uint<1> isSub,
    ap_int<getExtShiftSize<N, bankSize>()> shiftedSignificand,
    typename enable_if<(spread < 1)>::type* dummy = 0
)
{


	// -1-1 when banksize is 16
	// -1 when banksize is 32
	
	static constexpr int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositDim<N>::ProdSignificandSize,bankSize>::val;

	if((stageIndex>=(stageSelect+getMantSpread<N, bankSize>()-BANKS_FOR_USELESS_BITS) && (inputSign ^ isSub)) || ((stageIndex<stageSelect) && isSub )){
		return -1;
	}
	else{
		return 0;
	}
}

template <int N, int bankSize, int spread>
ap_uint<bankSize> getToAddRec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<PositDim<N>::ProdExpSize - getShiftSize<bankSize>()> stageSelect,
    ap_uint<1> inputSign,
    ap_uint<1> isSub,
    ap_int<getExtShiftSize<N, bankSize>()> shiftedSignificand,
    typename enable_if<(spread >= 1)>::type* dummy = 0
    )
{	
	// fprintf(stderr, "stageselect: ");
	// printApUint(stageSelect);
	// fprintf(stderr, "spread: %d\n", spread);
	// fprintf(stderr, "if %d == %d(select) \n", (int)stageIndex, (int)(stageSelect+spread-1-1));
	// fprintf(stderr, "\t return ");
	// ap_uint<bankSize> tmp = shiftedSignificand.range((spread*bankSize)-1,(spread-1)*bankSize);
	// printApUint(tmp);
	// fprintf(stderr, "else \n"); 
  
	// fprintf(stderr, "spread: %d\n", spread);

	static constexpr int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositDim<N>::ProdSignificandSize,bankSize>::val;

    if (stageIndex == (stageSelect+spread-1-BANKS_FOR_USELESS_BITS)) {
    	// int start = (spread*bankSize)-1;
    	// int end = (spread-1)*bankSize;
    	// fprintf(stderr, "start: %d, end: %d @ stage %d, returns:\n", start, end, (int) stageIndex);
    	// ap_uint<bankSize> tmp = shiftedSignificand.range((spread*bankSize)-1,(spread-1)*bankSize);
    	// printApUint(tmp);
        return shiftedSignificand.range((spread*bankSize)-1,(spread-1)*bankSize);
    } else {
        return getToAddRec<N, bankSize, (spread-1)>(stageIndex, stageSelect, inputSign, isSub, shiftedSignificand);
    }
}


template <int N, int bankSize, int spread>
ap_uint<bankSize> getToAdd(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<PositDim<N>::ProdExpSize - getShiftSize<bankSize>()> stageSelect,
    ap_uint<1> inputSign,
    ap_uint<1> isSub,
    ap_int<getExtShiftSize<N, bankSize>()> shiftedSignificand
		){
	return getToAddRec<N, bankSize, spread>(stageIndex, stageSelect, inputSign, isSub, shiftedSignificand);
}

template<int N, int bankSize>
ap_uint<bankSize+1> add_sub_quire_stage(SegmentedQuire<N, bankSize> quire, 
										ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
										ap_uint<PositDim<N>::ProdExpSize - getShiftSize<bankSize>()> stageSelect,
										ap_uint<1> inputSign, 
										ap_uint<1> isSub,
										ap_int<getExtShiftSize<N, bankSize>()> shiftedSignificand)
{
	// fprintf(stderr, "=== stage index ===\n");
	// printApUint(stageIndex);

	ap_uint<bankSize+1> quireBank = quire.getBank(stageIndex);
	
	ap_uint<1> quireCarry = (stageIndex==0) ? isSub : quire.getCarry(stageIndex-1);

	ap_uint<bankSize> toAdd = getToAdd<N, bankSize, getMantSpread<N, bankSize>()>(stageIndex, stageSelect, inputSign, isSub, shiftedSignificand);

	// fprintf(stderr, "sum = quirebank + toAdd + quireCarry @ stage %d\n", (int)stageIndex);
	// printApUint(quireBank);
	// printApUint(toAdd);
	// printApUint(quireCarry);

	ap_uint<bankSize+1> sum = quireBank + toAdd + quireCarry;
	return sum;
}


template<int N, int bankSize>
SegmentedQuire<N, bankSize> segmented_add_sub_quire(SegmentedQuire<N, bankSize> quire, 
													PositProd<N> input,
													ap_uint<1> isSub)
{	

	static constexpr int SHIFTED_SIGNIFICAND_SIZE = PositDim<N>::ProdSignificandSize+1 + (1<<getShiftSize<bankSize>());
	static constexpr int BANKS_FOR_USELESS_BITS = Static_Ceil_Div<PositDim<N>::ProdSignificandSize,bankSize>::val;
	static constexpr int ENCODING_BITS_FOR_USELESS_BANKS = get2Power(BANKS_FOR_USELESS_BITS);
	static constexpr int padding = bankSize*BANKS_FOR_USELESS_BITS - PositDim<N>::ProdSignificandSize;
	static constexpr int LOG2_SHIFT_SIZE = Static_Val<getExtShiftSize<N, bankSize>()>::_log2;

	// fprintf(stderr, "padding: %d\n", padding);

	ap_uint<PositDim<N>::ProdExpSize+ENCODING_BITS_FOR_USELESS_BANKS> prodExp = input.getExp()+padding;
	ap_int<PositDim<N>::ProdSignificandSize+1> inputSignificand = input.getSignedSignificand();
	ap_uint<1> sign = inputSignificand[PositDim<N>::ProdSignificandSize];
	
	ap_int<PositDim<N>::ProdSignificandSize+1> complementedInputIfIsSub;
	#pragma HLS UNROLL
	for (int i=0; i<PositDim<N>::ProdSignificandSize+1; i++){
		complementedInputIfIsSub[i] = input[i] ^ isSub;
	}

	// fprintf(stderr, "=== complemented if sub ===\n");
	// printApInt(complementedInputIfIsSub);

	ap_uint<getShiftSize<bankSize>()> shiftValue = prodExp.range(getShiftSize<bankSize>()-1,0);

	ap_int<1<<LOG2_SHIFT_SIZE> ext = complementedInputIfIsSub;
	ap_int<(1<<LOG2_SHIFT_SIZE)> shiftedInput = shifter<LOG2_SHIFT_SIZE>(ext,shiftValue,isSub);
	ap_int<getExtShiftSize<N, bankSize>()> shiftedInputShrinked = shiftedInput.range(getExtShiftSize<N, bankSize>()-1, 0);


	// ap_int<getExtShiftSize<N, bankSize>()> shiftedSignificand = (ap_int<getExtShiftSize<N, bankSize>()>) complementedInputIfIsSub << shiftValue;
	// fprintf(stderr, "=== shifted significand ===\n");
	// printApInt(shiftedSignificand);
	// fprintf(stderr, "=== shiftedInputShrinked ===\n");
	// printApInt(shiftedInputShrinked);
	// fprintf(stderr, "Significand size: %d\n", getExtShiftSize<N, bankSize>());

	ap_uint<PositDim<N>::ProdExpSize - getShiftSize<bankSize>()> stageSelect = prodExp.range(PositDim<N>::ProdExpSize-1, getShiftSize<bankSize>());
	// fprintf(stderr, "=== prod exp ===\n");
	// printApUint(prodExp);	

	// fprintf(stderr, "=== stage select ===\n");
	// printApUint(stageSelect);
	// fprintf(stderr, "=== shift value ===\n");
	// printApUint(shiftValue);

	SegmentedQuire<N, bankSize> fullQuire = SegmentedQuire<N, bankSize>(0);
	// ap_uint<1> carry = 0;
	// int tmp =getNbStages<N, bankSize>();
	// fprintf(stderr, "nbstages: %d \n", tmp);
	#pragma HLS UNROLL
	for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
		ap_uint<bankSize+1> stageResult = add_sub_quire_stage<N,bankSize>(quire, i, stageSelect, sign, isSub, shiftedInputShrinked);
		fullQuire[i] = stageResult[bankSize];
		// fprintf(stderr, "I : %d, result:\n", i);
		// printApUint(stageResult);
		// fprintf(stderr, "carry = %d\n",(int) fullQuire[i]);
		// fprintf(stderr, "oldcarry = %d\n",(int) quire[i]);
		fullQuire.range(getIndex<bankSize>(i+1, 1)+getNbStages<N, bankSize>(), getIndex<bankSize>(i, 0)+getNbStages<N, bankSize>()) = stageResult.range(bankSize-1,0);
	}
	ap_uint<1> resultIsNaR = input.getIsNaR() || quire.getIsNaR();
	ap_uint<PositDim<N>::ExtQuireSize+getNbStages<N, bankSize>()-1> fullQuireWoNaR = fullQuire.range(PositDim<N>::ExtQuireSize+getNbStages<N, bankSize>()-1-1, 0);
	return SegmentedQuire<N, bankSize>(resultIsNaR.concat(fullQuireWoNaR));
}

template<int N, int bankSize>
Quire<N> propagateCarries(SegmentedQuire<N, bankSize> quire)
{	
	// fprintf(stderr, "=== input quire ===\n");
	// quire.printContent();

	SegmentedQuire<N, bankSize> fullQuire = quire;
	for(int j=0; j<getNbStages<N, bankSize>(); j++){
		#pragma HLS UNROLL
		for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
			ap_uint<bankSize+1> stageResult = add_sub_quire_stage<N,bankSize>(fullQuire, i, 0, 0, 0, 0);
			fullQuire[i] = stageResult[bankSize];
			// printApUint(stageResult);
			// int start = getIndex<bankSize>(i+1, 1)+getNbStages<N, bankSize>();
			// int end = getIndex<bankSize>(i, 0)+getNbStages<N, bankSize>();
			// fprintf(stderr, "START: %d, END:%d, length: %d\n", start, end, start-end+1);
			fullQuire.range(getIndex<bankSize>(i+1, 1)+getNbStages<N, bankSize>(), getIndex<bankSize>(i, 0)+getNbStages<N, bankSize>()) = stageResult.range(bankSize-1,0);
		}
		// fprintf(stderr, "=== quire at step %d ===\n", j);
		// fullQuire.printContent();
		
	}
	// fprintf(stderr, "=== Before returning quire ===\n");
	// fullQuire.printContent();
	return fullQuire.getAsQuireWoCarries();
}
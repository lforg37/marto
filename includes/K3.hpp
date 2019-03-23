
#include "kulisch_dim.hpp"
using namespace std;
#include "shifter.hpp"

// template<int N, int bankSize>
// static constexpr int ext_shift_size(){
// 	return Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val * bankSize;
// }

template<int bankSize> 
constexpr int acc_2CK3_getIndex(int index, bool isUpper)
{
	return bankSize*index - isUpper;
}

template <int N, int bankSize, int spread>
ap_uint<bankSize> acc_2CK3_to_add_Rec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)> shiftedSignificand,
    typename enable_if<(spread < 1)>::type* dummy = 0
)
{
	#pragma HLS INLINE
	
	if((stageIndex>=(stageSelect+getMantSpread<N, bankSize>()) && inputSign)){
		return -1;
	}
	else{
		return 0;
	}
}

template <int N, int bankSize, int spread>
ap_uint<bankSize> acc_2CK3_to_add_Rec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)> shiftedSignificand,
    typename enable_if<(spread >= 1)>::type* dummy = 0
    )
{	
	#pragma HLS INLINE
    if (stageIndex == (stageSelect+spread-1)) {
        ap_uint<bankSize> tmp =shiftedSignificand.range(((spread)*bankSize)-1,(spread-1)*bankSize); 
        return tmp;
    } else {
        return acc_2CK3_to_add_Rec<N, bankSize, (spread-1)>(stageIndex, stageSelect, inputSign, shiftedSignificand);
    }
}


template <int N, int bankSize, int spread>
ap_uint<bankSize> acc_2CK3_to_add_(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)	> shiftedSignificand
		){
	#pragma HLS INLINE
	return acc_2CK3_to_add_Rec<N, bankSize, spread>(stageIndex, stageSelect, inputSign, shiftedSignificand);
}

template<int N, int bankSize>
ap_uint<bankSize+1> acc_2CK3_acc_stage(SegmentedKulischAcc<N, bankSize> acc, 
										ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
										ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
										ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)> shiftedSignificand,
										ap_uint<1> sign,
										ap_uint<1> carry0 =0)
{
	#pragma HLS INLINE
	ap_uint<bankSize+1> bank = acc.getBank(stageIndex);
	// fprintf(stderr, "main input with size %d: ", (int) (1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2));
	// printApInt(shiftedSignificand);
	
	ap_uint<1> accCarry = (stageIndex==0) ? carry0 : acc.getCarry(stageIndex-1);

	ap_uint<bankSize> toAdd = acc_2CK3_to_add_<N, bankSize, getMantSpread<N, bankSize>()>(stageIndex, stageSelect, sign, shiftedSignificand);

	ap_uint<bankSize+1> sum = bank + toAdd + accCarry;
	return sum;
}




template<int N, int bankSize>
SegmentedKulischAcc<N, bankSize> acc_2CK3(
													SegmentedKulischAcc<N, bankSize> acc, 
													FPProd<N> prod
													)
{	
	#pragma HLS INLINE

	static constexpr int spread = Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val+1;
	static constexpr int bits_for_shift = Static_Val<spread*bankSize>::_log2;

	ap_uint<FPDim<N>::WE+1 +1> prodExp = (ap_uint<FPDim<N>::WE+1 +1>)(prod.getExp());

	ap_uint<2*FPDim<N>::WF+2> inputSignificand = prod.getSignificand();
	ap_uint<2*FPDim<N>::WF+3> inputSignificandComplemented;
	if(prod.getSignBit() == 1){
		inputSignificandComplemented = ~((ap_uint<2*FPDim<N>::WF+3>)inputSignificand) + 1;
	}
	else{
		inputSignificandComplemented = inputSignificand;
	}

	ap_int<2*FPDim<N>::WF+2+2> inputSignificandWithSign = prod.getSignBit().concat(inputSignificandComplemented);

	Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val *bankSize ;


	ap_uint<Static_Val<bankSize>::_log2> shiftValue = prodExp.range(Static_Val<bankSize>::_log2-1,0);
	
	ap_int<(1<<bits_for_shift)> ext = inputSignificandWithSign;

	ap_int<(1<<bits_for_shift)> shiftedInput = shifter<bits_for_shift>(ext,shiftValue,0);

	ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect = prodExp.range(FPDim<N>::WE+1 +1 -1, Static_Val<bankSize>::_log2);
	SegmentedKulischAcc<N, bankSize> fullAcc = SegmentedKulischAcc<N, bankSize>(0, 0);

	for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
		#pragma HLS UNROLL
		ap_uint<bankSize+1> stageResult = acc_2CK3_acc_stage<N,bankSize>(acc, i, stageSelect, shiftedInput, prod.getSignBit());
		fullAcc[i] = stageResult[bankSize];
		fullAcc.range(acc_2CK3_getIndex<bankSize>(i+1, 1)+getNbStages<N, bankSize>(), acc_2CK3_getIndex<bankSize>(i, 0)+getNbStages<N, bankSize>()) = stageResult.range(bankSize-1,0);
	}

	return fullAcc;
}

template<int N, int bankSize>
KulischAcc<N> acc_2CK3_propagate_carries(SegmentedKulischAcc<N, bankSize> acc)
{	
	#pragma HLS INLINE
	SegmentedKulischAcc<N, bankSize> fullacc = acc;
	for(int j=0; j<getNbStages<N, bankSize>()+1; j++){
		for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
			#pragma HLS UNROLL
			ap_uint<bankSize+1> stageResult = acc_2CK3_acc_stage<N, bankSize>(fullacc, (ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> ) i, (ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1>) 0, (ap_uint<Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val * bankSize>) 0, 0);
			// ap_uint<bankSize+1> stageResult = acc_2CK3_acc_stage<N, bankSize>(fullacc, (ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> ) i, (ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1>) 0, (ap_uint<Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val * bankSize>) 0, (ap_uint<1>)((j==0) and (i==0) and acc.isNeg()));
			fullacc[i] = stageResult[bankSize];
			fullacc.range(acc_2CK3_getIndex<bankSize>(i+1, 1)+getNbStages<N, bankSize>(), acc_2CK3_getIndex<bankSize>(i, 0)+getNbStages<N, bankSize>()) = stageResult.range(bankSize-1,0);
		}	
	}
	// fullacc.printContent();
	return fullacc.range(FPDim<N>::ACC_SIZE + getNbStages<N, bankSize>()-1, getNbStages<N, bankSize>());
}
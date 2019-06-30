
#include "kulisch_dim.hpp"
using namespace std;
#include "marto/bitvector.hpp"

using hint::Static_Val;
using hint::Static_Ceil_Div;

template <int N, int bankSize, int spread>
ap_uint<bankSize> add_2CK3_to_add_Rec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
	ap_int<(1<<Static_Val<spread*bankSize>::_log2)>,
	typename enable_if<(spread < 1)>::type* = 0
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
ap_uint<bankSize> add_2CK3_to_add_Rec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)> shiftedSignificand,
	typename enable_if<(spread >= 1)>::type* = 0
    )
{	
	#pragma HLS INLINE
    if (stageIndex == (stageSelect+spread-1)) {
        ap_uint<bankSize> tmp = shiftedSignificand.range(((spread)*bankSize)-1,(spread-1)*bankSize); 
        return tmp;
    } else {
        return add_2CK3_to_add_Rec<N, bankSize, (spread-1)>(stageIndex, stageSelect, inputSign, shiftedSignificand);
    }
}


template <int N, int bankSize, int spread>
ap_uint<bankSize> add_2CK3_to_add_(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)	> shiftedSignificand
		){
	#pragma HLS INLINE
	return add_2CK3_to_add_Rec<N, bankSize, spread>(stageIndex, stageSelect, inputSign, shiftedSignificand);
}

template<int N, int bankSize>
ap_uint<bankSize+1> add_2CK3_acc_stage(acc_2CK3<N, bankSize> acc, 
										ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
										ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
										ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)> shiftedSignificand,
										ap_uint<1> sign,
										ap_uint<1> carry0 =0)
{
	#pragma HLS INLINE
	ap_uint<bankSize+1> bank = acc.getBank(stageIndex);
	ap_uint<1> accCarry = (stageIndex==0) ? carry0 : acc.getCarry(stageIndex-1);

	ap_uint<bankSize> toAdd = add_2CK3_to_add_<N, bankSize, getMantSpread<N, bankSize>()>(stageIndex, stageSelect, sign, shiftedSignificand);

	ap_uint<bankSize+1> op1 = bank;
	ap_uint<bankSize+1> op2 = toAdd;
	ap_uint<1> c = accCarry;
	ap_uint<bankSize+1> sum = op1 + op2 + 1;
	#pragma HLS RESOURCE variable=sum core=AddSub
	return sum;
}




template<int N, int bankSize>
acc_2CK3<N, bankSize> add_2CK3(
													acc_2CK3<N, bankSize> acc, 
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



	ap_uint<Static_Val<bankSize>::_log2> shiftValue = prodExp.range(Static_Val<bankSize>::_log2-1,0);
	
	ap_int<(1<<bits_for_shift)> ext = inputSignificandWithSign;

	ap_int<(1<<bits_for_shift)> shiftedInput = shifter<bits_for_shift>(ext,shiftValue,0);

	ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect = prodExp.range(FPDim<N>::WE+1 +1 -1, Static_Val<bankSize>::_log2);
	acc_2CK3<N, bankSize> fullAcc = acc_2CK3<N, bankSize>(0, 0);

	for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
		#pragma HLS UNROLL
		ap_uint<bankSize+1> stageResult = add_2CK3_acc_stage<N,bankSize>((acc_2CK3<N, bankSize>)acc, (ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> ) i, (ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1>) stageSelect, (ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)>) shiftedInput, (ap_uint<1>) prod.getSignBit());
		fullAcc.setCarry(i, ap_uint<1>{stageResult[bankSize]});
		fullAcc.setBank(i, stageResult.range(bankSize-1,0));
	}

	return fullAcc;
}

template<int N, int bankSize>
KulischAcc<N> propagate_carries_2CK3(acc_2CK3<N, bankSize> acc)
{	
	#pragma HLS INLINE
	acc_2CK3<N, bankSize> fullAcc = acc;
	for(int j=0; j<getNbStages<N, bankSize>()+1; j++){
		#pragma HLS PIPELINE II=1
		for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
			ap_uint<bankSize+1> stageResult = add_2CK3_acc_stage<N, bankSize>(fullAcc, (ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> ) i, (ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1>) 0, (ap_uint<Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val * bankSize>) 0, 0);
			fullAcc.setCarry(i, stageResult[bankSize]);
			fullAcc.setBank(i, stageResult.range(bankSize-1,0));
		}	
	}
	// fullacc.printContent();
	return fullAcc.getAcc();
}















// template<int bankSize> 
// constexpr int add_SMK3_getIndex(int index, bool isUpper)
// {
// 	return bankSize*index - isUpper;
// }


template <int N, int bankSize, int spread>
ap_uint<bankSize> add_SMK3_to_add_Rec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)> shiftedSignificand,
    typename enable_if<(spread < 1)>::type* dummy = 0
)
{
	#pragma HLS INLINE
	return 0;
}

template <int N, int bankSize, int spread>
ap_uint<bankSize> add_SMK3_to_add_Rec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)> shiftedSignificand,
    typename enable_if<(spread >= 1)>::type* dummy = 0
    )
{	
	#pragma HLS INLINE
	ap_uint<bankSize> toAdd;
    if (stageIndex == (stageSelect+spread-1)) {
    	if(inputSign==1){
    		toAdd=0;
    	}
    	else{
    		toAdd=shiftedSignificand.range(((spread)*bankSize)-1,(spread-1)*bankSize);
    	}
        return toAdd;
    } else {
        return add_SMK3_to_add_Rec<N, bankSize, (spread-1)>(stageIndex, stageSelect, inputSign, shiftedSignificand);
    }
}



template <int N, int bankSize, int spread>
ap_uint<bankSize> add_SMK3_to_sub_Rec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)> shiftedSignificand,
    typename enable_if<(spread < 1)>::type* dummy = 0
)
{
	#pragma HLS INLINE
	return 0;
}


template <int N, int bankSize, int spread>
ap_uint<bankSize> add_SMK3_to_sub_Rec(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)> shiftedSignificand,
    typename enable_if<(spread >= 1)>::type* dummy = 0
    )
{	
	#pragma HLS INLINE
	ap_uint<bankSize> toSub;
    if (stageIndex == (stageSelect+spread-1)) {
    	if(inputSign==1){
    		toSub=shiftedSignificand.range(((spread)*bankSize)-1,(spread-1)*bankSize); 
    	}
    	else{
    		toSub=0;
    	}
        return toSub;
    } else {
        return add_SMK3_to_sub_Rec<N, bankSize, (spread-1)>(stageIndex, stageSelect, inputSign, shiftedSignificand);
    }
}



template <int N, int bankSize, int spread>
ap_uint<bankSize> add_SMK3_to_add_(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)	> shiftedSignificand
		){
	#pragma HLS INLINE
	return add_SMK3_to_add_Rec<N, bankSize, spread>(stageIndex, stageSelect, inputSign, shiftedSignificand);
}




template <int N, int bankSize, int spread>
ap_uint<bankSize> add_SMK3_to_sub_(
    ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
    ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
    ap_uint<1> inputSign,
    ap_int<(1<<Static_Val<spread*bankSize>::_log2)	> shiftedSignificand
		){
	#pragma HLS INLINE
	return add_SMK3_to_sub_Rec<N, bankSize, spread>(stageIndex, stageSelect, inputSign, shiftedSignificand);
}



template<int N, int bankSize>
ap_uint<bankSize+1> add_SMK3_acc_sub(acc_SMK3<N, bankSize> acc, 
										ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
										ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
										ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)> shiftedSignificand,
										ap_uint<1> sign,
										ap_uint<1> borrow0 =0)
{
	#pragma HLS INLINE
	ap_uint<bankSize+1+1> bank = acc.getBank(stageIndex);
	ap_uint<bankSize> toSub = add_SMK3_to_sub_<N, bankSize, getMantSpread<N, bankSize>()>(stageIndex, stageSelect, sign, shiftedSignificand);

	ap_uint<bankSize+1> sub_op1{bank};
	ap_uint<bankSize+1> sub_op2{toSub};
	ap_uint<1> accBorrow = (stageIndex==0) ? borrow0 : acc.getBorrow(stageIndex-1);
	ap_uint<bankSize+1>  sub = sub_op1 - sub_op2 - accBorrow;

	return sub;	
}

template<int N, int bankSize>
ap_uint<bankSize+1> add_SMK3_acc_add(acc_SMK3<N, bankSize> acc, 
										ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
										ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
										ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)> shiftedSignificand,
										ap_uint<1> sign,
										ap_uint<bankSize> subed_to_add,
										ap_uint<1> carry0 =0)
{
	#pragma HLS INLINE
	ap_uint<bankSize> toAdd = add_SMK3_to_add_<N, bankSize, getMantSpread<N, bankSize>()>(stageIndex, stageSelect, sign, shiftedSignificand);
	ap_uint<1> accCarry = (stageIndex==0) ? carry0 : acc.getCarry(stageIndex-1);

	ap_uint<bankSize+1> sum_op1{subed_to_add};
	ap_uint<bankSize+1> sum_op2{toAdd};

	ap_uint<bankSize+1> sum = sum_op1 + sum_op2 + accCarry;

	return sum;	
}

template<int N, int bankSize>
ap_uint<bankSize+2> add_SMK3_acc_stage(acc_SMK3<N, bankSize> acc, 
										ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> stageIndex, 
										ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
										ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)> shiftedSignificand,
										ap_uint<1> sign,
										ap_uint<1> carry0 =0,
										ap_uint<1> borrow0 =0)
{
	// #pragma HLS INLINE
	// #pragma HLS PIPELINE
	
	ap_uint<bankSize+1> subed = add_SMK3_acc_sub(acc, stageIndex, stageSelect, shiftedSignificand, sign, borrow0);
	ap_uint<1> borrow = subed[bankSize];
	ap_uint<bankSize> subed_to_add = subed.range(bankSize-1,0);

	ap_uint<bankSize+1> added = add_SMK3_acc_add(acc, stageIndex, stageSelect, shiftedSignificand, sign, subed_to_add, carry0);
	ap_uint<1> carry = added[bankSize];
	ap_uint<bankSize> res = added.range(bankSize-1,0);
	

	ap_uint<1> finalCarry{carry and not(borrow)};
	ap_uint<bankSize+1> resWithCarry{finalCarry.concat(res)};
	ap_uint<1> finalBorrow{borrow and not(carry)};
	ap_uint<bankSize+2> resWithCarryBorrow{borrow.concat(resWithCarry)};

	return resWithCarryBorrow;
}


template<int N, int bankSize, int index>
void add_SMK3_acc_stage_top(acc_SMK3<N, bankSize> * acc, 
										ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
										ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)> shiftedSignificand,
										ap_uint<1> sign,
										ap_uint<1> carry0 =0,
										ap_uint<1> borrow0 =0,
    									typename enable_if<(index == 0)>::type* dummy = 0
										)
{
		ap_uint<bankSize+2> stage_val = add_SMK3_acc_stage<N,bankSize>((acc_SMK3<N, bankSize>)(*acc), stageSelect, shiftedSignificand, sign, carry0, borrow0);
		acc->setCarry(index, stage_val[bankSize]);
		acc->setBorrow(index, stage_val[bankSize+1]);
		acc->setBank(index, stage_val.range(bankSize-1,0));
}

template<int N, int bankSize, int index>
void add_SMK3_acc_stage_top(acc_SMK3<N, bankSize>*  acc, 
										ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect,
										ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)> shiftedSignificand,
										ap_uint<1> sign,
										ap_uint<1> carry0 =0,
										ap_uint<1> borrow0 =0,
    									typename enable_if<(index >= 1)>::type* dummy = 0
										)
{
		ap_uint<bankSize+2> stage_val = add_SMK3_acc_stage<N,bankSize>((acc_SMK3<N, bankSize>)(*acc), (ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2>)index, stageSelect, shiftedSignificand, sign, carry0, borrow0);
		acc->setCarry(index, stage_val[bankSize]);
		acc->setBorrow(index, stage_val[bankSize+1]);
		acc->setBank(index, stage_val.range(bankSize-1,0));
		add_SMK3_acc_stage_top<N, bankSize, index-1>(acc, stageSelect, shiftedSignificand, sign, carry0, borrow0);
}




template<int N, int bankSize>
acc_SMK3<N, bankSize> add_SMK3(
													acc_SMK3<N, bankSize> acc, 
													FPProd<N> prod
													)
{	
	#pragma HLS INLINE

	static constexpr int spread = Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val+1;
	static constexpr int bits_for_shift = Static_Val<spread*bankSize>::_log2;

	ap_uint<FPDim<N>::WE+1 +1> prodExp = (ap_uint<FPDim<N>::WE+1 +1>)(prod.getExp());
	ap_uint<2*FPDim<N>::WF+2> inputSignificand = prod.getSignificand();
	ap_uint<Static_Val<bankSize>::_log2> shiftValue = prodExp.range(Static_Val<bankSize>::_log2-1,0);
	ap_uint<(1<<bits_for_shift)> ext = inputSignificand;
	ap_uint<(1<<bits_for_shift)> shiftedInput = shifter<bits_for_shift>(ext,shiftValue,0);

	ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1> stageSelect = prodExp.range(FPDim<N>::WE+1 +1 -1, Static_Val<bankSize>::_log2);
	acc_SMK3<N, bankSize> fullAcc = acc_SMK3<N, bankSize>(0, 0);


	add_SMK3_acc_stage_top<N, bankSize, getNbStages<N, bankSize>()-1>((acc_SMK3<N, bankSize>*)(&acc), (ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1>) stageSelect, (ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)>) shiftedInput, (ap_uint<1>) prod.getSignBit());

	// for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
	// 	#pragma HLS UNROLL
	// 	ap_uint<bankSize+2> stage_val = add_SMK3_acc_stage<N,bankSize>((acc_SMK3<N, bankSize>)acc, (ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2>)i, (ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1>) stageSelect, (ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)>) shiftedInput, (ap_uint<1>) prod.getSignBit());
	// 	fullAcc.setCarry(i, stage_val[bankSize]);
	// 	fullAcc.setBorrow(i, stage_val[bankSize+1]);
	// 	fullAcc.setBank(i, stage_val.range(bankSize-1,0));
	// }

	return acc;
}

template<int N, int bankSize>
KulischAcc<N> propagate_carries_SMK3(acc_SMK3<N, bankSize> acc)
{	

	#pragma HLS INLINE
	// acc_SMK3<N, bankSize> fullAcc = acc;
  	for(int j=0; j<getNbStages<N, bankSize>()+1; j++){
		#pragma HLS PIPELINE II=1
		
		add_SMK3_acc_stage_top<N, bankSize, getNbStages<N, bankSize>()-1>((acc_SMK3<N, bankSize>*)(&acc), (ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1>) 0, (ap_int<(1<<Static_Val<getMantSpread<N, bankSize>()*bankSize>::_log2)>) 0, (ap_uint<1>) 0);

		// for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){

		// 	ap_uint<bankSize+2> stage_val = add_SMK3_acc_stage<N, bankSize>(fullAcc, (ap_uint<Static_Val<getNbStages<N, bankSize>()>::_log2> ) i, (ap_uint<FPDim<N>::WE+1 +1 - Static_Val<bankSize>::_log2 +1>) 0, (ap_uint<Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val * bankSize>) 0, 0);
		// 	fullAcc.setCarry(i, stage_val[bankSize]);
		// 	fullAcc.setBorrow(i, stage_val[bankSize+1]);
		// 	fullAcc.setBank(i, stage_val.range(bankSize-1,0));
		// }	
	}
	// fullAcc.printContent();
	return acc.getAcc();
}

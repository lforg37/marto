#include "kulisch_dim.hpp"

template<int N>
KulischAcc<N> add_2CK1(
		KulischAcc<N> acc, 
		FPProd<N> prod
){
	#pragma HLS INLINE

	ap_int<2*FPDim<N>::WF+3> ext_significand = prod.getSignificand();

	if(prod.getSignBit() == 1){
		ext_significand = ~((ap_int<2*FPDim<N>::WF+3>)prod.getSignificand()) + 1;
	}
	else{
		ext_significand = prod.getSignificand();
	}

	ap_uint<FPDim<N>::ACC_SIZE> shifted = ((ap_uint<FPDim<N>::ACC_SIZE>) ext_significand) << (prod.getExp());
	ap_uint<FPDim<N>::ACC_SIZE> r_acc = shifted + acc;
	return r_acc;
}


template<int N>
SignedKulischAcc<N> add_SMK1(
		SignedKulischAcc<N> acc, 
		FPProd<N> prod
){
	#pragma HLS INLINE

	ap_int<2*FPDim<N>::WF+3> ext_significand = prod.getSignificand();

	ap_uint<FPDim<N>::ACC_SIZE> summand = ((ap_int<FPDim<N>::ACC_SIZE>) ext_significand) << (prod.getExp());

	ap_uint<1> acc_sign = acc[FPDim<N>::ACC_SIZE];

	ap_uint<1> select;
	if((prod.getSignBit() && acc_sign) || ((!prod.getSignBit()) && (!acc_sign))){
		select = 0;
	}
	else{
		select = 1;
	}

	ap_int<FPDim<N>::ACC_SIZE> operand2;

	if(select == 1){
		operand2 = (~acc);
	}
	else{
		operand2 = acc;
	}

	ap_uint<FPDim<N>::ACC_SIZE> sum;
	ap_uint<FPDim<N>::ACC_SIZE> sub;
	sum = summand + operand2 +(ap_int<FPDim<N>::ACC_SIZE>)select;
	sub = acc - summand; 

	ap_uint<1>select2 = select && (!sub[FPDim<N>::ACC_SIZE-1]);

	SignedKulischAcc<N> res;

	if(select2==1){
		res = sub;
	}
	else{
		res = sum;
	}

	if(select && (!select2)){
		res[FPDim<N>::ACC_SIZE] = !acc_sign;
	}
	else{
		res[FPDim<N>::ACC_SIZE] = acc_sign;
	}


	return res;
}








template<int bankSize> 
constexpr int add_segmented_2CK1_getIndex(int index, bool isUpper)
{
	return bankSize*index - isUpper;
}

template<int N, int bankSize>
acc_segmented_2CK1<N, bankSize> add_segmented_2CK1(
													acc_segmented_2CK1<N, bankSize> acc, 
													FPProd<N> prod
													)
{	
	#pragma HLS INLINE

	ap_int<2*FPDim<N>::WF+3> ext_significand = prod.getSignificand();

	if(prod.getSignBit() == 1){
		ext_significand = ~((ap_int<2*FPDim<N>::WF+3>)prod.getSignificand()) + 1;
	}
	else{
		ext_significand = prod.getSignificand();
	}

	ap_uint< getSegmentedAccSize<N, bankSize>() > shifted = ((ap_uint< getSegmentedAccSize<N, bankSize>() >) ext_significand) << (prod.getExp());



	acc_segmented_2CK1<N, bankSize> fullAcc = acc_segmented_2CK1<N, bankSize>(0, 0);

	for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
		#pragma HLS UNROLL
		ap_uint<1> carry = (i==0) ? (ap_uint<1> )0 : acc.getCarry(i-1);
		ap_uint<bankSize+1> stageResult = acc.getBank(i)
										+ shifted.range((i+1)*bankSize-1,i*bankSize)
										+ carry;
		fullAcc.setCarry(i, stageResult[bankSize]);
		fullAcc.setBank(i, stageResult.range(bankSize-1,0));
	}
	return fullAcc;
}

template<int N, int bankSize>
KulischAcc<N> propagate_carries_segmented_2CK1(acc_segmented_2CK1<N, bankSize> acc)
{	
	#pragma HLS INLINE
	acc_segmented_2CK1<N, bankSize> fullAcc = acc;
	for(int j=0; j<getNbStages<N, bankSize>()+1; j++){
		for(int i=getNbStages<N, bankSize>()-1; i>=0; i--){
			#pragma HLS UNROLL
			ap_uint<1> carry = (i==0) ? (ap_uint<1> )0 : fullAcc.getCarry(i-1);
			ap_uint<bankSize+1> stageResult = fullAcc.getBank(i)
											+ carry;
			fullAcc.setCarry(i, stageResult[bankSize]);
			fullAcc.setBank(i, stageResult.range(bankSize-1,0));
		}	
	}
	// fullacc.printContent();
	return fullAcc.getAcc();;
}


#include "kulisch_dim.hpp"

template<int N>
KulischAcc<N> acc_2CK1(
		KulischAcc<N> acc, 
		FPProd<N> prod
){
	#pragma HLS INLINE

	// cout << "Prod :\n";
	// prod.printContent();

	ap_int<2*FPDim<N>::WF+3> ext_significand = prod.getSignificand();

	if(prod.getSignBit() == 1){
		ext_significand = ~((ap_int<2*FPDim<N>::WF+3>)prod.getSignificand()) + 1;
	}
	else{
		ext_significand = prod.getSignificand();
	}

	ap_uint<FPDim<N>::ACC_SIZE> shifted = ((ap_uint<FPDim<N>::ACC_SIZE>) ext_significand) << (prod.getExp());
	ap_uint<FPDim<N>::ACC_SIZE> r_acc = shifted + acc;

	// cout << "Acc " << FPDim<N>::ACC_SIZE << ":\n";
	// printApUint(r_acc);

	return r_acc;
}
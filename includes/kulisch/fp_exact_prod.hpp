#include "kulisch_dim.hpp"
#include "marto/bitvector.hpp"

template<int N>
FPProd<N> exact_prod(ap_uint<N> in1,
		ap_uint<N> in2
){
	#pragma HLS INLINE
	ap_uint<FPDim<N>::WF+1> m1, m2;
	ap_uint<FPDim<N>::WE> e1, e2;
	ap_uint<1> s1, s2;
	ap_uint<2*FPDim<N>::WF+2> mult_m;
	ap_uint<FPDim<N>::WE+1> mult_e;

	m1 = in1.range(FPDim<N>::WF-1, 0);
	e1 = in1.range(FPDim<N>::WF+FPDim<N>::WE-1, FPDim<N>::WF);
	m2 = in2.range(FPDim<N>::WF-1, 0);
	e2 = in2.range(FPDim<N>::WF+FPDim<N>::WE-1, FPDim<N>::WF);
	s1 = in1[FPDim<N>::WF+FPDim<N>::WE];
	s2 = in2[FPDim<N>::WF+FPDim<N>::WE];

	if(e1 == 0){
		m1[FPDim<N>::WF] = 0;
	}
	else{
		m1[FPDim<N>::WF] = 1;
	}

	if(e2 == 0){
		m2[FPDim<N>::WF] = 0;
	}
	else{
		m2[FPDim<N>::WF] = 1;
	}

	ap_uint<1> mult_s = s1 ^ s2;

	mult_e = (ap_uint<FPDim<N>::WE+1>)e1 + (ap_uint<FPDim<N>::WE+1>)e2;
	ap_uint<2*FPDim<N>::WF+2> mul_op1 = m1;
	ap_uint<2*FPDim<N>::WF+2> mul_op2 = m2;
	mult_m = mul_op1*mul_op2;
	// #pragma HLS RESOURCE variable=mult_m core=Mulns latency=8
	return FPProd<N>(	mult_s,
						mult_e,
						mult_m	
					);
}

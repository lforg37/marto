#pragma once
#include "ap_int.h"
#include "posit_dim.tpp"

template<size_t N>
Quire<N> add_sub_quire(
		Quire<N> quire, 
		PositProd<N> input,
	   	ap_uint<1> is_sub
){
	return 0;
}

#pragma once
#include "posit_dim.hpp"
#include "lzoc_shifter.hpp"

template<int N>
PositValue<N> posit_decoder(PositEncoding<N> positN)
{
	//Sign bit
	ap_uint<1> s = positN[N-1];
	//First regime bit
	ap_uint<1> count_type = positN[N-2];
	//Remainder
	ap_uint<N-2> remainder = positN.range(N-3, 0);

	

	return PositValue<N>(0);
}

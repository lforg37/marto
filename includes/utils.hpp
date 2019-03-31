#pragma once
#include "ap_int.h"
#include <stdio.h>
#include <type_traits>

#include <cstdio>

template<int N>
void printApUint(ap_uint<N> toPrint){
	for(int i = N-1 ; i>=0; i--){
		fprintf(stderr, "%d", (int)toPrint[i]);
	}
	fprintf(stderr, "\n");
}

template<int N>
void printApInt(ap_int<N> toPrint){
	for(int i = N-1 ; i>=0; i--){
		fprintf(stderr, "%d", (int)toPrint[i]);
	}
	fprintf(stderr, "\n");
}

//In case we need both the variable and its reversed value
//(vivado reverse change the bit order of the value on which its called)
template <int N>
ap_uint<N> backward(
        ap_uint<N> input,
        typename std::enable_if<(N>1)>::type* = 0
    )
{
    #pragma HLS INLINE
    return static_cast<ap_uint<1> >(input[0]).concat(
                backward(
                    static_cast<ap_uint<N-1> >(
                        input.range(N-1, 1)
                        )
                    )
                );
}

template <int N>
ap_uint<N> backward(
        ap_uint<N> input,
        typename std::enable_if<(N==1)>::type* = 0
    )
{
    #pragma HLS INLINE
    return input;
}

#pragma once
#include "ap_int.h"

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

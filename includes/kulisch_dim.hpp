#ifndef FP_DIM_TPP
#define FP_DIM_TPP
#define AP_INT_MAX_W 16000

#include "ap_int.h"
#include "static_math.hpp"


#include "../includes/utils.hpp"

constexpr int ceilLog2(int N, uint8_t remains = 0)
{
	return (N <= 1) ? remains : 1 + ceilLog2(N>>1, remains | (N%2));
}

constexpr int ceil2Power(int N)
{
	return 1 << ceilLog2(N);
}

template<int N>
class Static_Val
{
	public:
		static constexpr int _log2 = ceilLog2(N);
		static constexpr int _2pow = ceil2Power(N);
};

template<int N>
class FPDim {
	public:
	static constexpr int WE = (N==32) ? 8 : 11; 
	static constexpr int WF = (N==32) ? 23 : 52;
	static constexpr int BIAS = (1<<(WE-1))-1;
	static constexpr int ACC_SIZE = (1<<(WE+1)) -1 + 2*WF+2;
	static constexpr int SUBNORMAL_LIMIT = (1<<(WE))+BIAS;
};


template<int N>
using FPProdSize = ap_uint<1+FPDim<N>::WE+1+2*FPDim<N>::WF+2>;

template<int N>
class FPProd : public FPProdSize<N>
{
	public:
		FPProd(	
				ap_uint<1> mult_s,
				ap_uint<FPDim<N>::WE+1> mult_e,
				ap_uint<2*FPDim<N>::WF+2> mult_m	
				)
		{
			ap_uint<1+FPDim<N>::WE+1> signed_exp = mult_s.concat(mult_e);
			FPProdSize<N>::operator=(signed_exp.concat(mult_m));		
		};

		FPProd(ap_uint<1+FPDim<N>::WE+1+2*FPDim<N>::WF+2> val):FPProdSize<N>(val){}

		ap_uint<2*FPDim<N>::WF+2> getSignificand()
		{
			#pragma HLS INLINE
			return FPProdSize<N>::range(2*FPDim<N>::WF+2-1, 0);
		}
		

		ap_uint<1> getSignBit()
		{
			#pragma HLS INLINE
			return (*this)[FPDim<N>::WE+1+2*FPDim<N>::WF+2];
		}

		ap_uint<FPDim<N>::WE+1> getExp()
		{
			#pragma HLS INLINE
			return FPProdSize<N>::range(
					FPDim<N>::WE+1+2*FPDim<N>::WF+2-1, 
					2*FPDim<N>::WF+2
				);
		}


		void printContent(){
			
			fprintf(stderr, "sign: ");
			printApUint(getSignBit());

			fprintf(stderr, "exp: ");
			printApUint(this->getExp());


			fprintf(stderr, "significand: ");
			printApUint(this->getSignificand());

		}
};

template<int N>
using KulischAcc = ap_uint<FPDim<N>::ACC_SIZE>;

template<int N, int bankSize>
static constexpr int getNbStages(){
	return Static_Ceil_Div<FPDim<N>::ACC_SIZE, bankSize>::val;
} 

template<int N, int bankSize>
static constexpr int getSegmentedAccSize(){
	return getNbStages<N, bankSize>() * bankSize;
} 

template<int N, int bankSize>
using SegmentedKulischAccSize = ap_uint<getSegmentedAccSize<N, bankSize>() + getNbStages<N, bankSize>()>;

template<int N, int bankSize>
static constexpr int getMantSpread(){
	return Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val+1;
} 


template<int N, int bankSize>
class SegmentedKulischAcc : public SegmentedKulischAccSize<N, bankSize>
{
	public:
		SegmentedKulischAcc(	
				KulischAcc<N> acc,
				ap_uint<getNbStages<N, bankSize>()> carries = 0
				)
		{
			ap_int<getSegmentedAccSize<N, bankSize>()> acc_ext = (ap_int<FPDim<N>::ACC_SIZE>) acc;
			SegmentedKulischAccSize<N, bankSize>::operator=(acc_ext.concat(carries));		
		};

		SegmentedKulischAcc(ap_uint<getSegmentedAccSize<N, bankSize>() + getNbStages<N, bankSize>()> val):SegmentedKulischAccSize<N, bankSize>(val){}

		ap_uint<bankSize> getBank(int index){
			#pragma HLS INLINE
			return (*this).range(getNbStages<N, bankSize>() + (index+1)*bankSize -1,getNbStages<N, bankSize>() + index*bankSize);
		}

		ap_uint<1> getCarry(int index){
			#pragma HLS INLINE
			return (*this)[index];
		}

		ap_uint<1> isNeg(){
			#pragma HLS INLINE
			return (*this)[getSegmentedAccSize<N, bankSize>() + getNbStages<N, bankSize>()-1];
		}

		void printContent(){
			fprintf(stderr, "Acc size: %d\n", getSegmentedAccSize<N, bankSize>());
			fprintf(stderr, "Nb stages: %d\n", getNbStages<N, bankSize>());
			ap_uint<getSegmentedAccSize<N, bankSize>() - FPDim<N>::ACC_SIZE> uselessbits = this->range(getSegmentedAccSize<N, bankSize>() + getNbStages<N, bankSize>()-1,FPDim<N>::ACC_SIZE + getNbStages<N, bankSize>());
			ap_uint<FPDim<N>::ACC_SIZE> tmp = this->range(FPDim<N>::ACC_SIZE + getNbStages<N, bankSize>() -1, getNbStages<N, bankSize>());
			printApUint(uselessbits);
			printApUint(tmp);
			ap_uint<getNbStages<N, bankSize>()> c = this->range(getNbStages<N, bankSize>()-1,0);
			printApUint(c);
		}
};


#endif
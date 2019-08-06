#ifndef FP_DIM_TPP
#define FP_DIM_TPP

#include <cstdint>
#include <type_traits>

#include "ap_int.h"
#include "tools/static_math.hpp"

using namespace std;

using hint::Static_Val;
using hint::Static_Ceil_Div;

template<int N>
class FPDim {
	public:
	static constexpr int WE = (N==32) ? 8 : ((N==64) ? 11 : 5); 
	static constexpr int WF = (N==32) ? 23 : ((N==64) ? 52 : 10);
	static constexpr int BIAS = (1<<(WE-1))-1;
	static constexpr int ACC_SIZE = (1<<(WE+1)) -1 + 2*WF+2;
	static constexpr int ACC_MID = (1<<(WE)) -1 + 2*WF+2;
	static constexpr int FP_SPREAD = (1<<(WE-1));
	static constexpr int PROD_FP_SPREAD = (1<<(WE));
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
		}

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

template<int N>
using SignedKulischAcc = ap_uint<FPDim<N>::ACC_SIZE+1>;


template<int N, int bankSize>
static constexpr int getNbStages(){
	return Static_Ceil_Div<FPDim<N>::ACC_SIZE, bankSize>::val;
} 

template<int N, int bankSize>
static constexpr int getSegmentedAccSize(){
	return getNbStages<N, bankSize>() * bankSize;
} 

// template<int N, int bankSize>
// using acc_2CK3Size = ap_uint<getSegmentedAccSize<N, bankSize>() + getNbStages<N, bankSize>()>;

template<int N, int bankSize>
static constexpr int getMantSpread(){
	return Static_Ceil_Div<2*FPDim<N>::WF+2,bankSize>::val+1;
} 


template<int N, int bankSize, int stage>
ap_uint<bankSize*(stage+1)> concatAccBanksRec(ap_uint<bankSize> acc[getNbStages<N, bankSize>()],
	typename enable_if<(stage == 0)>::type* = 0
	){
	return acc[stage];
}

template<int N, int bankSize, int stage>
ap_uint<bankSize*(stage+1)> concatAccBanksRec(ap_uint<bankSize> acc[getNbStages<N, bankSize>()],
	typename enable_if<(stage >= 1)>::type* = 0
	){
	ap_uint<bankSize*(stage)> part = concatAccBanksRec<N, bankSize, stage-1>(acc);
	ap_uint<bankSize> top = acc[stage];
	ap_uint<bankSize*(stage+1)> res = top.concat(part);
	return res;
}

template<int N, int bankSize>
KulischAcc<N> concatAccBanks(ap_uint<bankSize> acc[getNbStages<N, bankSize>()]){
	ap_uint<bankSize*getNbStages<N, bankSize>()> res = concatAccBanksRec<N, bankSize, getNbStages<N, bankSize>()-1>(acc);
	return res.range(FPDim<N>::ACC_SIZE-1, 0);
}


template<int N, int bankSize>
class acc_2CK3 
{
	private:
		ap_uint<bankSize> banks[getNbStages<N, bankSize>()];
		ap_uint<1> carries[getNbStages<N, bankSize>()];
	public:
		acc_2CK3(	
				KulischAcc<N> acc,
				ap_uint<getNbStages<N, bankSize>()> carries_in = 0
				)
		{
			
			#pragma HLS array_partition variable=banks
			#pragma HLS array_partition variable=carries
			for(int i=0; i<getNbStages<N, bankSize>()-1; i++){
				#pragma HLS UNROLL
				banks[i] = acc.range((i+1)*bankSize-1,i*bankSize);
				carries[i] = carries_in[i];
			}	
			banks[getNbStages<N, bankSize>()-1] = (ap_int<FPDim<N>::ACC_SIZE-1-(getNbStages<N, bankSize>()-1)*bankSize+1>)acc.range(FPDim<N>::ACC_SIZE-1,(getNbStages<N, bankSize>()-1)*bankSize);
			carries[getNbStages<N, bankSize>()-1] = carries_in[getNbStages<N, bankSize>()-1];	
		}

		acc_2CK3(	
				const acc_2CK3<N, bankSize> & acc_in
				)
		{
			
			#pragma HLS array_partition variable=banks
			#pragma HLS array_partition variable=carries
			for(int i=0; i<getNbStages<N, bankSize>(); i++){
				#pragma HLS UNROLL
				banks[i] = acc_in.getBank(i);
				carries[i] = acc_in.getCarry(i);
			}	
		}

		//acc_2CK3(ap_uint<getSegmentedAccSize<N, bankSize>() + getNbStages<N, bankSize>()> val):acc_2CK3Size<N, bankSize>(val){}

		KulischAcc<N> getAcc(){
			#pragma HLS INLINE
			return concatAccBanks<N, bankSize>(banks);
			// return (*this).range(FPDim<N>::ACC_SIZE + getNbStages<N, bankSize>()-1, getNbStages<N, bankSize>());
		}

		ap_uint<bankSize> getBank(int index) const{
			#pragma HLS INLINE
			return (*this).banks[index];
			// return (*this).range(getNbStages<N, bankSize>() + (index+1)*bankSize -1,getNbStages<N, bankSize>() + index*bankSize);
		}


		void setBank(int index, ap_uint<bankSize> bank){
			#pragma HLS INLINE
			(*this).banks[index] = bank;
			// (*this).range(getNbStages<N, bankSize>() + (index+1)*bankSize -1,getNbStages<N, bankSize>() + index*bankSize) = bank;
		}

		ap_uint<1> getCarry(int index) const{
			#pragma HLS INLINE
			return (*this).carries[index];
			// return (*this)[index];
		}

		void setCarry(int index, ap_uint<1> carry){
			#pragma HLS INLINE
			(*this).carries[index] = carry;
			// (*this)[index]=carry;
		}

		ap_uint<1> isNeg(){
			#pragma HLS INLINE
			return (*this).banks[getNbStages<N, bankSize>()-1][bankSize-1];
			// return (*this)[getSegmentedAccSize<N, bankSize>() + getNbStages<N, bankSize>()-1];
		}

		// void printContent(){
		// 	fprintf(stderr, "Acc size: %d\n", getSegmentedAccSize<N, bankSize>());
		// 	fprintf(stderr, "Nb stages: %d\n", getNbStages<N, bankSize>());
		// 	ap_uint<getSegmentedAccSize<N, bankSize>() - FPDim<N>::ACC_SIZE> uselessbits = this->range(getSegmentedAccSize<N, bankSize>() + getNbStages<N, bankSize>()-1,FPDim<N>::ACC_SIZE + getNbStages<N, bankSize>());
		// 	ap_uint<FPDim<N>::ACC_SIZE> tmp = this->range(FPDim<N>::ACC_SIZE + getNbStages<N, bankSize>() -1, getNbStages<N, bankSize>());
		// 	printApUint(uselessbits);
		// 	printApUint(tmp);
		// 	ap_uint<getNbStages<N, bankSize>()> c = this->range(getNbStages<N, bankSize>()-1,0);
		// 	printApUint(c);
		// }
};


template<int N, int bankSize>
using acc_segmented_2CK1 = acc_2CK3<N, bankSize>;

template<int N, int bankSize>
using acc_SMK3Size = ap_uint<getSegmentedAccSize<N, bankSize>() + getNbStages<N, bankSize>()+ getNbStages<N, bankSize>()>;

template<int N, int bankSize>
class acc_SMK3 : public acc_SMK3Size<N, bankSize>
{
	private:
		ap_uint<bankSize> banks[getNbStages<N, bankSize>()];
		ap_uint<1> carries[getNbStages<N, bankSize>()];
		ap_uint<1> borrows[getNbStages<N, bankSize>()];

	public:
		acc_SMK3(	
				KulischAcc<N> acc,
				ap_uint<getNbStages<N, bankSize>()> carries_in = 0,
				ap_uint<getNbStages<N, bankSize>()> borrows_in = 0
				)
		{


			#pragma HLS array_partition variable=banks
			#pragma HLS array_partition variable=carries
			#pragma HLS array_partition variable=borrows
			for(int i=0; i<getNbStages<N, bankSize>()-1; i++){
				#pragma HLS UNROLL
				banks[i] = acc.range((i+1)*bankSize-1,i*bankSize);
				carries[i] = carries_in[i];
				borrows[i] = borrows_in[i];
			}	
			banks[getNbStages<N, bankSize>()-1] = (ap_int<FPDim<N>::ACC_SIZE-1-(getNbStages<N, bankSize>()-1)*bankSize+1>)acc.range(FPDim<N>::ACC_SIZE-1,(getNbStages<N, bankSize>()-1)*bankSize);
			carries[getNbStages<N, bankSize>()-1] = carries_in[getNbStages<N, bankSize>()-1];	
			borrows[getNbStages<N, bankSize>()-1] = borrows_in[getNbStages<N, bankSize>()-1];		
		}


		acc_SMK3(	
				const acc_SMK3<N, bankSize> & acc_in
				)
		{
			
			#pragma HLS array_partition variable=banks
			#pragma HLS array_partition variable=carries
			#pragma HLS array_partition variable=borrows
			for(int i=0; i<getNbStages<N, bankSize>(); i++){
				#pragma HLS UNROLL
				banks[i] = acc_in.getBank(i);
				carries[i] = acc_in.getCarry(i);
				borrows[i] = acc_in.getBorrow(i);
			}	
		}

		KulischAcc<N> getAcc(){
			#pragma HLS INLINE
			return concatAccBanks<N, bankSize>(banks);
			// return (*this).range(FPDim<N>::ACC_SIZE + getNbStages<N, bankSize>()-1, getNbStages<N, bankSize>());
		}
		
		ap_uint<bankSize> getBank(int index) const{
			#pragma HLS INLINE
			return (*this).banks[index];
		}

		void setBank(int index, ap_uint<bankSize> bank){
			#pragma HLS INLINE
			(*this).banks[index] = bank;
		}

		ap_uint<1> getCarry(int index) const{
			#pragma HLS INLINE
			return (*this).carries[index];
			// return (*this)[index];
		}

		void setCarry(int index, ap_uint<1> carry){
			#pragma HLS INLINE
			(*this).carries[index] = carry;
		}

		ap_uint<1> getBorrow(int index) const{
			#pragma HLS INLINE
			return (*this).borrows[index];
			// return (*this)[index];
		}
		
		void setBorrow(int index, ap_uint<1> carry){
			#pragma HLS INLINE
			(*this).borrows[index] = carry;
		}

		ap_uint<1> isNeg(){
			#pragma HLS INLINE
			return (*this).banks[getNbStages<N, bankSize>()-1][bankSize-1];
			// return (*this)[getSegmentedAccSize<N, bankSize>() + getNbStages<N, bankSize>()-1];
		}

};
#endif

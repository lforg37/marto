#ifndef POSIT_DIM_TPP
#define POSIT_DIM_TPP

using namespace std;

#include "ap_int.h"
#include "utils.hpp"

//N.B.: We are using int instead of size_t because of ap_uint is templatized
//=====  with int

constexpr int get2Power(int N)
{
	return (N<=1) ? 0 : 1 + get2Power(N >> 1); 
}

constexpr int ceilLog2(int N, uint8_t remains = 0)
{
	return (N <= 1) ? remains : 1 + ceilLog2(N>>1, remains | (N%2));
}

constexpr int ceil2Power(int N)
{
	return 1 << ceilLog2(N);
}

template<int N>
static ap_uint<N> positiveMaxPosit() {
    return (((ap_uint<N>)1)<<(N-1))-2;
}

template<int N>
static ap_uint<N> negativeMaxPosit() {
    return (((ap_uint<N>)1)<<(N-1))+1;
}

template<int N>
static ap_uint<N> positiveMinPosit() {
    return ((ap_uint<N>)1);
}

template<int N>
static ap_uint<N> negativeMinPosit() {
    return ((ap_uint<N>)-1);
}

template<int N>
class PositDim {
	public:
	static constexpr int WES = get2Power(N>>3); 
	static constexpr int WE = get2Power(N) + WES + 1; 
	static constexpr int WF = N - (WES+3);
	//Quire dimension
	static constexpr int WQ = (N*N) >> 1;

	// The "2" is for the guard bit and the sticky 
	static constexpr int ValSize = 2 + 3 + WE + WF;
	static constexpr int ExtQuireSize = WQ + 1; //+1 for NaR bit
	static constexpr int ProdSignificandSize = 2*WF + 2; 
	//implicit bit twice 
	
	static constexpr int ProdExpSize = WE + 1;
	static constexpr int ProdSize = 2 + ProdExpSize + ProdSignificandSize; // + sign bit and isNaR

	static constexpr int EXP_BIAS = (N-2) * (1 << WES) + 1; //Add one because negative mantissa have an exponent shift of one compared to their opposite due to sign bit

	static constexpr bool HAS_ES = (WES > 0);

	// static constexpr int maxpos = maxPosit();
	// static constexpr ap_uint<N> minpos = ZERO_AND_ZEROS.concat(ONE_ONE);
	// static constexpr ap_uint<N> minusMaxpos = ONE_AND_ZEROS.concat(ONE_ONE);
	// static constexpr ap_uint<N> minusMinpos = ONE_AND_ONES.concat(ONE_ONE);
};

//One bit is NaR + quire
template<int N>
using QuireSizedAPUint = ap_uint<PositDim<N>::ExtQuireSize>;

template<int N>
class Quire : public QuireSizedAPUint<N>
{
	//Storage :
	// isNar Sign Carry 2sCompValue
	public:
		Quire(ap_uint<PositDim<N>::ExtQuireSize> val):QuireSizedAPUint<N>(val){}

		ap_uint<PositDim<N>::WQ-N> getQuireValue()
		{
			return QuireSizedAPUint<N>::range(PositDim<N>::WQ-N -1, 0);
		}

		ap_uint<N-1> getCarry()
		{
			return QuireSizedAPUint<N>::range(PositDim<N>::WQ-1 -1, 
					PositDim<N>::WQ-N);
		}

		ap_uint<PositDim<N>::ExtQuireSize-1> getQuireWithoutNaR()
		{
			return QuireSizedAPUint<N>::range(PositDim<N>::ExtQuireSize-1 -1, 0);
		}

		ap_uint<1> getSignBit()
		{
			return (*this)[PositDim<N>::ExtQuireSize-1 -1];
		}

		ap_uint<1> getIsNaR()
		{
			return (*this)[PositDim<N>::ExtQuireSize -1];
		}

		static constexpr int PositRangeOffset = ((N*N) >> 3) - (N >> 2); 
};


template<int N>
using PositProdSizedAPUint = ap_uint<PositDim<N>::ProdSize>;

// One bit isNar + WE+1 + 2(WF+1) 
template<int N>
class PositProd : public PositProdSizedAPUint<N>
{
	//Storage :
	// isNar Exp Signed_significand
	public:
		PositProd(
				ap_uint<1> isNar, 
				ap_uint<PositDim<N>::ProdExpSize> exp,
				ap_uint<1> sign,
				ap_uint<PositDim<N>::ProdSignificandSize> fraction
			) {
			ap_uint<PositDim<N>::ProdSignificandSize + 1> signed_frac = sign.concat(fraction);
			ap_uint<PositDim<N>::ProdExpSize + PositDim<N>::ProdSignificandSize + 1> prod =
				exp.concat(signed_frac);
			PositProdSizedAPUint<N>::operator=(isNar.concat(prod));
		}

		PositProd(ap_uint<PositDim<N>::ProdSize> val):PositProdSizedAPUint<N>(val){}

		ap_uint<PositDim<N>::ProdSignificandSize> getSignificand()
		{
			return PositProdSizedAPUint<N>::range(PositDim<N>::ProdSignificandSize - 1, 0);
		}
		
		ap_int<PositDim<N>::ProdSignificandSize + 1> getSignedSignificand()
		{
			return PositProdSizedAPUint<N>::range(PositDim<N>::ProdSignificandSize, 0);
		}

		ap_uint<1> getSignBit()
		{
			return (*this)[PositDim<N>::ProdSignificandSize];
		}

		ap_uint<PositDim<N>::ProdExpSize> getExp()
		{
			return PositProdSizedAPUint<N>::range(
					PositDim<N>::ProdSignificandSize + 1 + PositDim<N>::ProdExpSize-1, 
					PositDim<N>::ProdSignificandSize + 1
				);
		}

		ap_uint<1> getIsNaR()
		{
			return (*this)[PositDim<N>::ProdSize - 1];
		}

		void printContent(){

			fprintf(stderr, "isNaR: %d\n", (int) this->getIsNaR());
			
			fprintf(stderr, "biased exp: ");
			printApUint(this->getExp());

			fprintf(stderr, "sign: ");
			printApUint(getSignBit());

			fprintf(stderr, "significand: ");
			printApUint(this->getSignificand());

		}

	private:
};

template<int N>
using PositEncoding = ap_uint<N>;

template<int N>
using PositValSizedAPUint = ap_uint<PositDim<N>::ValSize>;

template<int N>
class PositValue : public PositValSizedAPUint<N>
{
	//Storage :
	// Guard Sticky isNar Exp Sign ImplicitBit Fraction
	public:
		PositValue(
				ap_uint<1> isNar,
				ap_uint<PositDim<N>::WE> exp, //Warning : biased exp
				ap_uint<1> sign,
				ap_uint<1> implicit_bit,
				ap_uint<PositDim<N>::WF> fraction)
		{
			ap_uint<1+PositDim<N>::WE> tmp = isNar.concat(exp); 
			ap_uint<2> frac_lead = sign.concat(implicit_bit);
			ap_uint<2+PositDim<N>::WF> full_frac = frac_lead.concat(fraction);
			PositValSizedAPUint<N>::operator=(tmp.concat(full_frac));
		}

		PositValue(
				ap_uint<1> guard,
				ap_uint<1> sticky,
				ap_uint<1> isNar,
				ap_uint<PositDim<N>::WE> exp,
				ap_uint<1> sign,
				ap_uint<1> implicit_bit,
				ap_uint<PositDim<N>::WF> fraction)
		{
			ap_uint<2> guardAndSticky = guard.concat(sticky);
			ap_uint<1+PositDim<N>::WE> tmp = isNar.concat(exp); 
			ap_uint<2> frac_lead = sign.concat(implicit_bit);
			ap_uint<2+PositDim<N>::WF> full_frac = frac_lead.concat(fraction);
			ap_uint<2+PositDim<N>::WF+1+PositDim<N>::WE> allWoGuardAndSticky = tmp.concat(full_frac);
			PositValSizedAPUint<N>::operator=(guardAndSticky.concat(allWoGuardAndSticky));
		}
				
		PositValue(ap_uint<PositDim<N>::ValSize> val):PositValSizedAPUint<N>(val){}

		ap_uint<1> getGuardBit()
		{
			return (*this)[PositDim<N>::ValSize-1];
		}

		ap_uint<1> getStickyBit()
		{
			return (*this)[PositDim<N>::ValSize-1 -1];
		}


		ap_uint<PositDim<N>::WF+1> getSignificand()
		{
			return PositValSizedAPUint<N>::range(PositDim<N>::WF, 0);
		}

		ap_uint<1> getImplicitBit()
		{
			return (*this)[PositDim<N>::WF];
		}

		ap_uint<PositDim<N>::WF> getSignificandWoImp()
		{
			return PositValSizedAPUint<N>::range(PositDim<N>::WF-1, 0);
		}

		ap_uint<1> getSignBit()
		{
			return (*this)[PositDim<N>::WF + 1];
		}

		ap_uint<PositDim<N>::WE> getExp()
		{
			return PositValSizedAPUint<N>::range(PositDim<N>::WF + 1 + PositDim<N>::WE, PositDim<N>::WF+2);
		}

		ap_uint<1> getIsNaR()
		{
			return (*this)[PositDim<N>::WF+2+PositDim<N>::WE];
		}

		ap_int<PositDim<N>::WF+2> getSignedSignificand()
		{
			ap_uint<1> sign = getSignBit();
			return (ap_int<PositDim<N>::WF+2>) sign.concat(getSignificand());
		}

		ap_uint<1> getBit(unsigned int i)
		{
			return (*this)[i];
		}

		ap_uint<1> isZero()
		{
			ap_uint<1> isExponentNull = not(getExp().or_reduce());
			return isExponentNull and (not getSignBit());
		}

		void printContent(){
			fprintf(stderr, "guard: %d\n", (int)this->getGuardBit());
			fprintf(stderr, "sticky: %d\n", (int)this->getStickyBit());

			fprintf(stderr, "isNaR: %d\n", (int) this->getIsNaR());
			
			fprintf(stderr, "biased exp: ");
			printApUint(this->getExp());

			fprintf(stderr, "sign: %d\n", (int) this->getSignBit());

			fprintf(stderr, "significand: %d.", (int) (this->getSignificand())[PositDim<N>::WF+1 -1]);
			printApUint(this->getSignificandWoImp());

			double temp = getSignedSignificand().to_int();
			double exp = pow(2, getExp().to_int() - PositDim<N>::WF - PositDim<N>::EXP_BIAS);
			fprintf(stderr, "Value : %1.20f\n", temp*exp);
		}

		static PositValue getMaxPos() 
		{ //isNar Exp Sign Implicit Frac
			return PositValue(
					0, //isNar
					2*PositDim<N>::EXP_BIAS - 1, //Biased Exp
					0, //sign
					1, //implicit bit
					0 //fraction
				);
		}	

		static PositValue getMinPos()
		{
			return PositValue(
					0, //isNar
					1, // Biased exp
					0, // sign
					1, // implicit bit
					0 // fraction
				);
		}

		static PositValue getMaxNeg()
		{
			return PositValue(
					0,
					2*PositDim<N>::EXP_BIAS - 2, // Biased Exp
					1, //sign
					0, //implicit bit
					0  //fraction
				);
		}

		static PositValue getMinNeg()
		{
			return PositValue(
				0, //isNar
				0, // Biased Exp
				1, // sign
				0, // implicit bit
				0 // fraction
			);
		}
};

#endif

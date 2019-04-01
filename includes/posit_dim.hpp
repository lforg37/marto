#pragma once

#include <cstdint>
#include <cstdio>

#include "ap_int.h"
#include "utils.hpp"
#include "static_math.hpp"
//N.B.: We are using int instead of size_t because of ap_uint is templatized
//=====  with int
// #include <boost/integer/static_log2.hpp>

// #define get2Power(N) boost::static_log2<N>::value
// constexpr int get2Power(int N)
// {
// 	// unsigned int result = 0;
// 	// while(N>1){
// 	// 	result+=1;
// 	// 	N=N>>1;
// 	// }
// 	// return result;
// 	return ::boost::static_log2<N>::value;
// 	// return (N<=1) ? 0 : 1 + get2Power(N >> 1); 
// }

using namespace std;

constexpr int ceilLog2(int N, uint8_t remains = 0)
{
	return (N <= 1) ? remains : 1 + ceilLog2(N>>1, remains | (N%2));
}

constexpr int ceil2Power(int N)
{
	return 1 << ceilLog2(N);
}

constexpr int log2(int N)
{
	 return ((N<2) ? 1 : 1+log2(N>>1));
}

constexpr int r2pow(int N)
{
	 return 1 << log2(N);
}



template<int N>
class Static_Val
{
	public:
		static constexpr int _rlog2 = log2(N);
		static constexpr int _r2pow = r2pow(N);
		static constexpr int _log2 = ceilLog2(N);
		static constexpr int _2pow = ceil2Power(N);
};

template<int N>
static ap_uint<N> positiveMaxPosit() {
    return (ap_uint<N>{1} << (N-1)) - 2;
}

template<int N>
static ap_uint<N> negativeMaxPosit() {
    return (ap_uint<N>{1} << (N-1))+1;
}

template<int N>
static ap_uint<N> positiveMinPosit() {
    return ap_uint<N>{1};
}

template<int N>
static ap_uint<N> negativeMinPosit() {
    return ap_uint<N>{-1};
}




template<int N, int WES_Val>
class PositDim {
	public:
    static constexpr int WES = WES_Val;
	// get2Power((N>>3)); 
    static constexpr int WE = Static_Val<N>::_log2 + WES_Val + 1;
    static constexpr int WF = N - (WES_Val+3);

	// The "2" is for the guard bit and the sticky 
	static constexpr int ValSize = 2 + 3 + WE + WF;
    static constexpr int EMax = (N-2) * (1 << WES_Val);
	static constexpr int ProdSignificandSize = 2*WF + 2; 
	//implicit bit twice 
	
	static constexpr int ProdExpSize = WE + 1;
	static constexpr int ProdSize = 2 + ProdExpSize + ProdSignificandSize; // + sign bit and isNaR

    static constexpr int EXP_BIAS = EMax + 1; //Add one because negative mantissa have an exponent shift of one compared to their opposite due to sign bit

    static constexpr bool HAS_ES = (WES_Val > 0);
};

template <int N>
using StandardPositDim = PositDim<N, Static_Val<(N>>3)>::_log2>;

//One bit is NaR + quire
template<int N, int WES, int NB_CARRY>
using QuireSizedAPUint = ap_uint<PositDim<N, WES>::EMax * 4 + 3 + NB_CARRY>; // + 3 : Sign, isNar, 0 exp

template <int N, int WES>
class PositValue;
template <int N, int WES>
class PositProd;

template<int N, int WES, int NB_CARRY = 1>
class Quire : public QuireSizedAPUint<N, WES, NB_CARRY>
{
	//Storage :
	// isNar Sign Carry 2sCompValue
	public:
        static constexpr int ProductRangeSize = PositDim<N, WES>::EMax * 4 + 1; //Zero exp and
        static constexpr int Size =  ProductRangeSize + 2 + NB_CARRY; //isNaR bit + sign bit
        Quire(ap_uint<Size> val):QuireSizedAPUint<N, WES, NB_CARRY>(val){}

        ap_uint<ProductRangeSize> getQuireValue() const
		{
			#pragma HLS INLINE
            return QuireSizedAPUint<N, WES, NB_CARRY>::range(ProductRangeSize - 1, 0);
		}

        ap_uint<NB_CARRY> getCarry() const
		{
			#pragma HLS INLINE
            return QuireSizedAPUint<N, WES, NB_CARRY>::range(ProductRangeSize + NB_CARRY - 1,
                    ProductRangeSize);
		}

        ap_uint<Size-1> getQuireWithoutNaR() const
		{
			#pragma HLS INLINE
            return QuireSizedAPUint<N, WES, NB_CARRY>::range(Size-2, 0);
		}

        ap_uint<1> getSignBit() const
		{
			#pragma HLS INLINE
            return (*this)[Size - 2];
		}

        ap_uint<1> getIsNaR() const
		{
			#pragma HLS INLINE
            return (*this)[Size - 1];
		}

        void printContent() const{
            printApUint(*this);
		}

        operator PositValue<N, WES>() const;

        /*
         * Quire organisation
         *
         *  _______________________________________________________
         * !NaR|s| overflow        |  valid posits    | underflow  !
         * --------------------------------------------------------
         * Overflow  : weights > max pos and carry bits
         * Underflow : weights < minPos
         */



        // Width of underflow region
        static constexpr int PositRangeOffset = PositDim<N, WES>::EMax;
        // Width of valid posit region
		static constexpr int PositExpRange = 2*PositRangeOffset + 1;
        // Position of first overflow bit
        static constexpr int PositOverflowOffset = PositRangeOffset + PositExpRange;

};

template <int N, int WES, int NB_CARRY>
Quire<N, WES, NB_CARRY> operator+(
        Quire<N, WES, NB_CARRY> const & lhs,
        PositProd<N, WES> const & rhs
    );
template <int N, int WES, int NB_CARRY>
Quire<N, WES, NB_CARRY> operator-(
        Quire<N, WES, NB_CARRY> const & lhs,
        PositProd<N, WES> const & rhs
    );

template <int N>
using StandardQuire = Quire<N, Static_Val<(N>>3)>::_log2, N-2>;

template<int N, int WES, int NB_CARRY, int bankSize>
static constexpr int getNbStages(){
    return Static_Ceil_Div<Quire<N, WES, NB_CARRY>::Size - 1, bankSize>::val;
} 

template<int bankSize>
static constexpr int getShiftSize(){
	return Static_Val<bankSize>::_log2;
} 

template<int N, int WES, int bankSize>
static constexpr int getExtShiftSize(){
    return (Static_Ceil_Div<PositDim<N, WES>::ProdSignificandSize+1+(1<<getShiftSize<bankSize>()),bankSize>::val)*bankSize ;
} 

template<int N, int WES, int bankSize>
static constexpr int getMantSpread(){
    return Static_Ceil_Div<getExtShiftSize<N, WES, bankSize>(),bankSize>::val;
} 

// Stored as normal quire then carry bits from most significant to less
template<int N, int WES, int NB_CARRY, int bankSize>
class SegmentedQuire : public ap_uint<Quire<N, WES, NB_CARRY>::Size+getNbStages<N, WES, NB_CARRY, bankSize>()>
{
	//Storage :
	// isNar Sign Carry 2sCompValue
	public:
        static constexpr int Nb_stages = getNbStages<N, WES, NB_CARRY, bankSize>();
        static constexpr int Size = Quire<N, WES, NB_CARRY>::Size + Nb_stages;
        static constexpr int PositRangeOffset = ((N*N) >> 3) - (N >> 2);

        SegmentedQuire(ap_uint<Size> val):ap_uint<Size>(val){}

        ap_uint<bankSize> getBank(int index) const{
			#pragma HLS INLINE
			// fprintf(stderr, "== quire bank (%d, %d) ==\n", getNbStages<N, bankSize>() + (index+1)*bankSize -1, getNbStages<N, bankSize>() + index*bankSize);
			// fprintf(stderr, "From getNbStages<N, bankSize>()> (%d) + (index(%d)+1)*bankSize (%d) -1\n", getNbStages<N, bankSize>()  , index, bankSize);
			// fprintf(stderr, "To getNbStages<N, bankSize>()> + index*bankSize\n");
			// ap_uint<bankSize> tmp = (*this).range(getNbStages<N, bankSize>() + (index+1)*bankSize -1,getNbStages<N, bankSize>() + index*bankSize);
			// printApUint(tmp);
            return (*this).range(Nb_stages + (index+1)*bankSize -1, Nb_stages + index*bankSize);
		}

        ap_uint<1> getCarry(int index) const{
			#pragma HLS INLINE
			return (*this)[index];
		}

        ap_uint<1> getIsNaR() const{
			#pragma HLS INLINE
            return (*this)[Size-1];
		}

        ap_uint<Quire<N, WES, NB_CARRY>::Size> getAsQuireWoCarries() const {
			#pragma HLS INLINE
            return (*this).range(Size-1, Nb_stages);
		}

        ap_uint<Nb_stages> getAllCarries() const
        {
            return this->range(Nb_stages-1,0);
        }

        void printContent()  const {
            fprintf(stderr, "Quire size: %d\n", Quire<N, WES, NB_CARRY>::Size);
            fprintf(stderr, "Nb stages: %d\n", Nb_stages);
            printApUint(getAsQuireWoCarries());
            printApUint(getAllCarries());
		}

        operator PositValue<N, WES>() const;

};

template<int N, int WES, int NB_CARRY, int bankSize>
SegmentedQuire<N, WES, NB_CARRY, bankSize> operator+(
        SegmentedQuire<N, WES, NB_CARRY, bankSize> const & lhs,
        PositProd<N, WES> const & rhs
    );
template<int N, int WES, int NB_CARRY, int bankSize>
SegmentedQuire<N, WES, NB_CARRY, bankSize> operator-(
        SegmentedQuire<N, WES, NB_CARRY, bankSize> const & lhs,
        PositProd<N, WES> const & rhs
    );

template<int N, int banksize>
using StandardSegmentedQuire = SegmentedQuire<N, Static_Val<(N>>3)>::_log2, N-2, banksize>;

template <int N, int WES>
class PositEncoding;

template<int N, int WES>
using PositProdSizedAPUint = ap_uint<PositDim<N, WES>::ProdSize>;

// One bit isNar + WE+1 + 2(WF+1) 
template<int N, int WES>
class PositProd : public PositProdSizedAPUint<N, WES>
{
	//Storage :
	// isNar Exp Signed_significand
	public:
        static constexpr int Size = PositDim<N, WES>::ProdSize;
        static constexpr int ExpSize = PositDim<N, WES>::ProdExpSize;
        static constexpr int SignificandSize = PositDim<N, WES>::ProdSignificandSize;
		PositProd(
				ap_uint<1> isNar, 
                ap_uint<ExpSize> exp,
				ap_uint<1> sign,
                ap_uint<SignificandSize> fraction
			) {
            ap_uint<SignificandSize + 1> signed_frac = sign.concat(fraction);
            ap_uint<ExpSize + SignificandSize + 1> prod = exp.concat(signed_frac);
            PositProdSizedAPUint<N, WES>::operator=(isNar.concat(prod));
		}

        PositProd(ap_uint<Size> val):PositProdSizedAPUint<N, WES>(val){}

        ap_uint<SignificandSize> getSignificand() const
		{
			#pragma HLS INLINE
            return this->range(SignificandSize - 1, 0);
		}
		
        ap_int<SignificandSize + 1> getSignedSignificand() const
		{
			#pragma HLS INLINE
            return this->range(SignificandSize, 0);
		}

        ap_uint<1> getSignBit() const
		{
			#pragma HLS INLINE
            return (*this)[SignificandSize];
		}

        ap_uint<ExpSize> getExp() const
		{
			#pragma HLS INLINE
            return this->range(
                    SignificandSize + ExpSize,
                    SignificandSize + 1
				);
		}

        ap_uint<1> getIsNaR() const
		{
			#pragma HLS INLINE
            return (*this)[Size - 1];
		}

        operator PositValue<N, WES>() const;
        operator PositEncoding<N, WES>() const;

        void printContent() const {
            fprintf(stderr, "isNaR: %d\n", static_cast<int>(this->getIsNaR()));
			
			fprintf(stderr, "biased exp: ");
            printApUint(getExp());

			fprintf(stderr, "sign: ");
			printApUint(getSignBit());

			fprintf(stderr, "significand: ");
            printApUint(getSignificand());
		}
};
template<int N>
using StandardPositProd = PositProd<N, Static_Val<(N>>3)>::_log2>;

template<int N, int WES>
using PositValSizedAPUint = ap_uint<PositDim<N, WES>::ValSize>;

template<int N, int WES>
class PositValue : public PositValSizedAPUint<N, WES>
{
	//Storage :
	// Guard Sticky isNar Exp Sign ImplicitBit Fraction
	public:
        static constexpr int Size = PositDim<N, WES>::ValSize;
        static constexpr int ExpSize = PositDim<N, WES>::WE;
        static constexpr int FractionSize = PositDim<N, WES>::WF;
		PositValue(
				ap_uint<1> isNar,
                ap_uint<ExpSize> exp, //Warning : biased exp
				ap_uint<1> sign,
				ap_uint<1> implicit_bit,
                ap_uint<FractionSize> fraction)
		{
            ap_uint<1+ExpSize> tmp = isNar.concat(exp);
			ap_uint<2> frac_lead = sign.concat(implicit_bit);
            ap_uint<2+FractionSize> full_frac = frac_lead.concat(fraction);
            PositValSizedAPUint<N, WES>::operator=(tmp.concat(full_frac));
		}

		PositValue(
				ap_uint<1> guard,
				ap_uint<1> sticky,
				ap_uint<1> isNar,
                ap_uint<ExpSize> exp,
				ap_uint<1> sign,
				ap_uint<1> implicit_bit,
                ap_uint<FractionSize> fraction)
		{
			ap_uint<2> guardAndSticky = guard.concat(sticky);
            ap_uint<1+ExpSize> tmp = isNar.concat(exp);
			ap_uint<2> frac_lead = sign.concat(implicit_bit);
            ap_uint<2+FractionSize> full_frac = frac_lead.concat(fraction);
            ap_uint<2+FractionSize+1+ExpSize> allWoGuardAndSticky = tmp.concat(full_frac);
            PositValSizedAPUint<N, WES>::operator=(guardAndSticky.concat(allWoGuardAndSticky));
		}
				
        PositValue(ap_uint<Size> val):PositValSizedAPUint<N, WES>(val){}

        ap_uint<1> getGuardBit() const
		{
			#pragma HLS INLINE
            return (*this)[Size-1];
		}

        ap_uint<1> getStickyBit() const
		{
			#pragma HLS INLINE
            return (*this)[Size-2];
		}


        ap_uint<FractionSize + 1> getSignificand() const //Implicit bit + fractional part
		{
			#pragma HLS INLINE
            return this->range(FractionSize, 0);
		}

        ap_uint<1> getImplicitBit() const
		{
			#pragma HLS INLINE
            return (*this)[FractionSize];
		}

        ap_uint<FractionSize> getFraction() const
		{
			#pragma HLS INLINE
            return this->range(FractionSize-1, 0);
		}

        ap_uint<1> getSignBit() const
		{
			#pragma HLS INLINE
            return (*this)[FractionSize + 1];
		}

        ap_uint<ExpSize> getExp() const
		{
			#pragma HLS INLINE
            return this->range(FractionSize + 1 + ExpSize, FractionSize+2);
		}

        ap_uint<1> getIsNaR() const
		{
			#pragma HLS INLINE
            return (*this)[FractionSize+2+ExpSize];
		}

        ap_int<FractionSize+2> getSignedSignificand() const
		{
			#pragma HLS INLINE
			ap_uint<1> sign = getSignBit();
            return static_cast<ap_int<FractionSize+2> >(sign.concat(getSignificand()));
		}

        ap_uint<1> isZero() const
		{
			#pragma HLS INLINE
			ap_uint<1> isExponentNull = not(getExp().or_reduce());
			return isExponentNull and (not getSignBit());
		}

        void printContent() const {
            fprintf(stderr, "guard: %d\n", static_cast<int>(this->getGuardBit()));
            fprintf(stderr, "sticky: %d\n", static_cast<int>(this->getStickyBit()));

            fprintf(stderr, "isNaR: %d\n", static_cast<int>(this->getIsNaR()));
			
			fprintf(stderr, "biased exp: ");
			printApUint(this->getExp());

            fprintf(stderr, "sign: %d\n", static_cast<int>(this->getSignBit()));

            fprintf(stderr, "significand: %d.", static_cast<int>(this->getImplicitBit()));
            printApUint(this->getFraction());

			double temp = getSignedSignificand().to_int();
            double exp = pow(2, getExp().to_int() - FractionSize - PositDim<N, WES>::EXP_BIAS);
			fprintf(stderr, "Value : %1.20f\n", temp*exp);
		}

        operator PositProd<N, WES>() const;
        operator PositEncoding<N, WES>() const;

		static PositValue getMaxPos() 
		{ //isNar Exp Sign Implicit Frac
			return PositValue(
					0, //isNar
                    2*PositDim<N, WES>::EXP_BIAS - 1, //Biased Exp
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
                    2*PositDim<N, WES>::EXP_BIAS - 2, // Biased Exp
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

template <int N>
using StandardPositValue = PositValue<N, Static_Val<(N>>3)>::_log2>;


template<int N, int WES>
class PositEncoding : public ap_uint<N>
{
public:
    PositEncoding(ap_uint<N> val):ap_uint<N>{val}{}

    operator PositValue<N, WES>() const;
    operator PositProd<N, WES>() const;
//    template<int NB_CARRY>
//    operator Quire<N, WES, NB_CARRY>() const;
//    template<int NB_CARRY, int banksize>
//    operator SegmentedQuire<N, WES, NB_CARRY, banksize>() const;
};

template<int N, int WES>
PositEncoding<N, WES> operator+(
        PositEncoding<N, WES> const & lhs,
        PositEncoding<N, WES> const & rhs
   );

template<int N, int WES>
PositEncoding<N, WES> operator-(
        PositEncoding<N, WES> const & lhs,
        PositEncoding<N, WES> const & rhs
    );
template<int N, int WES>
PositEncoding<N, WES> operator-(
        PositEncoding<N, WES> const & val
    )
{
	return  PositEncoding<N, WES> ( ~(static_cast< ap_uint<N> >(val))+1 );
};

template<int N, int WES>
PositProd<N, WES> operator*(
        PositEncoding<N, WES> const & lhs,
        PositEncoding<N, WES> const & rhs
    );

template <int N>
using StandardPositEncoding = PositEncoding<N, Static_Val<(N>>3)>::_log2>;

#include "conversions.inc"

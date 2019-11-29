#pragma once
#include <cstdint>
#include "tools/static_math.hpp"

template<unsigned int N, template<unsigned int, bool> class Wrapper>
static Wrapper<N, false> positiveMaxPosit() {
	return Wrapper<N, false>{(1 << (N-1)) - 2};
}

template<unsigned int N, template<unsigned int, bool> class Wrapper>
static Wrapper<N, false> negativeMaxPosit() {
	return Wrapper<N, false>{(1 << (N-1))+1};
}

template<unsigned int N, template<unsigned int, bool> class Wrapper>
static Wrapper<N, false> positiveMinPosit() {
	return Wrapper<N, false>{1};
}

template<unsigned int N, template<unsigned int, bool> class Wrapper>
static Wrapper<N, false> negativeMinPosit() {
	return Wrapper<N, false>::generateSequence(Wrapper<1, false>{1});
}

template<unsigned int N, unsigned int WES_Val>
class PositDim {
	public:
	static constexpr unsigned int WES = WES_Val;
	// get2Power((N>>3));
	static constexpr unsigned int WE = hint::Static_Val<N>::_log2 + WES_Val + 1;
	static constexpr unsigned int WF = N - (WES_Val+3);

	// The "2" is for the guard bit and the sticky
	static constexpr unsigned int ValSize = 2 + 3 + WE + WF;
	static constexpr unsigned int EMax = (N-2) * (1 << WES_Val);
	static constexpr unsigned int ProdSignificandSize = 2*WF + 2;
	//implicit bit twice

	static constexpr unsigned int ProdExpSize = WE + 1;
	static constexpr unsigned int ProdSize = 2 + ProdExpSize + ProdSignificandSize; // + sign bit and isNaR

	static constexpr unsigned int EXP_BIAS = EMax + 1; //Add one because negative mantissa have an exponent shift of one compared to their opposite due to sign bit

	static constexpr bool HAS_ES = (WES_Val > 0);
};

template <unsigned int N>
using StandardPositDim = PositDim<N, hint::Static_Val<(N>>3)>::_log2>;


template<unsigned int N, unsigned int WES, unsigned int NB_CARRY>
struct QuireDim
{
	static constexpr unsigned int ProductRangeSize = PositDim<N, WES>::EMax * 4 + 1;
	static constexpr unsigned int Size = ProductRangeSize + 2 + NB_CARRY;
};


template <unsigned int N>
using StandardQuireDim = QuireDim<N, hint::Static_Val<(N>>3)>::_log2, N-2>;

//One bit is NaR + quire
template<int N, int WES, int NB_CARRY, template<unsigned int, bool> class Wrapper>
using QuireSizedHint = Wrapper<QuireDim<N, WES, NB_CARRY>::Size, false>; // + 3 : Sign, isNar, 0 exp

template <unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, bool isExact>
class PositIntermediateFormat;
template <unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
class PositProd;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY = 1>
class Quire : public QuireSizedHint<N, WES, NB_CARRY, Wrapper>
{
	//Storage :
	// isNar Sign Carry 2sCompValue
		using _storage = QuireSizedHint<N, WES, NB_CARRY, Wrapper>;
	public:
		static constexpr int ProductRangeSize = QuireDim<N, WES, NB_CARRY>::ProductRangeSize; //Zero exp and
		static constexpr int Size = QuireDim<N, WES, NB_CARRY>::Size; //isNaR bit + sign bit

		Quire(_storage val):_storage{val}{}
		Quire():_storage{}{}

		inline Wrapper<ProductRangeSize, false> getQuireValue() const
		{
			return _storage::template slice<ProductRangeSize - 1, 0>();
		}

		inline Wrapper<NB_CARRY, false> getCarry() const
		{
			return _storage::template slice<ProductRangeSize + NB_CARRY - 1,
					ProductRangeSize>();
		}

		inline Wrapper<Size-1, false> getQuireWithoutNaR() const
		{
			return _storage::template slice<Size-2, 0>();
		}

		inline Wrapper<1, false> getSignBit() const
		{
			return _storage::template get<Size - 2>();
		}

		inline Wrapper<1, false> getIsNaR() const
		{
			return _storage::template get<Size - 1>();
		}

		operator PositIntermediateFormat<N, WES, Wrapper, false>() const;

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

template <int N, int WES, int NB_CARRY, template<unsigned int, bool> class Wrapper>
Quire<N, WES, Wrapper, NB_CARRY> operator+(
		Quire<N, WES, Wrapper, NB_CARRY> const & lhs,
		PositProd<N, WES, Wrapper> const & rhs
	);
template <int N, int WES, int NB_CARRY, template<unsigned int, bool> class Wrapper>
Quire<N, WES, Wrapper, NB_CARRY> operator-(
		Quire<N, WES, Wrapper, NB_CARRY> const & lhs,
		PositProd<N, WES, Wrapper> const & rhs
	);

template <unsigned int N, template<unsigned int, bool> class Wrapper>
using StandardQuire = Quire<N, hint::Static_Val<(N>>3)>::_log2, Wrapper, N-2>;


template<unsigned int N, unsigned int WES, unsigned int NB_CARRY, unsigned int bankSize>
static constexpr unsigned int getNbStages(){
	return hint::Static_Ceil_Div<QuireDim<N, WES, NB_CARRY>::Size - 1, bankSize>::val;
}

template<unsigned int bankSize>
static constexpr unsigned int getShiftSize(){
	return hint::Static_Val<bankSize>::_log2;
}

template<unsigned int N, unsigned int WES, unsigned int bankSize>
static constexpr unsigned int getExtShiftSize(){
	return (hint::Static_Ceil_Div<PositDim<N, WES>::ProdSignificandSize+1+(1<<getShiftSize<bankSize>()),bankSize>::val)*bankSize ;
}

template<unsigned int N, unsigned int WES, unsigned int bankSize>
static constexpr unsigned int getMantSpread(){
	return hint::Static_Ceil_Div<getExtShiftSize<N, WES, bankSize>(), bankSize>::val;
}

template<unsigned int N, unsigned int WES, unsigned int NB_CARRY, unsigned int bankSize, template<unsigned int, bool> class Wrapper>
using SegmentedQuireSizedHint = Wrapper<QuireDim<N, WES, NB_CARRY>::Size + getNbStages<N, WES, NB_CARRY, bankSize>(), false>;

// Stored as normal quire then carry bits from most significant to less
template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int bankSize>
class SegmentedQuire : public SegmentedQuireSizedHint<N, WES, NB_CARRY, bankSize, Wrapper>
{
	//Storage :
	// isNar Sign Carry 2sCompValue
	using _storage = SegmentedQuireSizedHint<N, WES, NB_CARRY, bankSize, Wrapper>;
	public:
		static constexpr unsigned int Nb_stages = getNbStages<N, WES, NB_CARRY, bankSize>();
		static constexpr unsigned int Size = QuireDim<N, WES, NB_CARRY>::Size + Nb_stages;
		static constexpr unsigned int PositRangeOffset = ((N*N) >> 3) - (N >> 2);

		SegmentedQuire(_storage val):_storage{val}{}
		SegmentedQuire():_storage{}{}

		template<unsigned int BankIdx>
		inline Wrapper<bankSize, false> getBank() const{
			static_assert (BankIdx < Nb_stages, "Trying to access a bank with an out of range index");
			constexpr unsigned int low_idx = Nb_stages + BankIdx*bankSize;
			constexpr unsigned int up_idx = low_idx + bankSize - 1;
			return _storage::template slice<up_idx, low_idx>();
		}

		template<unsigned int BankIdx>
		inline Wrapper<1, false> getCarry() const{
			static_assert (BankIdx < Nb_stages, "Trying to access a carry with an out of range index");
			return _storage::template get<BankIdx>();
		}

		inline Wrapper<1, false> getIsNaR() const{
			return _storage::template get<Size-1>();
		}

		inline Wrapper<QuireDim<N, WES, NB_CARRY>::Size, false> getAsQuireWoCarries() const {
			return _storage::template slice<Size-1, Nb_stages>();
		}

		inline Wrapper<Nb_stages, false> getAllCarries() const
		{
			return _storage::template slice<Nb_stages-1,0>();
		}

		operator PositIntermediateFormat<N, WES, Wrapper, false>() const;

};

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int bankSize>
SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> operator+(
		SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> const & lhs,
		PositProd<N, WES, Wrapper> const & rhs
	);
template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper, unsigned int NB_CARRY, unsigned int bankSize>
SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> operator-(
		SegmentedQuire<N, WES, Wrapper, NB_CARRY, bankSize> const & lhs,
		PositProd<N, WES, Wrapper> const & rhs
	);

template<unsigned int N, unsigned int banksize, template<unsigned int, bool> class Wrapper>
using StandardSegmentedQuire = SegmentedQuire<N, hint::Static_Val<(N>>3)>::_log2, Wrapper, N-2, banksize>;


template <unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
class PositEncoding;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
using PositProdSizedHint = Wrapper<PositDim<N, WES>::ProdSize, false>;

// One bit isNar + WE+1 + 2(WF+1)
template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
class PositProd : public PositProdSizedHint<N, WES, Wrapper>
{
	//Storage :
	// isNar Exp Signed_significand
	public:
		typedef PositProdSizedHint<N, WES, Wrapper> hint_type;
		template<unsigned int W>
		using wrapper_helper = Wrapper<W, false>;
		static constexpr unsigned int Size = PositDim<N, WES>::ProdSize;
		static constexpr unsigned int ExpSize = PositDim<N, WES>::ProdExpSize;
		static constexpr unsigned int SignificandSize = PositDim<N, WES>::ProdSignificandSize;
		PositProd(
				wrapper_helper<1> isNar,
				wrapper_helper<ExpSize> exp,
				wrapper_helper<1> sign,
				wrapper_helper<SignificandSize> fraction
			) : hint_type{isNar.concatenate(exp.concatenate(sign.concatenate(fraction)))}
		{}

		PositProd(hint_type val):hint_type{val}{}
		// PositProd(PositIntermediateFormat<N, WES> val):PositProdSizedAPUint<N, WES>(val){}

		inline wrapper_helper<SignificandSize> getSignificand() const
		{
			return hint_type::template slice<SignificandSize-1, 0>();
		}

		inline wrapper_helper<SignificandSize + 1> getSignedSignificand() const
		{
			return hint_type::template slice<SignificandSize, 0>();
		}

		inline wrapper_helper<1> getSignBit() const
		{
			return hint_type::template get<SignificandSize>();
		}

		inline wrapper_helper<ExpSize> getExp() const
		{
			return hint_type::template slice<ExpSize+SignificandSize, SignificandSize + 1>();
		}

		inline wrapper_helper<1> getIsNaR() const
		{
			return hint_type::template get<Size - 1>();
		}

		operator PositIntermediateFormat<N, WES, Wrapper, false>() const;
		operator PositEncoding<N, WES, Wrapper>() const;

		/*void printContent() const {
			fprintf(stderr, "isNaR: %d\n", static_cast<int>(this->getIsNaR()));

			fprintf(stderr, "biased exp: ");
			getExp().print();

			fprintf(stderr, "sign: ");
			getSignBit().print();

			fprintf(stderr, "significand: ");
			getSignificand().print();
		}*/
};

template<unsigned int N, template<unsigned int, bool> class Wrapper>
using StandardPositProd = PositProd<N, hint::Static_Val<(N>>3)>::_log2, Wrapper>;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
using PositValSizedHint = Wrapper<PositDim<N, WES>::ValSize, false>;

template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
class PositIntermediateFormat<N, WES, Wrapper, true> : public PositValSizedHint<N, WES, Wrapper>
{
	//Storage :
	// isNar Exp Sign ImplicitBit Fraction
	public:
		typedef PositValSizedHint<N, WES, Wrapper> hint_type;
		template<unsigned int W>
		using wrapper_helper = Wrapper<W, false>;
		static constexpr unsigned int Size = PositDim<N, WES>::ValSize;
		static constexpr unsigned int ExpSize = PositDim<N, WES>::WE;
		static constexpr unsigned int FractionSize = PositDim<N, WES>::WF;
		PositIntermediateFormat(
				wrapper_helper<1> isNar,
				wrapper_helper<ExpSize> exp, //Warning : biased exp
				wrapper_helper<1> sign,
				wrapper_helper<1> implicit_bit,
				wrapper_helper<FractionSize> fraction):hint_type{isNar.concatenate(exp.concatenate(sign)).concatenate(implicit_bit.concatenate(fraction)).template leftpad<Size>()}
		{
		}


		// PositIntermediateFormat(hint<Size> val):PositValSizedAPUint<N, WES>(val){}
		PositIntermediateFormat(PositEncoding<N, WES, Wrapper> val):hint_type{val}{}

		PositIntermediateFormat(hint_type val):hint_type{val}{}


		inline wrapper_helper<FractionSize + 1> getSignificand() const //Implicit bit + fractional part
		{
			return hint_type::template slice<FractionSize, 0>();
		}

		inline wrapper_helper<1> getImplicitBit() const
		{
			return hint_type::template get<FractionSize>();
		}

		inline wrapper_helper<FractionSize> getFraction() const
		{
			return hint_type::template slice<FractionSize-1, 0>();
		}

		inline wrapper_helper<1> getSignBit() const
		{
			return hint_type::template get<FractionSize + 1>();
		}

		inline wrapper_helper<ExpSize> getExp() const
		{
			return hint_type::template slice<ExpSize+FractionSize+1, FractionSize+2>();
		}

		inline wrapper_helper<1> getIsNaR() const
		{
			return hint_type::template get<FractionSize+2+ExpSize>();
		}

		inline wrapper_helper<FractionSize+2> getSignedSignificand() const
		{
			return getSignBit().concatenate(getSignificand());
		}

		inline wrapper_helper<1> isZero() const
		{
			auto isExponentNull = getExp().or_reduction().invert();
			return isExponentNull.bitwise_and(getSignBit().invert());
		}

        operator PositEncoding<N, WES, Wrapper>() const;
        operator PositIntermediateFormat<N, WES, Wrapper, false>() const;
        operator PositProd<N, WES, Wrapper>() const;

		static PositIntermediateFormat getMaxPos()
		{ //isNar Exp Sign Implicit Frac
			return PositIntermediateFormat(
					{0}, //isNar
					{2*PositDim<N, WES>::EXP_BIAS - 1}, //Biased Exp
					{0}, //sign
					{1}, //implicit bit
					{0} //fraction
				);
		}

		static PositIntermediateFormat getMinPos()
		{
			return PositIntermediateFormat(
					{0}, //isNar
					{1}, // Biased exp
					{0}, // sign
					{1}, // implicit bit
					{0} // fraction
				);
		}

		static PositIntermediateFormat getMaxNeg()
		{
			return PositIntermediateFormat(
					{0},
					{2*PositDim<N, WES>::EXP_BIAS - 2}, // Biased Exp
					{1}, //sign
					{0}, //implicit bit
					{0}  //fraction
				);
		}

		static PositIntermediateFormat getMinNeg()
		{
			return PositIntermediateFormat(
					{0}, //isNar
					{0}, // Biased Exp
					{1}, // sign
					{0}, // implicit bit
					{0} // fraction
				);
		}
};


template<unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
class PositIntermediateFormat<N, WES, Wrapper, false> : public PositValSizedHint<N, WES, Wrapper>
{
	//Storage :
	// Guard Sticky isNar Exp Sign ImplicitBit Fraction
	public:
		typedef PositValSizedHint<N, WES, Wrapper> hint_type;
		template<unsigned int W>
		using wrapper_helper = Wrapper<W, false>;
		static constexpr unsigned int Size = PositDim<N, WES>::ValSize;
		static constexpr unsigned int ExpSize = PositDim<N, WES>::WE;
		static constexpr unsigned int FractionSize = PositDim<N, WES>::WF;
		PositIntermediateFormat(
				wrapper_helper<1> isNar,
				wrapper_helper<ExpSize> exp, //Warning : biased exp
				wrapper_helper<1> sign,
				wrapper_helper<1> implicit_bit,
				wrapper_helper<FractionSize> fraction):hint_type{isNar.concatenate(exp.concatenate(sign)).concatenate(implicit_bit.concatenate(fraction)).template leftpad<Size>()}
		{
		}

		PositIntermediateFormat(
				wrapper_helper<1> guard,
				wrapper_helper<1> sticky,
				wrapper_helper<1> isNar,
				wrapper_helper<ExpSize> exp,
				wrapper_helper<1> sign,
				wrapper_helper<1> implicit_bit,
				wrapper_helper<FractionSize> fraction):
			hint_type{(guard.concatenate(sticky)).concatenate(isNar.concatenate(exp)).concatenate(sign.concatenate(implicit_bit.concatenate(fraction)))}
		{}

		// PositIntermediateFormat(hint<Size> val):PositValSizedAPUint<N, WES>(val){}
		PositIntermediateFormat(PositEncoding<N, WES, Wrapper> val):hint_type{val}{}

		PositIntermediateFormat(hint_type val):hint_type{val}{}

		inline wrapper_helper<1> getGuardBit() const
		{
			return hint_type::template get<Size-1>();
		}

		inline wrapper_helper<1> getStickyBit() const
		{
			return hint_type::template get<Size-2>();
		}


		inline wrapper_helper<FractionSize + 1> getSignificand() const //Implicit bit + fractional part
		{
			return hint_type::template slice<FractionSize, 0>();
		}

		inline wrapper_helper<1> getImplicitBit() const
		{
			return hint_type::template get<FractionSize>();
		}

		inline wrapper_helper<FractionSize> getFraction() const
		{
			return hint_type::template slice<FractionSize-1, 0>();
		}

		inline wrapper_helper<1> getSignBit() const
		{
			return hint_type::template get<FractionSize + 1>();
		}

		inline wrapper_helper<ExpSize> getExp() const
		{
			return hint_type::template slice<ExpSize+FractionSize+1, FractionSize+2>();
		}

		inline wrapper_helper<1> getIsNaR() const
		{
			return hint_type::template get<FractionSize+2+ExpSize>();
		}

		inline wrapper_helper<FractionSize+2> getSignedSignificand() const
		{
			return getSignBit().concatenate(getSignificand());
		}

		inline wrapper_helper<1> isZero() const
		{
			auto isExponentNull = getExp().or_reduction().invert();
			return isExponentNull.bitwise_and(getSignBit().invert());
		}

        operator PositEncoding<N, WES, Wrapper>() const;
		
};

template <unsigned int N, template <unsigned int, bool> class Wrapper, bool isExact>
using StandardPIF = PositIntermediateFormat<N, hint::Static_Val<(N>>3)>::_log2, Wrapper, isExact>;


template<unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
class PositEncoding : public Wrapper<N, false>
{
public:
		typedef Wrapper<N, false> hint_type;
		template<unsigned int W>
		using wrapper_helper = Wrapper<W, false>;
	PositEncoding(hint_type const & val):hint_type{val}{}

	operator PositIntermediateFormat<N, WES, Wrapper, true>() const;
	operator PositProd<N, WES, Wrapper>() const;
//    template<int NB_CARRY>
//    operator Quire<N, WES, NB_CARRY>() const;
//    template<int NB_CARRY, int banksize>
//    operator SegmentedQuire<N, WES, NB_CARRY, banksize>() const;
};

template<unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
PositEncoding<N, WES, Wrapper> operator+(
		PositEncoding<N, WES, Wrapper> const & lhs,
		PositEncoding<N, WES, Wrapper> const & rhs
   );

template<unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
PositEncoding<N, WES, Wrapper> operator-(
		PositEncoding<N, WES, Wrapper> const & lhs,
		PositEncoding<N, WES, Wrapper> const & rhs
	);

template<unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
PositEncoding<N, WES, Wrapper> operator-(
		PositEncoding<N, WES, Wrapper> const & val
	)
{
	auto one = Wrapper<1, false>{{1}}.template leftpad<N>();
	auto xor_seq = Wrapper<N, false>::generateSequence({1});
	auto twocomp = (val xor xor_seq).modularAdd(one);
	return  PositEncoding<N, WES, Wrapper>{twocomp};
};

template<unsigned int N, unsigned int WES, template <unsigned int, bool> class Wrapper>
PositProd<N, WES, Wrapper> operator*(
		PositEncoding<N, WES, Wrapper> const & lhs,
		PositEncoding<N, WES, Wrapper> const & rhs
	);

template <unsigned int N, template <unsigned int, bool> class Wrapper>
using StandardPositEncoding = PositEncoding<N, hint::Static_Val<(N>>3)>::_log2, Wrapper>;

#include "conversions.ipp"

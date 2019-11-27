#pragma once

#include "posit_dim.hpp"

using namespace std;

template <unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositProd<N, WES, Wrapper> posit_mul(PositIntermediateFormat<N, WES, Wrapper, true> in1, PositIntermediateFormat<N, WES, Wrapper, true> in2)
{
	auto isNar = in1.getIsNaR() | in2.getIsNaR();
	auto isZero = in1.isZero() | in2.isZero();
	auto mul_op1 = in1.getSignedSignificand().as_signed();
	auto mul_op2 = in2.getSignedSignificand().as_signed();
	// Compute the significand
	auto significand = (mul_op1 * mul_op2).as_unsigned();

	//     fprintf(stderr, "ICI\n");
	// significand.print();

	auto sign = significand.template get<PositDim<N, WES>::ProdSignificandSize + 1>();
	// sign.print();
	auto exbit = significand.template get<PositDim<N, WES>::ProdSignificandSize - 1>();
	auto neg_neg_2power = sign.bitwise_xor(significand.template get<PositDim<N, WES>::ProdSignificandSize>());

	auto needs_shift = exbit.bitwise_xor(sign).bitwise_or(neg_neg_2power);

	auto potential_first_bit = exbit.bitwise_or(neg_neg_2power);
	auto last_bits = significand.template slice<PositDim<N, WES>::ProdSignificandSize - 2, 0>();

	auto fin_significand = Wrapper<PositDim<N, WES>::ProdSignificandSize, false>::mux(
				needs_shift,
				potential_first_bit.concatenate(last_bits),
				last_bits.concatenate(Wrapper<1, false>{0})
			);
	/*if (needs_shift) {
		hint<1> first_bit = (significand[PositProd<N, WES>::SignificandSize - 1])
			or neg_neg_2power;
		fin_significand = first_bit.concatenate(last_bits);
	} else {
		fin_significand = last_bits.concatenate(hint<1>{0});
	}*/

	// Compute the exponent
	auto expSumVal = in1.getExp().addWithCarry(
				in2.getExp(),
				neg_neg_2power
			).template slice<PositDim<N, WES>::ProdExpSize-1, 0>().modularSub(
					needs_shift.invert().template leftpad<PositDim<N, WES>::ProdExpSize>()
			);
	auto exponent = Wrapper<PositDim<N, WES>::ProdExpSize, false>::mux(
				isZero,
				{0},
				expSumVal
			);

	return PositProd<N, WES, Wrapper>{
		isNar,
		exponent,
		sign,
		fin_significand
	};
}

template <unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositIntermediateFormat<N, WES, Wrapper, true>  posit_mul_in_place(PositIntermediateFormat<N, WES, Wrapper, true> in1, PositIntermediateFormat<N, WES, Wrapper, true> in2)
{

    // TODO
	cerr << "Use of empty function PositIntermediateFormat<N, WES, Wrapper, true>  posit_mul_in_place(PositIntermediateFormat<N, WES, Wrapper, true> in1, PositIntermediateFormat<N, WES, Wrapper, true> in2)" << endl;

    return PositIntermediateFormat<N, WES, Wrapper, true>{{0}};
}

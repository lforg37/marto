#pragma once

#include "posit_dim.hpp"

using namespace std;

#ifdef POSIT_MUL_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif


template <unsigned int N, unsigned int WES, template<unsigned int, bool> class Wrapper>
inline PositProd<N, WES, Wrapper> posit_mul(PositIntermediateFormat<N, WES, Wrapper, true> in1, PositIntermediateFormat<N, WES, Wrapper, true> in2)
{
	auto isNar = in1.getIsNaR() | in2.getIsNaR();
	auto isZero = in1.isZero() | in2.isZero();
	auto mul_op1 = in1.getSignedSignificand().as_signed();
	auto mul_op2 = in2.getSignedSignificand().as_signed();
	// Compute the significand
	auto significand = (mul_op1 * mul_op2).as_unsigned();

	auto sign = significand.template get<PositDim<N, WES>::ProdSignificandSize+1>();
	// sign.print();
	auto exbit = significand.template get<PositDim<N, WES>::ProdSignificandSize - 1>();
	auto neg_neg_2power = sign.invert() & significand.template get<PositDim<N, WES>::ProdSignificandSize>();

	auto needs_shift = exbit.bitwise_xor(sign).bitwise_or(neg_neg_2power);

	auto potential_first_bit = exbit.bitwise_or(neg_neg_2power);
	auto last_bits = significand.template slice<PositDim<N, WES>::ProdSignificandSize - 2, 0>();

	auto fin_significand = Wrapper<PositDim<N, WES>::ProdSignificandSize, false>::mux(
				needs_shift,
				potential_first_bit.concatenate(last_bits),
				last_bits.concatenate(Wrapper<1, false>{0})
			);

	// Compute the exponent
	auto exp1 = in1.getExp().as_signed();
	auto exp2 = in2.getExp().as_signed();
	auto expSumValWoShift = exp1.addWithCarry(
				exp2,
				needs_shift
			).as_unsigned();

	auto exp_seq = neg_neg_2power.template leftpad<PositDim<N, WES>::ProdExpSize>();

	auto expSumVal = expSumValWoShift.modularAdd(exp_seq);

	auto exponent = Wrapper<PositDim<N, WES>::ProdExpSize, false>::mux(
				isZero,
				{0},
				expSumVal
			);

#ifdef POSIT_MUL_DEBUG
	cerr << "=== POSIT_MUL ===" << endl;
	cerr << "isNar: " << to_string(isNar) << endl;
	cerr << "isZero: " << to_string(isZero) << endl;
	cerr << "mul_op1: " << to_string(mul_op1) << endl;
	cerr << "mul_op2: " << to_string(mul_op2) << endl;
	cerr << "significand: " << to_string(significand) << endl;
	cerr << "sign: " << to_string(sign) << endl;
	cerr << "exbit: " << to_string(exbit) << endl;
	cerr << "neg_neg_2power: " << to_string(neg_neg_2power) << endl;
	cerr << "needs_shift: " << to_string(needs_shift) << endl;
	cerr << "potential_first_bit: " << to_string(potential_first_bit) << endl;
	cerr << "last_bits: " << to_string(last_bits) << endl;
	cerr << "expSumValWoShift: " << to_string(expSumValWoShift) << endl;
	cerr << "exp_seq: " << to_string(exp_seq) << endl;
	cerr << "fin_significand: " << to_string(fin_significand) << endl;
	cerr << "expSumVal: " << to_string(expSumVal) << endl;
	cerr << "exponent: " << to_string(exponent) << endl;
	cerr << "=================" << endl;
#endif

	return PositProd<N, WES, Wrapper>{
		isNar,
		exponent,
		sign,
		fin_significand
	};
}

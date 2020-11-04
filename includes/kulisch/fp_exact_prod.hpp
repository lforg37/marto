#ifndef FP_EXACT_PROD_HPP
#define FP_EXACT_PROD_HPP
#include "ieeefloats/ieeetype.hpp"
#include "kulisch_dim.hpp"

#ifdef IEEE_MUL_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif


template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
inline FPProd<WE, WF, Wrapper> exact_prod(
		IEEENumber<WE, WF, Wrapper> const & in1,
		IEEENumber<WE, WF, Wrapper> const & in2
){

	using _dim = FPDim<WE, WF>;

	auto m1 = in1.getLeadBitVal().concatenate(in1.getFractionnalPart());
	auto m2 = in2.getLeadBitVal().concatenate(in2.getFractionnalPart());

	auto m1_is_zero = m1.nor_reduction();
	auto m2_is_zero = m2.nor_reduction();

	auto e1 = in1.getExponent();
	auto e2 = in2.getExponent();

	auto e1_is_n = e1.or_reduction();
	auto e2_is_n = e2.or_reduction();

	auto e1_is_f_o = e1.and_reduction();
	auto e2_is_f_o = e2.and_reduction();

	auto in1_is_zero = e1_is_n.invert() & m1_is_zero;
	auto in2_is_zero = e2_is_n.invert() & m2_is_zero;

	auto in1_is_inf = e1_is_f_o & m1_is_zero;
	auto in2_is_inf = e2_is_f_o & m2_is_zero;

	auto in1_is_nan = e1_is_f_o & m1_is_zero.invert();
	auto in2_is_nan = e2_is_f_o & m2_is_zero.invert();

	auto one_is_nan = in1_is_nan | in2_is_nan;
	auto one_is_inf = in1_is_inf | in2_is_inf;
	auto one_iszero = in1_is_zero | in2_is_zero;

	auto res_is_nan = one_is_nan | (one_iszero & one_is_inf);
	auto res_is_inf = one_is_inf & one_is_nan.invert() & one_iszero.invert();
	auto res_is_zero = one_iszero & one_is_nan.invert()& one_is_inf.invert();
	auto res_is_normal = one_is_inf.invert() & one_is_nan.invert() & one_iszero.invert();

	auto exp_mask = Wrapper<_dim::WE_Prod, false>::generateSequence(res_is_normal);
	auto m_mask = Wrapper<_dim::WFF_Prod, false>::generateSequence(res_is_normal);


	auto sub_two = e1_is_n & e2_is_n;
	auto sub_one = e1_is_n ^ e2_is_n;

	auto to_sub = sub_two.concatenate(sub_one).template leftpad<WE+1>();

	auto s1 = in1.getSign();
	auto s2 = in2.getSign();

	auto mult_s = (s1 ^ s2) & res_is_zero.invert();
	auto mult_e = ((e1 + e2).modularSub(to_sub)) & exp_mask;
	auto mult_m = (m1 * m2) & m_mask;

#ifdef IEEE_MUL_DEBUG
	cerr << "=== FP_EXACT_PROD ===" << endl;
	cerr << "m1: " << to_string(m1) << endl;
	cerr << "m2: " << to_string(m2) << endl;
	cerr << "e1: " << to_string(e1) << endl;
	cerr << "e2: " << to_string(e2) << endl;

	cerr << "e1_is_n: " << to_string(e1_is_n) << endl;
	cerr << "e2_is_n: " << to_string(e2_is_n) << endl;
	cerr << "sub_two: " << to_string(sub_two) << endl;
	cerr << "sub_one: " << to_string(sub_one) << endl;
	cerr << "to_sub: " << to_string(to_sub) << endl;

	cerr << "s1: " << to_string(s1) << endl;
	cerr << "s2: " << to_string(s2) << endl;
	cerr << "mult_s: " << to_string(mult_s) << endl;
	cerr << "mult_e: " << to_string(mult_e) << endl;
	cerr << "mult_m: " << to_string(mult_m) << endl;
	cerr << "=====================" << endl;
#endif


	return FPProd<WE, WF, Wrapper>(
			mult_s,
			mult_e,
			mult_m,
			res_is_nan,
			res_is_inf
		);
}
#endif

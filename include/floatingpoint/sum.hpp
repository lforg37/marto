#ifndef SUM_HPP
#define SUM_HPP

#include <primitives/lzoc_shifter.hpp>
#include <primitives/shifter.hpp>
#include <tools/static_math.hpp>

using hint::LZOC_shift;
using hint::Static_Val;
using hint::shifter;

#ifdef FPEXPR_SUM_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif


#include "fp_number.hpp"

template<typename Dim1, typename Dim2>
struct FPSumDim
{
	private:
		static constexpr int64_t op1maxE = Dim1::MAX_EXP;
		static constexpr int64_t op2maxE = Dim2::MAX_EXP;
		//Overestimation
		static constexpr int64_t maxExp = (op1maxE > op2maxE) ? op1maxE + 1 : op2maxE + 1;

		static constexpr int64_t op1minE = Dim1::MIN_EXP;
		static constexpr int64_t op2minE = Dim2::MIN_EXP;

		static constexpr int64_t op1WF{Dim1::WF};
		static constexpr int64_t op2WF{Dim2::WF};

		//Overestimation
		static constexpr int64_t minMinExp = (op1minE > op2minE) ? op2minE : op1minE;
		static constexpr int64_t op1lowbitExp = op1minE - op1WF;
		static constexpr int64_t op2lowbitExp = op2minE - op2WF;
		static constexpr bool possibleCancellation = (op1lowbitExp <= minMinExp) and (op2lowbitExp <= minMinExp);
		static constexpr int64_t maxLowBitExp = (op1lowbitExp > op2lowbitExp) ? op1lowbitExp : op2lowbitExp;
		static constexpr int64_t minExp = (possibleCancellation) ? maxLowBitExp : minMinExp;

		static constexpr vec_width range1 = op1maxE - op2minE + 1;
		static constexpr vec_width range2 = op2maxE - op1minE + 1;
		static constexpr vec_width WFmax = (range1 > range2) ? range1 : range2;

	public:
		using type = TightFixedFormat<WFmax, maxExp, minExp>;
};

template<vec_width TargetWF, typename Dim1, typename Dim2>
struct RoundedFPSum
{
	private:
		using ExactSumDim = typename FPSumDim<Dim1, Dim2>::type;
		static constexpr vec_width max_wf = (Dim1::WF > Dim2::WF) ? Dim1::WF : Dim2::WF;
		static constexpr vec_width full_sum_wf = ExactSumDim::WF;
		static constexpr bool isExact = (TargetWF >= full_sum_wf);
		static constexpr vec_width TruncatingWidth = isExact ?  full_sum_wf : TargetWF;
		static constexpr vec_width roundingbitOffset = (isExact) ? 0 : 1;
		static constexpr vec_width partialWF = TruncatingWidth + roundingbitOffset;
		static constexpr vec_width max_signif_width = max_wf + 1;
		//static_assert (not isExact, "exact  sum is buggy at the moment");

		using partial_sum_dim = FixedFormat<ExactSumDim::WE, partialWF, ExactSumDim::MAX_EXP, ExactSumDim::MIN_EXP>;
		static constexpr vec_width partial_we = partial_sum_dim::WE;
		using rounding_helper = RoundDimHelper<partial_sum_dim, TargetWF>;
		static constexpr vec_width sum_op_align = (TruncatingWidth > max_wf) ? TruncatingWidth : max_wf; //take at minimum maxWF bits
		static constexpr vec_width op_align_wf = sum_op_align + roundingbitOffset;
		static constexpr vec_width significand_align_width = 1 + op_align_wf;
		static constexpr vec_width lzoc_count = (Static_Val<max_signif_width + 1>::_2pow == 1 + max_signif_width) ? Static_Val<1 + max_signif_width>::_2pow : Static_Val<1 + max_signif_width>::_2pow - 1;
		static constexpr vec_width lzoc_size = Static_Val<lzoc_count>::_storage;

		template<template<unsigned int, bool> class Wrapper>
		static inline FixedNumber<partial_sum_dim, Wrapper> compute_truncated_sum(
				FixedNumber<Dim1, Wrapper> const & op1,
				FixedNumber<Dim2, Wrapper> const & op2
			)
		{
			auto op1IsZero = op1.isZero();
			auto op2IsZero = op2.isZero();

			auto isInf1 = op1.isInf();
			auto isInf2 = op2.isInf();

			auto exp1 = op1.getExponent().template sign_extend<partial_we>();
			auto exp2 = op2.getExponent().template sign_extend<partial_we>();

			auto diff1 = exp1.modularSub(exp2);
			auto diff2 = exp2.modularSub(exp1);

			auto signif1 = op1.getSignificand().template rightpad<max_signif_width>();
			auto signif2 = op2.getSignificand().template rightpad<max_signif_width>();

			auto negSignif1 = op1.getSignificand().invert().modularAdd({{1}}).template rightpad<max_signif_width>();
			auto negSignif2 = op2.getSignificand().invert().modularAdd({{1}}).template rightpad<max_signif_width>();

			auto exp1_greater = exp1 > exp2;
			auto exp_equals = (exp1 == exp2);

			auto frac1_greater = (signif1 > signif2);

			auto op1_greater = isInf2.invert() & (isInf1 | (op1IsZero.invert() & (exp1_greater | (exp_equals & frac1_greater) | op2IsZero)));
			auto exp_diff = Wrapper<partial_we, false>::mux(op1_greater, diff1, diff2);
			auto exp_diff_is_one = (exp_diff == Wrapper<partial_we, false>{{1}});

			auto s1 = op1.getSign();
			auto s2 = op2.getSign();

			auto eff_sub = s1 ^ s2;

			auto g_exp = Wrapper<partial_we, true>::mux(op1_greater, exp1, exp2);

			auto s_res = (s1 & op1_greater) | (s2 & op1_greater.invert());


			auto g_signif = Wrapper<max_signif_width, false>::mux(op1_greater, signif1, signif2);
			auto l_signif_pos = Wrapper<max_signif_width, false>::mux(op1_greater, signif2, signif1);
			auto l_signif_neg = Wrapper<max_signif_width, false>::mux(op1_greater, negSignif2, negSignif1);
			auto l_signif = Wrapper<max_signif_width, false>::mux(eff_sub, l_signif_neg, l_signif_pos);

			//------ Flags ----------------------------------------------------------------------------------------------------/

			auto bothInf = isInf1 & isInf2;
			auto isNaN = op1.isNaN() | op2.isNaN() | (bothInf & eff_sub);
			auto isInf = (isInf1 | isInf2) & isNaN.invert();
			auto isZeroFromFlags = op1IsZero & op2IsZero;

			//------ Close path ------------------------------------------------------------------------------------------------
			// Handle only the cancellation case
			// Two fracs are extended on signif_align_width + 1
			auto g_signif_close = g_signif.template rightpad<max_signif_width + 1>();
			auto p_l_close = l_signif.template rightpad<max_signif_width + 1>();
			auto l_close = (Wrapper<1, false>{{1}}.concatenate(p_l_close) >> (exp_equals.invert())).template slice <max_signif_width, 0>();
			auto close_add = g_signif_close.modularAdd(l_close);
			auto lzoc_shifted = LZOC_shift<max_signif_width + 1, lzoc_count>(close_add, {{0}});
			auto lzoc = lzoc_shifted.lzoc;
			auto cancellation_full_frac = lzoc_shifted.shifted.template slice<max_signif_width, 0>();
			auto cancellation_exp = g_exp.modularSub(lzoc.template leftpad<partial_we>().as_signed()).as_signed();
			auto cancel_to_zero = lzoc_shifted.shifted.template get<max_signif_width>().invert();
			auto up_cancellation_frac = cancellation_full_frac.template slice<max_signif_width-1, 0>();
			auto cancellation_frac = up_cancellation_frac.template rightpad<op_align_wf + 1>().template slice<op_align_wf, op_align_wf - partialWF + 1>();

			//------ Far Path ---------------------------------------------------------------------------------------------------
			auto normal_g_signif = g_signif.template rightpad<significand_align_width + 1>();
			auto shifted_low_frac = shifter<true>(l_signif.template rightpad<significand_align_width + 1>(), exp_diff, eff_sub);
			auto add_res = normal_g_signif.addWithCarry(shifted_low_frac, {{0}});
			auto overflowed = add_res.template get<significand_align_width + 1>() & eff_sub.invert();
			auto normal_pos = add_res.template get<significand_align_width>() & overflowed.invert();
			auto shiftval = overflowed.concatenate(normal_pos);
			auto frac_ext = add_res.template slice<significand_align_width, 0>() >> shiftval;
			auto fracRes = frac_ext.template slice<partialWF - 1, 0>();

			auto expMask = Wrapper<partial_we, false>::generateSequence(overflowed.invert()).as_signed();
			auto expNormal = g_exp.addWithCarry(expMask, overflowed | normal_pos).template slice<partial_we - 1, false>().as_signed();
			//-------------------------------------------------------------------------------------------------------------------

			auto select_close_path = (exp_equals | exp_diff_is_one) & eff_sub;
			auto final_frac = Wrapper<partialWF, false>::mux(select_close_path, cancellation_frac, fracRes);
			auto final_exp = Wrapper<partial_we, true>::mux(select_close_path, cancellation_exp, expNormal);

			auto isZero = isZeroFromFlags | (select_close_path & cancel_to_zero);

#ifdef FPEXPR_SUM_DEBUG
			cerr << "==== FPEXPR SUM ======" << endl;
			cerr << "op1IsZero: " << to_string(op1IsZero) << endl;
			cerr << "op2IsZero: " << to_string(op2IsZero) << endl;
			cerr << "isInf1: " << to_string(isInf1) << endl;
			cerr << "isInf2: " << to_string(isInf2) << endl;
			cerr << "exp1: " << to_string(exp1) << endl;
			cerr << "exp2: " << to_string(exp2) << endl;
			cerr << "diff1: " << to_string(diff1) << endl;
			cerr << "diff2: " << to_string(diff2) << endl;
			cerr << "signif1: " << to_string(signif1) << endl;
			cerr << "signif2: " << to_string(signif2) << endl;
			cerr << "negSignif1: " << to_string(negSignif1) << endl;
			cerr << "negSignif2: " << to_string(negSignif2) << endl;
			cerr << "exp1_greater: " << to_string(exp1_greater) << endl;
			cerr << "exp_equals: " << to_string(exp_equals) << endl;
			cerr << "frac1_greater: " << to_string(frac1_greater) << endl;
			cerr << "op1_greater: " << to_string(op1_greater) << endl;
			cerr << "exp_diff: " << to_string(exp_diff) << endl;
			cerr << "exp_diff_is_one: " << to_string(exp_diff_is_one) << endl;
			cerr << "s1: " << to_string(s1) << endl;
			cerr << "s2: " << to_string(s2) << endl;
			cerr << "eff_sub: " << to_string(eff_sub) << endl;
			cerr << "g_exp: " << to_string(g_exp) << endl;
			cerr << "s_res: " << to_string(s_res) << endl;
			cerr << "g_signif: " << to_string(g_signif) << endl;
			cerr << "l_signif_pos: " << to_string(l_signif_pos) << endl;
			cerr << "l_signif_neg: " << to_string(l_signif_neg) << endl;
			cerr << "l_signif: " << to_string(l_signif) << endl;
			cerr << "bothInf: " << to_string(bothInf) << endl;
			cerr << "isNaN: " << to_string(isNaN) << endl;
			cerr << "isInf: " << to_string(isInf) << endl;
			cerr << "isZeroFromFlags: " << to_string(isZeroFromFlags) << endl;
			cerr << "g_signif_close: " << to_string(g_signif_close) << endl;
			cerr << "p_l_close: " << to_string(p_l_close) << endl;
			cerr << "l_close: " << to_string(l_close) << endl;
			cerr << "close_add: " << to_string(close_add) << endl;
			cerr << "lzoc: " << to_string(lzoc) << endl;
			cerr << "cancellation_full_frac: " << to_string(cancellation_full_frac) << endl;
			cerr << "cancellation_exp: " << to_string(cancellation_exp) << endl;
			cerr << "cancel_to_zero: " << to_string(cancel_to_zero) << endl;
			cerr << "cancellation_frac: " << to_string(cancellation_frac) << endl;
			cerr << "normal_g_signif: " << to_string(normal_g_signif) << endl;
			cerr << "shifted_low_frac: " << to_string(shifted_low_frac) << endl;
			cerr << "add_res: " << to_string(add_res) << endl;
			cerr << "overflowed: " << to_string(overflowed) << endl;
			cerr << "normal_pos: " << to_string(normal_pos) << endl;
			cerr << "shiftval: " << to_string(shiftval) << endl;
			cerr << "frac_ext: " << to_string(frac_ext) << endl;
			cerr << "fracRes: " << to_string(fracRes) << endl;
			cerr << "expMask: " << to_string(expMask) << endl;
			cerr << "expNormal: " << to_string(expNormal) << endl;
			cerr << "select_close_path: " << to_string(select_close_path) << endl;
			cerr << "final_frac: " << to_string(final_frac) << endl;
			cerr << "final_exp: " << to_string(final_exp) << endl;
			cerr << "isZero: " << to_string(isZero) << endl;
			cerr << "============" << endl;
#endif

			return FixedNumber<partial_sum_dim, Wrapper>(
						final_frac,
						final_exp,
						s_res,
						isInf,
						isNaN,
						isZero
						);
		}

	public:
		using dim = typename rounding_helper::dim;
		template<template<unsigned int, bool> class Wrapper>
		static FixedNumber<dim, Wrapper> compute(FixedNumber<Dim1, Wrapper> const & op1, FixedNumber<Dim2, Wrapper> const & op2)
		{
			return Rounder<dim, partial_sum_dim>::compute(compute_truncated_sum(op1, op2));
		}
};

#endif // SUM_HPP

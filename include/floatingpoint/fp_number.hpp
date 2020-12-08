#ifndef FP_DIM_HPP
#define FP_DIM_HPP

#include <cstdint>
#include <type_traits>
#include <utility>

using std::conditional;
using std::enable_if;
using std::pair;

#include <tools/static_math.hpp>

using hint::Static_Val;

typedef unsigned int vec_width;
typedef int biasval;


#ifdef FPEXPR_ROUND_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif


template<int64_t maxexp, int64_t minexp>
struct TightWEHelper
{
		static_assert (maxexp > 0, "getTighhtWE expects maxexp to be strictly positive");
		static_assert (minexp < 0, "getTighhtWE expects minexp to be strictly negative");
		static constexpr int64_t absMinExp = -minexp;
		static constexpr bool minExprIsToPow = Static_Val<absMinExp>::_is2Pow;
		static constexpr vec_width minExpStorage = Static_Val<absMinExp>::_storage;
		static constexpr vec_width minExpRequiredWidth = (minExprIsToPow) ? minExpStorage : minExpStorage + 1;
		static constexpr vec_width maxExpStorage = Static_Val<maxexp>::_storage;
		static constexpr vec_width maxExpRequiredWidth = maxExpStorage + 1;
		static constexpr vec_width requiredWidth = (maxExpRequiredWidth > minExpRequiredWidth) ? maxExpRequiredWidth : minExpRequiredWidth;
};

template <vec_width _WE, vec_width _WF, int64_t _maxExp = 0, int64_t _minExp = 0>
struct FPDim {
		/**
		 * @brief WE Exponent field width
		 */
		static constexpr vec_width WE = _WE;

		/**
		 * @brief WF Fractional width
		 */
		static constexpr vec_width WF = _WF;

		/**
		 * @brief WFF Full significand width
		 */
		static constexpr vec_width WFS = WF +1;

		/**
		  * @brief MAX_STORABLE_EXP maximum positiv exponent storable on WE bits
		  */
		static constexpr int64_t MAX_STORABLE_EXP = (int64_t{1} << (WE - 1)) - 1;

		/**
		 * @brief MIN_STORABLE_EXP minimal exponant storable on WE bits
		 */
		static constexpr int64_t MIN_STORABLE_EXP = -MAX_STORABLE_EXP - 1;

		/**
		 * @brief EXP_CONSTRAINED true if there is a constraint on exp values
		 */
		static constexpr bool EXP_CONSTRAINED = (_maxExp != _minExp) || (_maxExp != 0);

		static_assert(not(EXP_CONSTRAINED and (_maxExp > MAX_STORABLE_EXP)), "Max exponent constraint is greater than max storable exp");
		static_assert(not(EXP_CONSTRAINED and (_minExp < MIN_STORABLE_EXP)), "Min exponent constraint is smaller than min storable exp");
		static_assert(not(EXP_CONSTRAINED and (_maxExp <= _minExp)), "Max exponent constraint smaller than min exp xonstraint");

		static constexpr int64_t MIN_EXP = (EXP_CONSTRAINED) ? (_minExp) : (MIN_STORABLE_EXP);
		static constexpr int64_t MAX_EXP = (EXP_CONSTRAINED) ? (_maxExp) : (MAX_STORABLE_EXP);
};

template<vec_width WF, int64_t maxExp, int64_t minExp>
using TightFPDim = FPDim<TightWEHelper<maxExp, minExp>::requiredWidth, WF, maxExp, minExp>;

template<typename sourceDim, vec_width targetWF>
struct RoundDimHelper {
		static constexpr bool CAN_ROUND = (targetWF < sourceDim::WF);
		static constexpr int64_t MAX_EXP = (CAN_ROUND) ? sourceDim::MAX_EXP + 1 : sourceDim::MAX_EXP;
		static constexpr vec_width WF = (CAN_ROUND) ? targetWF : sourceDim::WF;
		using dim = TightFPDim<WF, MAX_EXP, sourceDim::MIN_EXP>;
};

template <typename _dim, template<unsigned int, bool> class Wrapper>
class FPNumber {
	public:
		using dim = _dim;
	private:
		using sign_t = Wrapper<1, false>;
		using frac_t = Wrapper<dim::WF, false>;
		using significand_t = Wrapper<dim::WF + 1, false>;
		using exp_t = Wrapper<dim::WE, true>;
		using flag_t = Wrapper<1, false>;
		using normal_value_t = Wrapper<1 + dim::WE + dim::WF, false>;

		frac_t fraction;
		exp_t exponent;
		sign_t sign;
		flag_t NaN;
		flag_t Inf;
		flag_t nonZero;
	public:
		FPNumber(normal_value_t const & val):NaN{{0}}, Inf{{0}}, nonZero{{1}}
		{
			fraction = val.template slice<dim::WF-1, 0>();
			exponent = val.template slice<dim::WF + dim::WE -1, dim::WF>();
			sign = val.template get<dim::WF+dim::WE>();
		}

		FPNumber(
				frac_t const & frac,
				exp_t const & exp,
				sign_t const & s,
				flag_t const & isInf,
				flag_t const & isNaN,
				flag_t const & isZero):
			fraction{frac},
			exponent{exp},
			sign{s},
			Inf{isInf},
			NaN{isNaN}
		{
			nonZero = isZero.invert();
		}

		inline FPNumber opposite() const {
			return {fraction, exponent, sign.invert(), NaN, Inf, nonZero};
		}

		inline frac_t getFraction() const {return fraction;}
		inline sign_t getSign() const {return sign;}
		inline exp_t getExponent() const {return exponent;}
		inline flag_t isZero() const {return nonZero.invert();}
		inline flag_t isInf() const {return Inf;}
		inline flag_t isNaN() const {return NaN;}
		inline significand_t getSignificand() const {return nonZero.concatenate(fraction);}

		static inline FPNumber getZero(sign_t s = {{0}}) {
			return FPNumber{{{0}}, {{0}}, s, {{0}}, {{0}}, {{1}}};
		}

		static inline FPNumber getInf(sign_t s = {{0}}) {
			return FPNumber{{{0}}, {{0}}, s, {{1}}, {{0}}, {{0}}};
		}

		static inline FPNumber getNaN() {
			return FPNumber{{{0}}, {{0}}, {{0}}, {{0}}, {{1}}, {{0}}};
		}
};

template<typename dim, template<unsigned int, bool> class Wrapper>
Wrapper<1, false> operator==(FPNumber<dim, Wrapper> const & op1, FPNumber<dim, Wrapper> const & op2)
{
	auto exp_frac_ok = (op1.getFraction() == op2.getFraction()) & (op1.getExponent() == op2.getExponent());
	auto s_ok = (op1.getSign() == op2.getSign());
	auto spec_ok = (op1.isNaN() == op2.isNaN()) & (op1.isInf() == op2.isInf()) & (op1.isZero() == op2.isZero());
	auto res = (exp_frac_ok | op1.isInf() | op1.isNaN() | op1.isZero()) & s_ok & spec_ok;
	return res;
}

template <typename TargetDim, typename SourceDim, bool AllowPrecisionExtension = false>
struct Rounder {
	public:
		static_assert(AllowPrecisionExtension or (TargetDim::WF <= SourceDim::WF), "Rounding with precision extension");
		static constexpr bool CAN_ROUND = TargetDim::WF < SourceDim::WF;
		static constexpr int64_t max_rounded_exp = SourceDim::MAX_EXP + ((CAN_ROUND) ? 1 : 0);
		static constexpr bool EXP_CAN_OVERFLOW = (max_rounded_exp > TargetDim::MAX_STORABLE_EXP) or SourceDim::MIN_EXP < TargetDim::MIN_STORABLE_EXP;
		static constexpr vec_width overflowed_exp_width = (max_rounded_exp > SourceDim::MAX_STORABLE_EXP) ? SourceDim::WE + 1 : SourceDim::WE;
	private:
		//------------- Fraction rounding ----------------------------------------------------------------
		// Two cases that are splitted into four due to problem with zero length bit vector :
		// target WF >= source WF (no rounding) or
		// target WF < source WF
		template<vec_width target_wf, vec_width source_wf, template<unsigned int, bool> class Wrapper>
		static inline pair<Wrapper<TargetDim::WF, false>, Wrapper<1, false>> fp_round_frac(
				FPNumber<SourceDim, Wrapper> const & source,
				typename enable_if<(target_wf > source_wf)>::type* = 0
				)
		{
			return {source.getFraction().concatenate(Wrapper<TargetDim::WF-SourceDim::WF, false>{0}), {{0}}};
		}

		template<vec_width target_wf, vec_width source_wf, template<unsigned int, bool> class Wrapper>
		static inline pair<Wrapper<TargetDim::WF, false>, Wrapper<1, false>> fp_round_frac(
				FPNumber<SourceDim, Wrapper> const & source,
				typename enable_if<(target_wf == source_wf)>::type* = 0
				)
		{
			return {source.getFraction(), {{0}}};
		}

		template<vec_width target_wf, vec_width source_wf, template<unsigned int, bool> class Wrapper>
		static inline pair<Wrapper<TargetDim::WF, false>, Wrapper<1, false>> fp_round_frac(
				FPNumber<SourceDim, Wrapper> const & source,
				typename enable_if<(target_wf == (source_wf - 1))>::type* = 0
				)
		{
			auto frac = source.getFraction();
			auto rounded = frac.addWithCarry({{1}}, {{0}});
			auto overflow = rounded.template get<SourceDim::WF>();
			auto resFrac = rounded.template slice<SourceDim::WF-1, 1>();
			return {resFrac, overflow};
		}

		template<vec_width target_wf, vec_width source_wf, template<unsigned int, bool> class Wrapper>
		static inline pair<Wrapper<TargetDim::WF, false>, Wrapper<1, false>> fp_round_frac(
				FPNumber<SourceDim, Wrapper> const & source,
				typename enable_if<(target_wf < source_wf - 1)>::type* = 0
				)
		{
			auto frac = source.getFraction();
			auto top = frac.template slice<SourceDim::WF-1, SourceDim::WF-TargetDim::WF - 1>();
			auto rounded = top.addWithCarry({{1}}, {{0}});
			auto overflow = rounded.template get<TargetDim::WF+1>();
			auto resFrac = rounded.template slice<TargetDim::WF, 1>();
			return {resFrac, overflow};
		}

		//----------------------------------------------------------------End Fraction Rounding -------
		//---------------------- Exponent Rounding ----------------------------------------------------
		template <bool exp_can_overflow, template<unsigned int, bool> class Wrapper>
		static inline pair<Wrapper<TargetDim::WE, true>, Wrapper<2, false>> fp_round_exp(
				Wrapper<overflowed_exp_width, true> const & source_exp,
				typename enable_if<exp_can_overflow>::type* = 0)
		{
			Wrapper<overflowed_exp_width, true> maxReprExp{TargetDim::MAX_STORABLE_EXP};
			Wrapper<overflowed_exp_width, true> minReprExp{TargetDim::MIN_STORABLE_EXP};
			Wrapper<overflowed_exp_width, true> minReprExpLim{TargetDim::MIN_STORABLE_EXP - 1};
			Wrapper<TargetDim::WE, true> maxReprRes{TargetDim::MAX_STORABLE_EXP};
			Wrapper<TargetDim::WE, true> minReprRes{TargetDim::MIN_STORABLE_EXP};
			auto exp_sign = source_exp.template get<overflowed_exp_width - 1>();
			auto exp_greater = source_exp > maxReprExp;
			auto exp_smaller = source_exp < minReprExp;
			auto ov_und_flow = exp_greater | exp_smaller;
			auto exp_is_limit = (source_exp == minReprExpLim);
			auto overflow_res = Wrapper<TargetDim::WE, true>::mux(exp_sign, minReprRes, maxReprRes);
			auto normal_res = source_exp.template slice<TargetDim::WE-1, 0>().as_signed();
			auto res = Wrapper<TargetDim::WE, true>::mux(ov_und_flow, overflow_res, normal_res);
#ifdef FPEXPR_ROUND_DEBUG
			cerr << "--------> FPRoundExp<true> " << endl;
			cerr << "exp_sign: " << to_string(exp_sign) << endl;
			cerr << "exp_greater: " << to_string(exp_greater) << endl;
			cerr << "exp_smaller: " << to_string(exp_smaller) << endl;
			cerr << "ov_und_flow: " << to_string(ov_und_flow) << endl;
			cerr << "exp_is_limit: " << to_string(exp_is_limit) << endl;
			cerr << "overflow_res: " << to_string(overflow_res) << endl;
			cerr << "normal_res: " << to_string(normal_res) << endl;
			cerr << "res: " << to_string(res) << endl;
#endif
			return {res, ov_und_flow.concatenate(exp_is_limit)};
		}

		template <bool exp_can_overflow, template<unsigned int, bool> class Wrapper>
		static inline pair<Wrapper<TargetDim::WE, true>, Wrapper<2, false>> fp_round_exp(
				Wrapper<overflowed_exp_width, true> const & source_exp,
				typename enable_if<not exp_can_overflow>::type* = 0)
		{
			constexpr vec_width max_we = (TargetDim::WE > overflowed_exp_width) ? TargetDim::WE : overflowed_exp_width;
			auto ext_exp = source_exp.template sign_extend<max_we>().template slice<TargetDim::WE - 1, 0>().as_signed();
			return {ext_exp, {{0}}};
		}
		//---------------------------------------------------------------- End exponent rounding

		//---------------------- Exponent Overflow -----------------------
		template<bool can_round, template<unsigned int, bool> class Wrapper>
		static inline pair<Wrapper<overflowed_exp_width, true>, Wrapper<overflowed_exp_width, true>> getExpPossibilities(
				Wrapper<SourceDim::WE, true> const & inExp,
				typename enable_if<can_round>::type* = 0
				) {
			auto extExp = inExp.template sign_extend<overflowed_exp_width>();
			auto inc = extExp.modularAdd({{1}}).as_signed();
			return {extExp, inc};
		}

		template<bool can_round, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<overflowed_exp_width, true> getExpPossibilities(
				Wrapper<SourceDim::WE, true> const & inExp,
				typename enable_if<not can_round>::type* = 0
				) {
			return inExp;
		}

		template<bool can_round, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<overflowed_exp_width, true> selectExp(
				pair<Wrapper<overflowed_exp_width, true>, Wrapper<overflowed_exp_width, true>> const & choices,
				Wrapper<1, false> const & overflowed,
				typename enable_if<can_round>::type* = 0
				)
		{
			return Wrapper<overflowed_exp_width, true>::mux(overflowed, choices.second, choices.first);
		}

		template<bool can_round, template<unsigned int, bool> class Wrapper>
		static inline Wrapper<overflowed_exp_width, true> selectExp(
				Wrapper<overflowed_exp_width, true> const & choices,
				Wrapper<1, false> const & ,
				typename enable_if<not can_round>::type* = 0
				)
		{
			return choices;
		}
		//---------------------------------------------------------------------------------

	public:
		template<template<unsigned int, bool> class Wrapper>
		static inline FPNumber<TargetDim, Wrapper>
		compute(FPNumber<SourceDim, Wrapper> const & source)
		{
			auto roundFrac_ov = fp_round_frac<TargetDim::WF, SourceDim::WF>(source);
			auto roundedFrac = roundFrac_ov.first;
			auto fracOverflowed = roundFrac_ov.second;
			auto exp = source.getExponent();
			auto choices = getExpPossibilities<CAN_ROUND>(exp);
			auto ov_exp = selectExp<CAN_ROUND>(choices, fracOverflowed);
			auto roundExpInf = fp_round_exp<EXP_CAN_OVERFLOW>(ov_exp);
			auto ov_und_exp = roundExpInf.second;
			auto should_clear_frac = ov_und_exp.template get<1>();
			auto ov_und = ov_und_exp.template get<1>();
			auto erase = ov_und_exp.template get<0>() & fracOverflowed.invert();
			auto rounded_exp = roundExpInf.first;
			auto exp_sign = rounded_exp.template get<TargetDim::WE-1>();
			auto isNaN = source.isNaN();
			auto isZero = isNaN.invert() & (source.isZero() | (ov_und & exp_sign & erase.invert()));
			auto isInf = isNaN.invert() & (source.isInf() | (ov_und & exp_sign.invert()));
			auto finalFrac = roundedFrac & Wrapper<TargetDim::WF, false>::generateSequence(erase.invert());

#ifdef FPEXPR_ROUND_DEBUG
			cerr << "====== FPEXPR_ROUND ======" << endl;
			cerr << "roundedFrac: " << to_string(roundedFrac) << endl;
			cerr << "exp: " << to_string(exp) << endl;
			cerr << "ov_exp: " << to_string(ov_exp) << endl;
			cerr << "ov_und_exp: " << to_string(ov_und_exp) << endl;
			cerr << "should_clear_frac: " << to_string(should_clear_frac) << endl;
			cerr << "ov_und: " << to_string(ov_und) << endl;
			cerr << "erase: " << to_string(erase) << endl;
			cerr << "rounded_exp: " << to_string(rounded_exp) << endl;
			cerr << "exp_sign: " << to_string(exp_sign) << endl;
			cerr << "isNaN: " << to_string(isNaN) << endl;
			cerr << "isZero: " << to_string(isZero) << endl;
			cerr << "isInf: " << to_string(isInf) << endl;
			cerr << "finalFrac: " << to_string(finalFrac) << endl;
			cerr << "============" << endl;
#endif

			return {finalFrac, rounded_exp, source.getSign(), isInf, isNaN, isZero};
		}
};

template<vec_width targetWF, typename sourceDim>
struct strictRounderOp {
	private:
		using round_helper = RoundDimHelper<sourceDim, targetWF>;
	public:
		using dim = typename round_helper::dim;
		template<template<unsigned int, bool> class Wrapper>
		static inline FPNumber<dim, Wrapper>
		compute(FPNumber<sourceDim, Wrapper> const & source)
		{
			return Rounder<dim, sourceDim>::compute(source);
		}
};

template<typename sourceDim>
struct TightResize {
	private:
		using tightdim = TightFPDim<sourceDim::WF, sourceDim::MAX_EXP, sourceDim::MIN_EXP>;
		static_assert(tightdim::WE <= sourceDim::WE, "TightResize should only be used with complete exponent range fp dim");
		static constexpr bool isReduce = (tightdim::WE < sourceDim::WE);
	public:
		using dim = typename conditional<isReduce, tightdim, sourceDim>::type;
	private:
		template<bool reduction, template<unsigned int, bool> class Wrapper>
		static inline FPNumber<dim, Wrapper> do_compute(
				FPNumber<sourceDim, Wrapper> const & in,
				typename enable_if<reduction>::type* = 0
		) {
			auto exp = in.getExponent().template slice<dim::WE-1, 0>().as_signed();
			return {
				in.getFraction(),
				exp,
				in.getSign(),
				in.isInf(),
				in.isNaN(),
				in.isZero()
			};
		}

		template<bool reduction, template<unsigned int, bool> class Wrapper>
		static inline FPNumber<dim, Wrapper> do_compute(
				FPNumber<sourceDim, Wrapper> const & in,
				typename enable_if<not reduction>::type* = 0
		) {
			return in;
		}
	public:
		template<template<unsigned int, bool> class Wrapper>
		static inline FPNumber<dim, Wrapper>
		compute(FPNumber<sourceDim, Wrapper> const & source)
		{
			return do_compute<isReduce>(source);
		}
};

template<vec_width, typename sourceDim>
struct OppositeOp {
	public:
		using dim = sourceDim;
		template<template<unsigned int, bool> class Wrapper>
		static inline FPNumber<dim, Wrapper>
		compute(FPNumber<sourceDim, Wrapper> const & source)
		{
			return FPNumber<sourceDim, Wrapper>(
						source.getFraction(),
						source.getExponent(),
						source.getSign().invert(),
						source.isInf(),
						source.isNaN(),
						source.isZero()
						);
		}
};

template<vec_width, typename sourceDim>
struct IdentityOp {
	public:
		using dim = sourceDim;
		template<template<unsigned int, bool> class Wrapper>
		static inline FPNumber<dim, Wrapper>
		compute(FPNumber<sourceDim, Wrapper> const & source)
		{
			return source;
		}
};
#endif

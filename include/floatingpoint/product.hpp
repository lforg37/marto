#ifndef PRODUCT_HPP
#define PRODUCT_HPP

#include <type_traits>

#include "fp_number.hpp"

using std::enable_if;

template<typename DIM1, typename DIM2>
using FPProdDim = TightFPDim<
		DIM1::WF + DIM2::WF + 1,
		DIM1::MAX_EXP + DIM2::MAX_EXP,
		DIM1::MIN_EXP + DIM1::MIN_EXP
	>;

template<vec_width TargetWF, typename DIM1, typename DIM2>
struct RoundedFPProd {
	private:
		using ExactProdDim = FPProdDim<DIM1, DIM2>;
		using round_helper = RoundDimHelper<ExactProdDim, TargetWF>;
	public:
		using dim = typename round_helper::dim;

	private:
		template<template<unsigned int, bool> class Wrapper>
		static inline FPNumber<ExactProdDim, Wrapper> do_exact_prod(FPNumber<DIM1, Wrapper> const & op1, FPNumber<DIM2, Wrapper> const & op2)
		{
			auto exp_no_overflow = op1.getExponent()
									  .addWithCarry(
										op2.getExponent(),
											{{0}}
										).template slice<ExactProdDim::WE - 1, 0>().as_signed();
			auto exp_overflow = op1.getExponent()
									  .addWithCarry(
										op2.getExponent(),
											{{1}}
										).template slice<ExactProdDim::WE - 1, 0>().as_signed();

			auto signif1 = op1.getSignificand();
			auto signif2 = op2.getSignificand();
			auto res = signif1 * signif2;
			auto overflow = res.template get<ExactProdDim::WF>();
			auto exp = Wrapper<ExactProdDim::WE, true>::mux(overflow, exp_overflow, exp_no_overflow);
			auto signif = Wrapper<ExactProdDim::WF, false>::mux(overflow,
				res.template slice<ExactProdDim::WF-1, 0>(),
				res.template slice<ExactProdDim::WF-2, 0>().concatenate(Wrapper<1, false>{0})
			);
			auto isNan1 = op1.isNaN();
			auto isNaN2 = op2.isNaN();
			auto isZero1 = op1.isZero();
			auto isZero2 = op2.isZero();
			auto isInf1 = op1.isInf();
			auto isInf2 = op2.isInf();
			auto oneIsNan = isNaN2 | isNan1;
			auto oneIsInf = isInf1 | isInf2;
			auto oneIsZero = isZero1 | isZero2;

			auto isNaN = oneIsNan | (oneIsInf & oneIsZero);
			auto isZero = oneIsZero & isNaN.invert();
			auto isInf = oneIsInf & isNaN.invert();
			auto res_sign = op1.getSign() ^ op2.getSign();

			return FPNumber<ExactProdDim, Wrapper>(signif, exp, res_sign, isInf, isNaN, isZero);
		}

		template<bool can_round, template<unsigned int, bool> class Wrapper>
		static inline FPNumber<dim, Wrapper> do_compute(FPNumber<DIM1, Wrapper> const & op1, FPNumber<DIM2, Wrapper> const & op2,
											  typename enable_if<not can_round>::type* = 0)
		{
			return do_exact_prod(op1, op2);
		}

		template<bool can_round, template<unsigned int, bool> class Wrapper>
		static inline FPNumber<dim, Wrapper> do_compute(FPNumber<DIM1, Wrapper> const & op1, FPNumber<DIM2, Wrapper> const & op2,
											  typename enable_if<can_round>::type* = 0)
		{
			using rounder = Rounder<dim, ExactProdDim>;
			auto exactProd = do_exact_prod(op1, op2);
			return rounder::compute(exactProd);
		}

		public:
		template<template<unsigned int, bool> class Wrapper>
		static inline FPNumber<dim, Wrapper> compute(FPNumber<DIM1, Wrapper> const & op1, FPNumber<DIM2, Wrapper> const & op2)
		{
			return do_compute<round_helper::CAN_ROUND>(op1, op2);
		}

};
#endif // PRODUCT_HPP

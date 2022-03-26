#include <boost/test/unit_test.hpp>

namespace utf = boost::unit_test;

#ifdef MPFR
#include <mpfr.h>
#endif

#include <hint.hpp>

using hint::VivadoWrapper;

#include "floatingpoint/expression.hpp"

#ifdef MPFR
BOOST_AUTO_TEST_CASE(TestRoundedProductFull, *utf::disabled() * utf::label("long"))
{
	constexpr vec_width WE1 = 4;
	constexpr vec_width WF1 = 7;

	constexpr vec_width WE2 = 5;
	constexpr vec_width WF2 = 6;

	constexpr vec_width targetWF = WF1+3;

	using dim1 = FixedFormat<WE1, WF1>;
	using dim2 = FixedFormat<WE2, WF2>;
	constexpr int64_t min_significand{0}, max_ex_significand1{1 << WF1}, max_ex_significand2{1 << WF2};

	/**** Using expr ********/

	using fpnum1 = FixedNumber<dim1, VivadoWrapper>;
	using fpnum2 = FixedNumber<dim2, VivadoWrapper>;


	/**** Using mpfr ****/
	mpfr_t mpfr_op1, mpfr_op2, mpfr_prod, mpfr_round;
	mpfr_init2(mpfr_op1, WF1 + 1);
	mpfr_init2(mpfr_op2, WF2 + 1);
	mpfr_init2(mpfr_prod, targetWF + 2);
	mpfr_init2(mpfr_round, targetWF + 2);
	for (int64_t exp1 = dim1::MIN_EXP ; exp1 <= dim1::MAX_EXP ; ++exp1) {
		for (int64_t frac1 = min_significand; (frac1 < max_ex_significand1); ++frac1) {
			uint64_t signif1 = frac1 | max_ex_significand1;
			mpfr_set_ui_2exp(mpfr_op1, signif1, exp1 - WF1, MPFR_RNDN);
			for (int64_t exp2 = dim2::MIN_EXP ; exp2 <= dim2::MAX_EXP ; ++exp2) {
				for (int64_t frac2 = min_significand; frac2 < max_ex_significand2 ; ++frac2) {
					uint64_t signif2 = frac2 | max_ex_significand2;
					mpfr_set_ui_2exp(mpfr_op2, signif2, exp2 - WF2, MPFR_RNDN);
					for (int64_t s1 = 0 ; s1 <= 1 ; s1++) {
						for(int64_t s2 = 0 ; s2 <=1 ; s2++) {
							/*** Expressions ***/
							fpnum2 op2{{{frac2}}, {{exp2}}, {{s2}}, {{0}}, {{0}}, {{0}}};
							fpnum1 op1{{{frac1}}, {{exp1}}, {{s1}}, {{0}}, {{0}}, {{0}}};

							auto eop1 = to_expr(op1);
							auto eop2 = to_expr(op2);

							auto exprres = eop1 * eop2;
							auto res = exprres.template computeWithTargetPrecision<targetWF>();
							bool res_is_zero = res.isZero().unravel();

							/*** MPFR ***/
							mpfr_setsign(mpfr_op1, mpfr_op1, s1, MPFR_RNDN);
							mpfr_setsign(mpfr_op2, mpfr_op2, s2, MPFR_RNDN);
							// MPFR does not provide rounding to nearest tie to away at the moment :
							// Work around : perform round down with one bit extra precision
							mpfr_mul(mpfr_prod, mpfr_op1, mpfr_op2, MPFR_RNDZ);
							auto mpfr_exp = mpfr_get_exp(mpfr_prod) - 1;
							mpfr_mul_2si(mpfr_round, mpfr_prod, targetWF - mpfr_exp, MPFR_RNDN);
							mpfr_abs(mpfr_round, mpfr_round, MPFR_RNDN);
							uint64_t mpfr_extra_prec_rd = mpfr_get_ui(mpfr_round, MPFR_RNDA);
							if (mpfr_extra_prec_rd > ((uint64_t{1} << (targetWF + 1))  - 1)) {
								++mpfr_exp;
							}
							bool mpfr_is_zero = (mpfr_extra_prec_rd == 0);
							uint64_t mpfr_res_frac = mpfr_extra_prec_rd & ((uint64_t{1} << targetWF) - 1);
							auto mpfr_res_sign = mpfr_signbit(mpfr_prod);

							bool ok_zero = not(mpfr_is_zero xor res_is_zero);
							bool fraction_ok = (res.getFraction().unravel() == mpfr_res_frac);
							bool exp_ok = (res.getExponent().unravel() == mpfr_exp) || res_is_zero;
							bool sign_ok = (res.getSign().unravel() == mpfr_res_sign) || res_is_zero;
							bool isInfOk = (res.isInf().unravel() == 0);
							bool isNaNOk = (res.isNaN().unravel() == 0);

							auto res_bool = fraction_ok and exp_ok and sign_ok and isInfOk and isNaNOk and ok_zero;

							if (not res_bool) {
								cerr << "Error with <" << std::abs(s1) << ", " << exp1 << ", " << frac1 << "> + <" << std::abs(s2) <<
										", " << exp2 << ", " << frac2 << ">" << endl;
								mpfr_dump(mpfr_op1);
								mpfr_dump(mpfr_op2);
								mpfr_dump(mpfr_prod);
								mpfr_dump(mpfr_round);
								BOOST_REQUIRE(res_bool);
							}
						}
					}
				}
			}
		}
	}
	mpfr_clears(mpfr_op1, mpfr_op2, mpfr_prod, mpfr_round, static_cast<mpfr_ptr>(0));
}
#else
BOOST_AUTO_TEST_CASE(WarningIncompleteTestsProduct)
{
	BOOST_CHECK_MESSAGE(false, "The current test executable was compiled without mpfr support, rounded product component will not be tested.");
}
#endif

BOOST_AUTO_TEST_CASE(TestProduct)
{
	constexpr vec_width WF1 = 15, WF2=12;
	constexpr int64_t maxExp1 = 6, minExp1 = -2, maxExp2=1, minExp2=-2;

	using dim1 = TightFixedFormat<WF1, maxExp1, minExp1>;
	using dim2 = TightFixedFormat<WF2, maxExp2, minExp2>;

	constexpr vec_width WE1 = dim1::WE, WE2 = dim2::WE;

	using proddim = FPProdDim<dim1, dim2>;

	constexpr vec_width WFProd = proddim::WF;

	constexpr uint64_t signif_1 = (uint64_t{1} << (WF1+1)) - 1;
	constexpr uint64_t signif_2 = (uint64_t{1} << (WF2+1)) - 1;

	constexpr uint64_t frac1 = signif_1 >> 1;
	constexpr uint64_t frac2 = signif_2 >> 1;

	constexpr uint64_t prod = signif_1 * signif_2;
	constexpr uint64_t prod_frac = prod & ((uint64_t{1} << (WF1 + WF2 + 1)) - 1);

	constexpr int64_t exp1 = 2, exp2 = 0;
	constexpr int64_t resexp = exp1+exp2+1;

	VivadoWrapper<WF1, false> wfrac1{{frac1}};
	VivadoWrapper<WF2, false> wfrac2{{frac2}};

	VivadoWrapper<WE1, true> wexp1{{exp1}};
	VivadoWrapper<WE2, true> wexp2{{exp2}};

	FixedNumber<dim1, VivadoWrapper> op1 {wfrac1, wexp1, {{0}}, {{0}}, {{0}}, {{0}}};
	FixedNumber<dim2, VivadoWrapper> op2 {wfrac2, wexp2, {{1}}, {{0}}, {{0}}, {{0}}};

	auto eop1 = to_expr(op1);
	auto eop2 = to_expr(op2);

	auto eprod = eop1 * eop2;
	auto res = eprod.computeWithTargetPrecision<WFProd>();

	BOOST_REQUIRE(res.getExponent().unravel() == resexp);
	BOOST_REQUIRE(res.getFraction().unravel() == prod_frac);
	BOOST_REQUIRE(res.getSign().unravel() == 1);
	BOOST_REQUIRE(res.isInf().unravel() == 0);
	BOOST_REQUIRE(res.isNaN().unravel() == 0);
	BOOST_REQUIRE(res.isZero().unravel() == 0);
}

BOOST_AUTO_TEST_CASE(TestProductFull, *utf::disabled() * utf::label("long"))
{
	constexpr vec_width WE1 = 4;
	constexpr vec_width WF1 = 7;

	constexpr vec_width WE2 = 5;
	constexpr vec_width WF2 = 6;

	using dim1 = FixedFormat<WE1, WF1>;
	using dim2 = FixedFormat<WE2, WF2>;

	constexpr int64_t min_significand{0}, max_ex_significand1{1 << WF1}, max_ex_significand2{1 << WF2};

	using proddim = FPProdDim<dim1, dim2>;

	constexpr vec_width WFProd = proddim::WF;

	constexpr uint64_t prod_mask{1 << WFProd};

	/**** Using expr ********/

	using fpnum1 = FixedNumber<dim1, VivadoWrapper>;
	using fpnum2 = FixedNumber<dim2, VivadoWrapper>;

	bool keep_going = true;

	for (int64_t exp1 = dim1::MIN_EXP ; exp1 <= dim1::MAX_EXP ; ++exp1) {
		for (int64_t frac1 = min_significand; (frac1 < max_ex_significand1) && keep_going; ++frac1) {
			uint64_t signif1 = frac1 | max_ex_significand1;
			for (int64_t exp2 = dim2::MIN_EXP ; exp2 <= dim2::MAX_EXP ; ++exp2) {
				for (int64_t frac2 = min_significand; frac2 < max_ex_significand2 ; ++frac2) {
					uint64_t signif2 = frac2 | max_ex_significand2;
					for (int64_t s1 : {0, 1}) {
						for(int64_t s2 : {0, 1}) {
							/*** Expressions ***/
							fpnum2 op2{{{frac2}}, {{exp2}}, {{s2}}, {{0}}, {{0}}, {{0}}};
							fpnum1 op1{{{frac1}}, {{exp1}}, {{s1}}, {{0}}, {{0}}, {{0}}};

							auto eop1 = to_expr(op1);
							auto eop2 = to_expr(op2);

							auto exprres = eop1 * eop2;
							auto res = exprres.template computeWithTargetPrecision<WFProd>();

							uint64_t prod = signif1*signif2;
							int64_t exp = exp1 + exp2;
							int64_t prod_overflow = (prod & prod_mask) >> WFProd;
							exp += prod_overflow;
							prod <<= (1-prod_overflow);
							uint64_t prod_frac = prod & (prod_mask - 1);

							bool exp_ok = (exp == res.getExponent().unravel());
							bool frac_ok = (res.getFraction().unravel() == prod_frac);
							bool sok = (res.getSign().unravel() == s1 ^ s2);
							bool infok = (res.isInf().unravel() == 0);
							bool nanOk = (res.isNaN().unravel() == 0);
							bool zeroOk = (res.isZero().unravel() == 0);
							bool resOK = exp_ok and frac_ok and sok and infok and nanOk and zeroOk;
							if (not resOK) {
								keep_going = false;
							}
						}
					}
				}
			}
		}
	}

	BOOST_REQUIRE(keep_going);
}

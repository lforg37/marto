#include <boost/test/unit_test.hpp>

namespace utf = boost::unit_test;

#ifdef MPFR
#include <mpfr.h>
#endif

#include <hint.hpp>

using hint::VivadoWrapper;

#include "floatingpoint/expression.hpp"

BOOST_AUTO_TEST_CASE(TestSumCompilation)
{
	constexpr vec_width CONST_WE = 2;
	constexpr vec_width EXTERNAL_WF = 23;
	constexpr vec_width INTERNAL_WF = 42;
	constexpr vec_width CONST_WF = 23;

	using vardim = TightFPDim<EXTERNAL_WF, 20, -20>;
	using const_dim = FPDim<CONST_WE, CONST_WF>;

	using fpnum = FPNumber<vardim, VivadoWrapper>;
	using constnum = FPNumber<const_dim, VivadoWrapper>;

	auto x = to_expr(fpnum::getZero());
	auto y = to_expr(fpnum::getZero());

	auto a = to_expr(constnum::getZero());
	auto b = to_expr(constnum::getZero());
	auto one = to_expr(constnum::getZero());

	auto xn_expr = y + one;// - a * x * x;
	auto yn_expr = b * x;

	auto xn = xn_expr.roundTo<vardim, INTERNAL_WF>();
	auto yn = yn_expr.roundTo<vardim, INTERNAL_WF>();
}

#ifdef MPFR
template<typename Dim1, typename Dim2, vec_width targetWF>
inline bool check_addition_standard(int64_t s1, int64_t s2, int64_t exp1, int64_t exp2, int64_t frac1, int64_t frac2)
{
	/**** Using expr ********/

	using fpnum1 = FPNumber<Dim1, VivadoWrapper>;
	using fpnum2 = FPNumber<Dim2, VivadoWrapper>;

	fpnum1 op1{{{frac1}}, {{exp1}}, {{s1}}, {{0}}, {{0}}, {{0}}};
	fpnum2 op2{{{frac2}}, {{exp2}}, {{s2}}, {{0}}, {{0}}, {{0}}};

	auto eop1 = to_expr(op1);
	auto eop2 = to_expr(op2);

	auto exprres = eop1+eop2;
	auto res = exprres.template computeWithTargetPrecision<targetWF>();
	bool res_is_zero = res.isZero().unravel();

	/**** Using mpfr ****/
	mpfr_t mpfr_op1, mpfr_op2, mpfr_sum, mpfr_round;
	mpfr_init2(mpfr_op1, Dim1::WF + 1);
	mpfr_init2(mpfr_op2, Dim2::WF + 1);
	mpfr_init2(mpfr_sum, targetWF + 2);
	mpfr_init2(mpfr_round, targetWF + 2);

	constexpr uint64_t signif_mask1 = {1 << Dim1::WF};
	uint64_t signif1 = frac1 | signif_mask1;
	mpfr_set_ui_2exp(mpfr_op1, signif1, exp1 - Dim1::WF, MPFR_RNDN);

	constexpr uint64_t signif_mask2 = {1 << Dim2::WF};
	uint64_t signif2 = frac2 | signif_mask2;
	mpfr_set_ui_2exp(mpfr_op2, signif2, exp2 - Dim2::WF, MPFR_RNDN);


	mpfr_setsign(mpfr_op1, mpfr_op1, s1, MPFR_RNDN);
	mpfr_setsign(mpfr_op2, mpfr_op2, s2, MPFR_RNDN);
	// MPFR does not provide rounding to nearest tie to away at the moment :
	// Work around : perform round down with one bit extra precision
	mpfr_add(mpfr_sum, mpfr_op1, mpfr_op2, MPFR_RNDZ);
	auto mpfr_exp = mpfr_get_exp(mpfr_sum) - 1;
	mpfr_mul_2si(mpfr_round, mpfr_sum, targetWF - mpfr_exp, MPFR_RNDN);
	mpfr_abs(mpfr_round, mpfr_round, MPFR_RNDN);
	uint64_t mpfr_extra_prec_rd = mpfr_get_ui(mpfr_round, MPFR_RNDA);
	if (mpfr_extra_prec_rd > ((uint64_t{1} << (targetWF + 1))  - 1)) {
		++mpfr_exp;
	}
	bool mpfr_is_zero = (mpfr_extra_prec_rd == 0);
	uint64_t mpfr_res_frac = mpfr_extra_prec_rd & ((uint64_t{1} << targetWF) - 1);
	auto mpfr_res_sign = mpfr_signbit(mpfr_sum);

	bool ok_zero = not(mpfr_is_zero xor res_is_zero);
	BOOST_CHECK_MESSAGE(ok_zero, "Error in non-zero result in mpfr or expression");

	bool fraction_ok = (res.getFraction().unravel() == mpfr_res_frac);
	BOOST_CHECK_MESSAGE(fraction_ok, "mpfr and expression results fraction differs");
	if (not fraction_ok) {
		cerr << "Expected : " << mpfr_res_frac << endl;
		cerr << "Got : " << res.getFraction().unravel() << endl;

	}

	bool exp_ok = (res.getExponent().unravel() == mpfr_exp) || res_is_zero;
	BOOST_CHECK_MESSAGE(exp_ok, "mpfr and expression result differs on the exponent");

	bool sign_ok = (res.getSign().unravel() == mpfr_res_sign) || res_is_zero;
	BOOST_CHECK_MESSAGE(sign_ok, "mpfr and expression results differs on the sign");

	bool isInfOk = (res.isInf().unravel() == 0);
	BOOST_CHECK_MESSAGE(isInfOk, "Result inf flag is set, that should not happen");

	bool isNaNOk = (res.isNaN().unravel() == 0);
	BOOST_CHECK_MESSAGE(isNaNOk, "Result NaN flag is set, that should not happen");

	auto res_bool = fraction_ok and exp_ok and sign_ok and isInfOk and isNaNOk and ok_zero;
	if (!res_bool) {
		mpfr_dump(mpfr_op1);
		mpfr_dump(mpfr_op2);
		mpfr_dump(mpfr_sum);
		mpfr_dump(mpfr_round);
	}
	mpfr_clears(mpfr_op1, mpfr_op2, mpfr_sum, mpfr_round, static_cast<mpfr_ptr>(0));
	return res_bool;
}

BOOST_AUTO_TEST_CASE(TestAdditionNoOverlap)
{
	using dim1 = TightFPDim<15, 18, -5>;
	using dim2 = TightFPDim<8, 7, -26>;
	check_addition_standard<dim1, dim2, 15>(0, 0, 15, -2, 0, 0);
}

BOOST_AUTO_TEST_CASE(TestAdditionSimpleOverlap)
{
	using dim1 = TightFPDim<15, 18, -5>;
	using dim2 = TightFPDim<8, 7, -26>;
	check_addition_standard<dim1, dim2, 15>(0, 0, 6, 4, 0, 0);
}

BOOST_AUTO_TEST_CASE(AdditionOpposite)
{
	constexpr vec_width WE = 5;
	constexpr vec_width WF = 9;
	constexpr vec_width targetWF = WF;
	using dim = FPDim<WE, WF>;

	bool res = check_addition_standard<dim, dim, targetWF>(
				0, 1, -16, -16, 0, 0
				);

	BOOST_REQUIRE(res);
}

BOOST_AUTO_TEST_CASE(AdditionDoubleNeg)
{
	constexpr vec_width WE = 5;
	constexpr vec_width WF = 9;
	constexpr vec_width targetWF = WF;
	using dim = FPDim<WE, WF>;

	bool res = check_addition_standard<dim, dim, targetWF>(
				1, 1, -16, -16, 0, 0
				);

	BOOST_REQUIRE(res);
}

BOOST_AUTO_TEST_CASE(AddOpSignOneBitCancel)
{
	constexpr vec_width WE = 5;
	constexpr vec_width WF = 9;
	constexpr vec_width targetWF = WF;
	using dim = FPDim<WE, WF>;

	bool res = check_addition_standard<dim, dim, targetWF>(
				0, 1, -16, -14, 0, 0
				);

	BOOST_REQUIRE(res);
}

BOOST_AUTO_TEST_CASE(TestAdditionCompleteFormat, *utf::disabled() * utf::label("long"))
{
	constexpr vec_width WE = 5;
	constexpr vec_width WF = 7;
	constexpr vec_width targetWF = WF;
	using dim = FPDim<WE, WF>;
	constexpr int64_t min_significand{0}, max_ex_significand{1 << WF};

	/**** Using expr ********/

	using fpnum = FPNumber<dim, VivadoWrapper>;


	/**** Using mpfr ****/
	constexpr uint64_t signif_mask{1 << WF};
	mpfr_t mpfr_op1, mpfr_op2, mpfr_sum, mpfr_round;
	mpfr_init2(mpfr_op1, WF + 1);
	mpfr_init2(mpfr_op2, WF + 1);
	mpfr_init2(mpfr_sum, targetWF + 2);
	mpfr_init2(mpfr_round, targetWF + 2);
	for (int64_t exp1 = dim::MIN_EXP ; exp1 <= dim::MAX_EXP ; ++exp1) {
		for (int64_t frac1 = min_significand; frac1 < max_ex_significand ; ++frac1) {
			uint64_t signif1 = frac1 | signif_mask;
			mpfr_set_ui_2exp(mpfr_op1, signif1, exp1 - WF, MPFR_RNDN);
			for (int64_t exp2 = exp1 ; exp2 <= dim::MAX_EXP ; ++exp2) {
				for (int64_t frac2 = frac1; frac2 < max_ex_significand ; ++frac2) {
					uint64_t signif2 = frac2 | signif_mask;
					mpfr_set_ui_2exp(mpfr_op2, signif2, exp2 - WF, MPFR_RNDN);
					for (int64_t s1 : {0, 1}) {
						for(int64_t s2 : {0, 1}) {
							/*** Expressions ***/
							fpnum op2{{{frac2}}, {{exp2}}, {{s2}}, {{0}}, {{0}}, {{0}}};
							fpnum op1{{{frac1}}, {{exp1}}, {{s1}}, {{0}}, {{0}}, {{0}}};

							auto eop1 = to_expr(op1);
							auto eop2 = to_expr(op2);

							auto exprres = eop1+eop2;
							auto res = exprres.template computeWithTargetPrecision<targetWF>();
							bool res_is_zero = res.isZero().unravel();

							/*** MPFR ***/
							mpfr_setsign(mpfr_op1, mpfr_op1, s1, MPFR_RNDN);
							mpfr_setsign(mpfr_op2, mpfr_op2, s2, MPFR_RNDN);
							// MPFR does not provide rounding to nearest tie to away at the moment :
							// Work around : perform round down with one bit extra precision
							mpfr_add(mpfr_sum, mpfr_op1, mpfr_op2, MPFR_RNDZ);
							auto mpfr_exp = mpfr_get_exp(mpfr_sum) - 1;
							mpfr_mul_2si(mpfr_round, mpfr_sum, targetWF - mpfr_exp, MPFR_RNDN);
							mpfr_abs(mpfr_round, mpfr_round, MPFR_RNDN);
							uint64_t mpfr_extra_prec_rd = mpfr_get_ui(mpfr_round, MPFR_RNDA);
							if (mpfr_extra_prec_rd > ((uint64_t{1} << (targetWF + 1))  - 1)) {
								++mpfr_exp;
							}
							bool mpfr_is_zero = (mpfr_extra_prec_rd == 0);
							uint64_t mpfr_res_frac = mpfr_extra_prec_rd & ((uint64_t{1} << targetWF) - 1);
							auto mpfr_res_sign = mpfr_signbit(mpfr_sum);

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
								mpfr_dump(mpfr_sum);
								mpfr_dump(mpfr_round);
								BOOST_REQUIRE(false);
							}
						}
					}
				}
			}
		}
	}
	mpfr_clears(mpfr_op1, mpfr_op2, mpfr_sum, mpfr_round, static_cast<mpfr_ptr>(0));
}

BOOST_AUTO_TEST_CASE(TestAdditionCompleteInequalFormats, *utf::disabled() * utf::label("long"))
{
	constexpr vec_width WE1 = 5;
	constexpr vec_width WF1 = 7;

	constexpr vec_width WE2 = 6;
	constexpr vec_width WF2 = 5;

	constexpr vec_width targetWF = WF1+2;

	using dim1 = FPDim<WE1, WF1>;
	using dim2 = FPDim<WE2, WF2>;
	constexpr int64_t min_significand{0}, max_ex_significand1{1 << WF1}, max_ex_significand2{1 << WF2};

	/**** Using expr ********/

	using fpnum1 = FPNumber<dim1, VivadoWrapper>;
	using fpnum2 = FPNumber<dim2, VivadoWrapper>;


	/**** Using mpfr ****/
	mpfr_t mpfr_op1, mpfr_op2, mpfr_sum, mpfr_round;
	mpfr_init2(mpfr_op1, WF1 + 1);
	mpfr_init2(mpfr_op2, WF2 + 1);
	mpfr_init2(mpfr_sum, targetWF + 2);
	mpfr_init2(mpfr_round, targetWF + 2);
	for (int64_t exp1 = dim1::MIN_EXP ; exp1 <= dim1::MAX_EXP ; ++exp1) {
		for (int64_t frac1 = min_significand; frac1 < max_ex_significand1 ; ++frac1) {
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

							auto exprres = eop1+eop2;
							auto res = exprres.template computeWithTargetPrecision<targetWF>();
							bool res_is_zero = res.isZero().unravel();

							/*** MPFR ***/
							mpfr_setsign(mpfr_op1, mpfr_op1, s1, MPFR_RNDN);
							mpfr_setsign(mpfr_op2, mpfr_op2, s2, MPFR_RNDN);
							// MPFR does not provide rounding to nearest tie to away at the moment :
							// Work around : perform round down with one bit extra precision
							mpfr_add(mpfr_sum, mpfr_op1, mpfr_op2, MPFR_RNDZ);
							auto mpfr_exp = mpfr_get_exp(mpfr_sum) - 1;
							mpfr_mul_2si(mpfr_round, mpfr_sum, targetWF - mpfr_exp, MPFR_RNDN);
							mpfr_abs(mpfr_round, mpfr_round, MPFR_RNDN);
							uint64_t mpfr_extra_prec_rd = mpfr_get_ui(mpfr_round, MPFR_RNDA);
							if (mpfr_extra_prec_rd > ((uint64_t{1} << (targetWF + 1))  - 1)) {
								++mpfr_exp;
							}
							bool mpfr_is_zero = (mpfr_extra_prec_rd == 0);
							uint64_t mpfr_res_frac = mpfr_extra_prec_rd & ((uint64_t{1} << targetWF) - 1);
							auto mpfr_res_sign = mpfr_signbit(mpfr_sum);

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
								mpfr_dump(mpfr_sum);
								mpfr_dump(mpfr_round);
								BOOST_REQUIRE(res_bool);
							}
						}
					}
				}
			}
		}
	}
	mpfr_clears(mpfr_op1, mpfr_op2, mpfr_sum, mpfr_round, static_cast<mpfr_ptr>(0));
}

#else
BOOST_AUTO_TEST_CASE(WarningIncompleteTestsSum)
{
	BOOST_CHECK_MESSAGE(false, "The current test executable was compiled without mpfr support, Sum component will not be tested.");
}
#endif

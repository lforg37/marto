#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestFPExpr

#include <boost/test/unit_test.hpp>

#ifdef MPFR
#include <mpfr.h>
#endif

#include <hint.hpp>

using hint::VivadoWrapper;

#include "floatingpoint/expression.hpp"
BOOST_AUTO_TEST_CASE(TestDimComputation) {
	using dim = TightFPDim<13, 5, -8>;
	static_assert(dim::WF == 13, "Error when computing tight dimension WF");
	static_assert(dim::WE == 4, "Error when computing tight fp dimension WE");
	static_assert(dim::MAX_STORABLE_EXP==7, "Error when computing tight FP Max_Storable_Exp");
	static_assert(dim::MIN_STORABLE_EXP==-8, "Error when computing tight FP dim Min_storable_exp");

	using dim2 = TightFPDim<13, 8, -8>;
	static_assert(dim2::WF == 13, "Error when computing tight dimension WF");
	static_assert(dim2::WE == 5, "Error when computing tight fp dimension WE");
	static_assert(dim2::MAX_STORABLE_EXP==15, "Error when computing tight FP Max_Storable_Exp");
	static_assert(dim2::MIN_STORABLE_EXP==-16, "Error when computing tight FP dim Min_storable_exp");

	using dim3 = FPDim<5, 13>;
	static_assert (dim3::EXP_CONSTRAINED == false, "Error with unconstrained FPDim");
	static_assert (dim3::MAX_EXP == 15, "Error with unconstrained FPDim");
}

BOOST_AUTO_TEST_CASE(TestOpposite) {
	constexpr vec_width WF = 13;
	constexpr vec_width WE = 5;
	using dim = FPDim<WE, WF>;
	using fpexpr = FPExpr<dim, VivadoWrapper>;
	using fpnum = FPNumber<dim, VivadoWrapper>;
	VivadoWrapper<WE, true> exp{{3}};
	VivadoWrapper<WF, false> signif {{(1 << WF) - 1}};
	VivadoWrapper<1, false> zero{{0}};

	fpnum initval{signif, exp, zero, zero, zero, zero};
	fpexpr a{initval};
	auto b = -a;
	auto rounded = b.computeWithTargetPrecision<WF>();
	BOOST_REQUIRE((rounded.getFraction() == initval.getFraction()).unravel());
	BOOST_REQUIRE((rounded.getExponent() == initval.getExponent()).unravel());
	BOOST_REQUIRE(rounded.getSign().unravel() == 1);
	BOOST_REQUIRE((rounded.isNaN() == initval.isNaN()).unravel());
	BOOST_REQUIRE((rounded.isZero() == initval.isZero()).unravel());
	BOOST_REQUIRE((rounded.isInf() == initval.isInf()).unravel());
}

BOOST_AUTO_TEST_CASE(TestProduct)
{
	constexpr vec_width WF1 = 15, WF2=12;
	constexpr int64_t maxExp1 = 6, minExp1 = -2, maxExp2=1, minExp2=-2;

	using dim1 = TightFPDim<WF1, maxExp1, minExp1>;
	using dim2 = TightFPDim<WF2, maxExp2, minExp2>;

	constexpr vec_width WE1 = dim1::WE, WE2 = dim2::WE;

	using proddim = FPProdDim<dim1, dim2>;

	constexpr vec_width WFProd = proddim::WF;
	constexpr vec_width WEProd = proddim::WE;

	static_assert (WFProd == 28, "Error in ProdDim WF computation");
	static_assert (WEProd == WE1, "Error in ProdDim WE computation");

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

	FPNumber<dim1, VivadoWrapper> op1 {wfrac1, wexp1, {{0}}, {{0}}, {{0}}, {{0}}};
	FPNumber<dim2, VivadoWrapper> op2 {wfrac2, wexp2, {{1}}, {{0}}, {{0}}, {{0}}};

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

BOOST_AUTO_TEST_CASE(TestAdditionCompleteFormat)
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
					for (int64_t s1 = 0 ; s1 <= 1 ; s1++) {
						for(int64_t s2 = 0 ; s2 <=1 ; s2++) {
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

							if (not res_bool) {
								cerr << "Error with <" << std::abs(s1) << ", " << exp1 << ", " << frac1 << "> + <" << std::abs(s2) <<
										", " << exp2 << ", " << frac2 << ">" << endl;
								mpfr_dump(mpfr_op1);
								mpfr_dump(mpfr_op2);
								mpfr_dump(mpfr_sum);
								mpfr_dump(mpfr_round);
							}
							BOOST_REQUIRE(res_bool);
						}
					}
				}
			}
		}
	}
	mpfr_clears(mpfr_op1, mpfr_op2, mpfr_sum, mpfr_round, static_cast<mpfr_ptr>(0));
}

BOOST_AUTO_TEST_CASE(TestAdditionCompleteInequalFormats)
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

							if (not res_bool) {
								cerr << "Error with <" << std::abs(s1) << ", " << exp1 << ", " << frac1 << "> + <" << std::abs(s2) <<
										", " << exp2 << ", " << frac2 << ">" << endl;
								mpfr_dump(mpfr_op1);
								mpfr_dump(mpfr_op2);
								mpfr_dump(mpfr_sum);
								mpfr_dump(mpfr_round);
							}
							BOOST_REQUIRE(res_bool);
						}
					}
				}
			}
		}
	}
	mpfr_clears(mpfr_op1, mpfr_op2, mpfr_sum, mpfr_round, static_cast<mpfr_ptr>(0));
}

#else
BOOST_AUTO_TEST_CASE(WarningIncompleteTests)
{
	std::cerr << "Warning : the test executable was not linked against mpfr, so test on SUM and Round operator will not be performed" << std::endl;
	BOOST_WARN_MESSAGE(true, "The current test executable was compiled without mpfr support, Sum and rounding component will not be tested.");
}
#endif

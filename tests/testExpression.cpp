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
void check_addition_standard(int64_t exp1, int64_t exp2, int64_t sfrac1, int64_t sfrac2)
{
	int sign1 = (sfrac1 < 0) ? 1 : 0;
	int sign2 = (sfrac2 < 0) ? 1 : 0;
	uint64_t frac1 = (sfrac1 < 0) ? -sfrac1 : sfrac1;
	uint64_t frac2 = (sfrac2 < 0) ? -sfrac2 : sfrac2;

	/**** Using expr ********/

	using fpnum1 = FPNumber<Dim1, VivadoWrapper>;
	using fpnum2 = FPNumber<Dim2, VivadoWrapper>;

	fpnum1 op1{{{frac1}}, {{exp1}}, {{sign1}}, {{0}}, {{0}}, {{0}}};
	fpnum2 op2{{{frac2}}, {{exp2}}, {{sign2}}, {{0}}, {{0}}, {{0}}};

	auto eop1 = to_expr(op1);
	auto eop2 = to_expr(op2);

	auto exprres = eop1+eop2;
	auto res = exprres.template computeWithTargetPrecision<targetWF>();

	/**** Using mpfr ****/
	mpfr_t mpfr_op1, mpfr_op2, mpfr_res;
	mpfr_init2(mpfr_op1, Dim1::WF);
	mpfr_init2(mpfr_op2, Dim2::WF);
	mpfr_init2(mpfr_res, targetWF);

	uint64_t signif1 = frac1 | (uint64_t{1} << Dim1::WF);
	uint64_t signif2 = frac2 | (uint64_t{1} << Dim2::WF);

	mpfr_set_ui_2exp(mpfr_op1, signif1, exp1 - Dim1::WF, MPFR_RNDN);
	mpfr_set_ui_2exp(mpfr_op2, signif2, exp2 - Dim2::WF, MPFR_RNDN);
	mpfr_setsign(mpfr_op1, mpfr_op1, sign1, MPFR_RNDN);
	mpfr_setsign(mpfr_op2, mpfr_op2, sign2, MPFR_RNDN);
	// MPFR does not provide rounding to nearest tie to away at the moment :
	mpfr_add(mpfr_res, mpfr_op1, mpfr_op2, MPFR_RNDN);
	auto mpfr_exp = mpfr_get_exp(mpfr_res) - 1;
	mpfr_mul_2si(mpfr_res, mpfr_res, targetWF-mpfr_exp, MPFR_RNDN);
	uint64_t mpfr_res_signif = mpfr_get_ui(mpfr_res, MPFR_RNDN);
	uint64_t mpfr_res_frac = mpfr_res_signif & ((uint64_t{1} << targetWF) - 1);
	auto mpfr_res_sign = mpfr_signbit(mpfr_res);
	mpfr_clears(mpfr_op1, mpfr_op2, mpfr_res, static_cast<mpfr_ptr>(0));
	//----------------------------------------------------------------------------

	BOOST_REQUIRE_MESSAGE(res.getFraction().unravel() == mpfr_res_frac, "mpfr and expression results fraction differs");
	BOOST_REQUIRE_MESSAGE(res.getExponent().unravel() == mpfr_exp, "mpfr and expression result differs on the exponent");
	BOOST_REQUIRE_MESSAGE(res.getSign().unravel() == mpfr_res_sign, "mpfr and expression results differs on the sign");
	BOOST_REQUIRE_MESSAGE(res.isInf().unravel() == 0, "Result inf flag is set, that should not happen");
	BOOST_REQUIRE_MESSAGE(res.isNaN().unravel() == 0, "Result NaN flag is set, that should not happen");
}

BOOST_AUTO_TEST_CASE(TestAdditionNoOverlap)
{
	using dim1 = TightFPDim<15, 18, -5>;
	using dim2 = TightFPDim<8, 7, -26>;
	check_addition_standard<dim1, dim2, 15>(15, -2, 0, 0);
}
#endif

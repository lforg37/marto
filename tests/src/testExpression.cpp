#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestFPExpr

#include <boost/test/unit_test.hpp>

namespace utf = boost::unit_test;

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




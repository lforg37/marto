#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestFPExpr

#include <boost/test/unit_test.hpp>

#include <hint.hpp>

using hint::VivadoWrapper;

#include "floatingpoint/expression.hpp"

BOOST_AUTO_TEST_CASE(Test)
{
	using dim = TightFPDim<13, 5, -8>;
	using leaf_content_t = FPNumber<dim, VivadoWrapper>;
	using leaf_t = FPExpr<dim, VivadoWrapper>;

	leaf_t a{leaf_content_t::getZero()};
	leaf_t b{leaf_content_t::getZero()};

	auto expr = a*b;
	auto rounded = expr.template computeWithTargetPrecision<15>();
}

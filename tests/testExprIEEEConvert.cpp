#include <cmath>
#include <cstdint>

#include <omp.h>

#include <boost/test/unit_test.hpp>

namespace utf = boost::unit_test;

#include "hint.hpp"

using hint::VivadoWrapper;

#include "ieeefloats/ieeetype.hpp"

union u32float {
	float f;
	uint32_t val;
};

inline bool test_round_val(uint32_t repr) {
	constexpr unsigned int WE = 8;
	constexpr unsigned int WF = 23;

	using ieeenum = IEEENumber<WE, WF, VivadoWrapper>;
	using ieeerounder = IEEEtoFPNum<WE, WF>;
	u32float fu32;
	bool nanOK, zeroOK, infOK, sOK, fracOK, expOK;
	int exp;
	int expIsZero, expIsInf, expIsNaN;
	expIsInf = expIsNaN = expIsZero = 0;
	ieeenum val_marto = ieeenum{{repr}};
	auto rounded = ieeerounder::compute(val_marto);

	constexpr float scale = static_cast<float>(1 << (WF+1));

	fu32.val = repr;
	auto normalised = frexp(fu32.f, &exp);
	sOK = (signbit(fu32.val) == rounded.getSign().unravel());

	if (isnan(normalised)) {
		expIsNaN = 1;
		fracOK = true;
		expOK = true;
	} else if (isinf(normalised)) {
		expIsInf = 1;
		fracOK = true;
		expOK = true;
	} else if (iszero(normalised)) {
		expIsZero = 1;
		fracOK = true;
		expOK = true;
	} else {
		normalised = fabs(normalised);
		auto significand = static_cast<uint64_t>(normalised * scale);
		auto nexp = exp - 1;
		expOK = (nexp == rounded.getExponent().unravel());
		fracOK = (significand == rounded.getSignificand().unravel());
	}

	nanOK = (rounded.isNaN().unravel() == expIsNaN);
	zeroOK = (rounded.isZero().unravel() == expIsZero);
	infOK = (rounded.isInf().unravel() == expIsInf);
	bool ok = (nanOK and zeroOK and infOK and expOK and fracOK);
	if (not ok) {
		cerr << repr << endl;
		cerr << "===" << endl;
		cerr << nanOK << "\n" << zeroOK << "\n" << infOK << "\n" << expOK << "\n" << fracOK << endl;
	}
	return ok;
}

BOOST_AUTO_TEST_CASE(TestIEEEBinary32ToFPExpr, *utf::disabled() * utf::label("long"))
{
	bool total_ok = true;

	constexpr uint64_t limit = uint64_t{1} << 32;

	#pragma omp parallel
	{
		uint64_t inc;
		uint64_t start;
		#pragma omp critical
		{
			inc = omp_get_num_threads();
			start = omp_get_thread_num();
		}


		for (uint64_t repr_ext = start; repr_ext < limit && total_ok; repr_ext += inc) {
			uint32_t repr = static_cast<uint32_t>(repr_ext);

			auto ok = test_round_val(repr);
			if (!ok) {
				total_ok = false;
			}
		}
	}
	BOOST_REQUIRE(total_ok);
}

BOOST_AUTO_TEST_CASE(testRoundFP32SpecificVal) {
	BOOST_REQUIRE(test_round_val(2147483652));
	// TODO populate with interesting values
}

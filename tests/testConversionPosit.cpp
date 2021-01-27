#ifndef VIVADO_BACKEND
#define VIVADO_BACKEND
#endif

#include <boost/test/unit_test.hpp>

#include "hint.hpp"
#include "tools/printing.hpp"

#include "posit/posit_decoder.hpp"
#include "posit/posit_encoder.hpp"
#include "posit/posit_in_place_round.hpp"

using namespace std;
using hint::to_string;
namespace utf = boost::unit_test;
#define Wrapper hint::VivadoWrapper

inline void check_inplace_val(uint32_t value)
{
	PositIntermediateFormat<16, 1, Wrapper, false> unrounded_pif {{value}};
	auto inplace_rounded = in_place_rounder(unrounded_pif);
	auto direct_encoded = posit_encoder(unrounded_pif);
	auto decoded = posit_decoder(direct_encoded);

	if (unrounded_pif.getIsNaR().unravel()) {
		BOOST_REQUIRE_MESSAGE(inplace_rounded.getIsNaR().unravel(),
							  "Error with value " << value <<
							  "Expecting NaR, isNaR bit unset"
							  );
	} else {
		BOOST_REQUIRE_MESSAGE(decoded.unravel() == inplace_rounded.unravel(),
							  "Error with value " << value <<
							  "\ngot:\n" <<
							  inplace_rounded.unravel() <<
							  "\nInstead of:\n" <<
							  decoded.unravel()
							  );
	}
}

BOOST_AUTO_TEST_CASE(testSomeInplaceRoundP16)
{
	check_inplace_val(446464);
	check_inplace_val(24576);
	check_inplace_val(20481);
	check_inplace_val(45055);
	check_inplace_val(57343);
	check_inplace_val(434176);
	check_inplace_val(1069056);
}

BOOST_AUTO_TEST_CASE(testInplaceRoundP16, *utf::label("long"))
{
	constexpr size_t N = 16;
	constexpr size_t WE = 6;
	constexpr size_t WF = 12;
	constexpr size_t limit_frac = (1 << (WF));


	constexpr size_t PIF_SIZE = StandardPositDim<16>::UPIFSize;
	for (size_t isNar = 0 ; isNar <= 1 ; isNar++) {
		size_t isNarMask = isNar << (WF + 2 + WE);
		for (size_t gs = 0 ; gs < 4 ; gs++) {
			size_t gsmask = gs << (WF + WE + 3) | isNarMask;
			for (size_t exp = 1 ; exp < (N-2) << 1 ; exp++) {
				size_t exp_mask  = exp << (WF+2) | gsmask;
				for (size_t si = 1; si < 3 ; ++si) {
					size_t exp_si_mask = (si << WF) | exp_mask;
					for (size_t frac = 0 ; frac <  limit_frac ; ++frac) {
						size_t val = frac | exp_si_mask;
						check_inplace_val(val);
					}
				}
			}
		}
	}
}

#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE PositArithmeticOps

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/test/unit_test.hpp>
#include "softposit.h"

#include "posit_decoder.hpp"
#include "posit_encoder.hpp"
#include "posit_dim.hpp"
#include "posit_add.hpp"
#include "lzoc_shifter.hpp"

using namespace std;

BOOST_AUTO_TEST_CASE(TestAllSumPosit8) 
{
	uint16_t value = 0;
	auto positZero = posit_decoder(PositEncoding<16> (0));

	do {
		auto valueEncoding = PositEncoding<16> (value);
		auto decoded = posit_decoder(valueEncoding);
		auto sum = posit_add(decoded, positZero);
		auto encoded = posit_encoder(sum);
		if(!(encoded == valueEncoding)){
			fprintf(stderr, "\n\n\n\n");
			fprintf(stderr, "=== Expected result === \n");
			decoded.printContent();
			fprintf(stderr, "=== Computed result === \n");
			sum.printContent();
		}
		BOOST_REQUIRE_MESSAGE(encoded == valueEncoding, "Sum of " << value << " and 0 returned " << (unsigned int)encoded << " while it should have returned " << (unsigned int)valueEncoding);
		value++;
	} while (value != 0);
}

#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE PositArithmeticOps

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/test/unit_test.hpp>
#include "softposit.h"

#include "posit_decoder.hpp"
#include "posit_dim.hpp"
#include "lzoc_shifter.hpp"

using namespace std;

BOOST_AUTO_TEST_CASE(TestAllSumPosit8) 
{
	uint8_t left = 0;
	do {

		left++;
	} while (left != 0);
	BOOST_REQUIRE_MESSAGE(false, "Please, write the test");
}


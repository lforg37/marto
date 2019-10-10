#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestIEEEMultiplier

#include <fstream>
#include <string.h>
#include <iostream>
#include <stdio.h>

#include <boost/test/unit_test.hpp>

#include "hint.hpp"
#include "tools/printing.hpp"
#include "ieeefloats/ieee_multiplier.hpp"
#include <omp.h>

#include "ieeefloats/ieeetype.hpp"

#include "numeric_formats/ieee_small.hpp"

using namespace std;
using namespace hint;
#include <bitset>


BOOST_AUTO_TEST_CASE(TestIEEEMul0)
{
	constexpr unsigned int WE = 3;
	constexpr unsigned int WF = 4;

	using SIEEE = libnumform::SmallIEEENumber<WE, WF>;

	SIEEE a{0b00000011}, b{0b00001011};

	auto c = a*b;

	BOOST_REQUIRE_MESSAGE(c.getRepr()==0, "Error");
}


BOOST_AUTO_TEST_CASE(TestIEEMul_3_4)
{
	constexpr unsigned int WE = 3;
	constexpr unsigned int WF = 4;

	using SIEEE = libnumform::SmallIEEENumber<WE, WF>;


	// half_float::half var1, var2, half_prod;
	// uint16_t var2_ui16 = 0b0000000000000011;
	// uint16_t var1_ui16 = 0b0110000101010110;

	// var1.data_ = var1_ui16;
	// var2.data_ = var2_ui16;
	// std::bitset<16> x1(var1_ui16);
	// std::bitset<16> x2(var2_ui16);
	// cerr << "var1 " << x1 << endl;
	// cerr << "var2 " << x2 << endl;
	// half_prod = var1 * var2;
	// ap_uint<16> resgmp{half_prod.data_};
	// IEEENumber<WE, WF, hint::VivadoWrapper> ieee_res_gmp{resgmp};

	// ap_uint<16> var1gmp{var1_ui16};
	// ap_uint<16> var2gmp{var2_ui16};
	// IEEENumber<WE, WF, hint::VivadoWrapper> var1ieee{var1gmp};
	// IEEENumber<WE, WF, hint::VivadoWrapper> var2ieee{var2gmp};
	// auto marto_prod = ieee_product(var1ieee, var2ieee);

	// cerr << "half \t" << hint::to_string(static_cast<hint::VivadoWrapper<16, false> >(ieee_res_gmp) ) << endl;
	// cerr << "marto \t" << hint::to_string(static_cast<hint::VivadoWrapper<16, false> >(marto_prod) ) << endl;


	// omp_set_num_threads(16);
	constexpr unsigned int FORMAT_SIZE = 1 + WE + WF;
	constexpr uint32_t FORMAT_LIMIT = 1 << FORMAT_SIZE;
	SIEEE op1{0}, op2{0};

	//#pragma omp parallel for
	for (uint32_t count1 = 0 ; count1 < FORMAT_LIMIT ; ++count1) {
		uint32_t op1_repr = count1;
		mpz_class var1gmp{op1_repr};
		IEEENumber<WE, WF, hint::GMPWrapper> var1ieee{var1gmp};
		op1 = op1_repr;

		for (uint32_t count2=0 ; count2 < FORMAT_LIMIT ; count2++ ) {
			std::cerr << count1 << ", " << count2 << endl;
			uint32_t op2_repr = count2;
			mpz_class var2gmp{op2_repr};
			IEEENumber<WE, WF, hint::GMPWrapper> var2ieee{var2gmp};
			op2 = op2_repr;
			auto prod = op1 * op2;
			auto prod_repr = prod.getRepr();

			auto marto_prod = ieee_product(var1ieee, var2ieee);
			if (prod.isNaN()){//result is NaN
				bool marto_is_nan = (marto_prod.isNaN().unravel() == 1);
				//#pragma omp critical
				BOOST_REQUIRE_MESSAGE(marto_is_nan, "NAN CASE Error for \t" <<
										hint::to_string(static_cast<GMPWrapper<FORMAT_SIZE, false> >(var1ieee) ) <<
										" and " << hint::to_string(static_cast<GMPWrapper<FORMAT_SIZE, false> >(var2ieee) ));
			}
			else{
				bool must = (prod.getRepr()==marto_prod.unravel());
				// bool must = ((marto_prod.unravel()-resgmp) <= 1) and ((resgmp-marto_prod.unravel()) <= 1);
				if(not must){
					//#pragma omp critical
					BOOST_REQUIRE_MESSAGE(false, "Error for" <<
										  hint::to_string(static_cast<GMPWrapper<FORMAT_SIZE, false> >(var1ieee) ) <<
										  " and " << hint::to_string(static_cast<GMPWrapper<FORMAT_SIZE, false> >(var2ieee) ) <<
										  "\nexpecting\t" << hint::to_string(GMPWrapper<FORMAT_SIZE, false>{prod_repr} ) <<
										  "\ngot  \t" << hint::to_string(static_cast<GMPWrapper<FORMAT_SIZE, false> >(marto_prod) ));
				}
			}
		}
	}
}//*/

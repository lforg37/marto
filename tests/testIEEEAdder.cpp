#include <fstream>
#include <string.h>
#include <iostream>
#include <stdio.h>

#include <boost/test/unit_test.hpp>

#include "hint.hpp"
#include "tools/printing.hpp"
#include "ieeefloats/ieee_adder.hpp"

#include "ieeefloats/ieeetype.hpp"

#include "numeric_formats/ieee_small.hpp"

using namespace std;
using namespace hint;
namespace utf = boost::unit_test;
constexpr int PRINT_EVERY = 100;
#include <bitset>

#ifdef SOFTFLOAT
extern "C" {
#include "softfloat.h"
}

inline bool isNan(float16_t val)
{
	constexpr uint16_t nan_inf_flag = ((1 << 5) - 1) << 10;
	constexpr uint16_t absval_flag = (1 << 15) - 1;
	uint16_t repr = val.v & absval_flag;
	return ((repr & nan_inf_flag) == nan_inf_flag) && (repr ^ nan_inf_flag);
}

BOOST_AUTO_TEST_CASE(TestIEEEAdd_5_10_SP_RNTE, *utf::disabled() * utf::label("long"))
{
	constexpr unsigned int WE = 5;
	constexpr unsigned int WF = 10;

	using MartoIEEE = IEEENumber<WE, WF, hint::VivadoWrapper>;
	constexpr uint64_t FORMAT_SIZE = 1 + WE + WF;
	constexpr uint64_t FORMAT_LIMIT = 1 << FORMAT_SIZE;
	uint64_t counter = 0;
	softfloat_roundingMode = softfloat_round_near_even;
	auto marto_roundingmode = IEEERoundingMode::RoundNearestTieEven;
	for (uint64_t count1 = 0 ; count1 < FORMAT_LIMIT ; ++count1) {
		uint16_t op1_repr = count1;
		float16_t op1_sf{op1_repr};
		MartoIEEE op1_marto{{op1_repr}};

		for (uint32_t count2=0 ; count2 < FORMAT_LIMIT ; count2++ ) {
			uint16_t op2_repr = count2;
			float16_t op2_sf{op2_repr};
			MartoIEEE op2_marto{{op2_repr}};
			auto sum_sf = f16_add(op1_sf, op2_sf);
			auto sum_repr = sum_sf.v;
			auto sum_marto = ieee_add_sub_impl(op1_marto, op2_marto, marto_roundingmode);
			if (isNan(sum_sf)) {//result is NaN
				bool marto_is_nan = (sum_marto.isNaN().unravel() == 1);
				if(not marto_is_nan){
					#pragma omp critical
					BOOST_REQUIRE_MESSAGE(false, "NAN CASE Error for \t" <<
											hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op1_marto) ) <<
											" and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op2_marto) ));
				}
			}
			else{
				bool must = (sum_sf.v == sum_marto.unravel());
				// bool must = ((marto_prod.unravel()-resgmp) <= 1) and ((resgmp-marto_prod.unravel()) <= 1);
				if(not must){
					BOOST_REQUIRE_MESSAGE(false, "Error for " <<
										  hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op1_marto) ) <<
										  " and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op2_marto) ) <<
										  "\nexpecting\t" << hint::to_string(VivadoWrapper<FORMAT_SIZE, false>{sum_repr} ) <<
										  "\ngot  \t" << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(sum_marto) ));
				}
			}
		}
		counter++;
		if((counter % PRINT_EVERY) == 0){
			counter = 0;
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(FORMAT_LIMIT)*100, counter,FORMAT_LIMIT);
		}
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)\n", static_cast<double>(FORMAT_LIMIT)/static_cast<double>(FORMAT_LIMIT)*100, FORMAT_LIMIT,FORMAT_LIMIT);
}

BOOST_AUTO_TEST_CASE(TestIEEEAdd_5_10_SP_RNUp, *utf::disabled() * utf::label("long"))
{
	constexpr unsigned int WE = 5;
	constexpr unsigned int WF = 10;

	using MartoIEEE = IEEENumber<WE, WF, hint::VivadoWrapper>;
	constexpr uint64_t FORMAT_SIZE = 1 + WE + WF;
	constexpr uint64_t FORMAT_LIMIT = 1 << FORMAT_SIZE;
	uint64_t counter = 0;
	softfloat_roundingMode = softfloat_round_max;
	auto marto_roundingmode = IEEERoundingMode::RoundUp;
	for (uint64_t count1 = 0 ; count1 < FORMAT_LIMIT ; ++count1) {
		uint16_t op1_repr = count1;
		float16_t op1_sf{op1_repr};
		MartoIEEE op1_marto{{op1_repr}};

		for (uint32_t count2=0 ; count2 < FORMAT_LIMIT ; count2++ ) {
			uint16_t op2_repr = count2;
			float16_t op2_sf{op2_repr};
			MartoIEEE op2_marto{{op2_repr}};
			auto sum_sf = f16_add(op1_sf, op2_sf);
			auto sum_repr = sum_sf.v;
			auto sum_marto = ieee_add_sub_impl(op1_marto, op2_marto, marto_roundingmode);
			if (isNan(sum_sf)) {//result is NaN
				bool marto_is_nan = (sum_marto.isNaN().unravel() == 1);
				if(not marto_is_nan){
					BOOST_REQUIRE_MESSAGE(false, "NAN CASE Error for \t" <<
											hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op1_marto) ) <<
											" and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op2_marto) ));
				}
			}
			else{
				bool must = (sum_sf.v == sum_marto.unravel());
				// bool must = ((marto_prod.unravel()-resgmp) <= 1) and ((resgmp-marto_prod.unravel()) <= 1);
				if(not must){
					BOOST_REQUIRE_MESSAGE(false, "Error for " <<
										  hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op1_marto) ) <<
										  " and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op2_marto) ) <<
										  "\nexpecting\t" << hint::to_string(VivadoWrapper<FORMAT_SIZE, false>{sum_repr} ) <<
										  "\ngot  \t" << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(sum_marto) ));
				}
			}
		}
		counter++;
		if((counter % PRINT_EVERY) == 0){
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(FORMAT_LIMIT)*100, counter,FORMAT_LIMIT);
		}
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)\n", static_cast<double>(FORMAT_LIMIT)/static_cast<double>(FORMAT_LIMIT)*100, FORMAT_LIMIT,FORMAT_LIMIT);
}

BOOST_AUTO_TEST_CASE(TestIEEEAdd_5_10_SP_RNDown, *utf::disabled() * utf::label("long"))
{
	constexpr unsigned int WE = 5;
	constexpr unsigned int WF = 10;

	using MartoIEEE = IEEENumber<WE, WF, hint::VivadoWrapper>;
	constexpr uint64_t FORMAT_SIZE = 1 + WE + WF;
	constexpr uint64_t FORMAT_LIMIT = 1 << FORMAT_SIZE;
	uint64_t counter = 0;
	softfloat_roundingMode = softfloat_round_min;
	auto marto_roundingmode = IEEERoundingMode::RoundDown;
	for (uint64_t count1 = 0 ; count1 < FORMAT_LIMIT ; ++count1) {
		uint16_t op1_repr = count1;
		float16_t op1_sf{op1_repr};
		MartoIEEE op1_marto{{op1_repr}};

		for (uint32_t count2=0 ; count2 < FORMAT_LIMIT ; count2++ ) {
			uint16_t op2_repr = count2;
			float16_t op2_sf{op2_repr};
			MartoIEEE op2_marto{{op2_repr}};
			auto sum_sf = f16_add(op1_sf, op2_sf);
			auto sum_repr = sum_sf.v;
			auto sum_marto = ieee_add_sub_impl(op1_marto, op2_marto, marto_roundingmode);
			if (isNan(sum_sf)) {//result is NaN
				bool marto_is_nan = (sum_marto.isNaN().unravel() == 1);
				if(not marto_is_nan){
					BOOST_REQUIRE_MESSAGE(false, "NAN CASE Error for \t" <<
											hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op1_marto) ) <<
											" and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op2_marto) ));
				}
			}
			else{
				bool must = (sum_sf.v == sum_marto.unravel());
				// bool must = ((marto_prod.unravel()-resgmp) <= 1) and ((resgmp-marto_prod.unravel()) <= 1);
				if(not must){
					BOOST_REQUIRE_MESSAGE(false, "Error for " <<
										  hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op1_marto) ) <<
										  " and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op2_marto) ) <<
										  "\nexpecting\t" << hint::to_string(VivadoWrapper<FORMAT_SIZE, false>{sum_repr} ) <<
										  "\ngot  \t" << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(sum_marto) ));
				}
			}
		}
		counter++;
		if((counter % PRINT_EVERY) == 0){

			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(FORMAT_LIMIT)*100, counter,FORMAT_LIMIT);
		}
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)\n", static_cast<double>(FORMAT_LIMIT)/static_cast<double>(FORMAT_LIMIT)*100, FORMAT_LIMIT,FORMAT_LIMIT);
}

BOOST_AUTO_TEST_CASE(TestIEEEAdd_5_10_SP_RNDZero, *utf::disabled() * utf::label("long"))
{
	constexpr unsigned int WE = 5;
	constexpr unsigned int WF = 10;

	using MartoIEEE = IEEENumber<WE, WF, hint::VivadoWrapper>;
	constexpr uint64_t FORMAT_SIZE = 1 + WE + WF;
	constexpr uint64_t FORMAT_LIMIT = 1 << FORMAT_SIZE;
	uint64_t counter = 0;
	softfloat_roundingMode = softfloat_round_minMag;
	auto marto_roundingmode = IEEERoundingMode::RoundTowardZero;
	for (uint64_t count1 = 0 ; count1 < FORMAT_LIMIT ; ++count1) {
		uint16_t op1_repr = count1;
		float16_t op1_sf{op1_repr};
		MartoIEEE op1_marto{{op1_repr}};

		for (uint32_t count2=0 ; count2 < FORMAT_LIMIT ; count2++ ) {
			uint16_t op2_repr = count2;
			float16_t op2_sf{op2_repr};
			MartoIEEE op2_marto{{op2_repr}};
			auto sum_sf = f16_add(op1_sf, op2_sf);
			auto sum_repr = sum_sf.v;
			auto sum_marto = ieee_add_sub_impl(op1_marto, op2_marto, marto_roundingmode);
			if (isNan(sum_sf)) {//result is NaN
				bool marto_is_nan = (sum_marto.isNaN().unravel() == 1);
				if(not marto_is_nan){
					BOOST_REQUIRE_MESSAGE(false, "NAN CASE Error for \t" <<
											hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op1_marto) ) <<
											" and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op2_marto) ));
				}
			}
			else{
				bool must = (sum_sf.v == sum_marto.unravel());
				// bool must = ((marto_prod.unravel()-resgmp) <= 1) and ((resgmp-marto_prod.unravel()) <= 1);
				if(not must){
					BOOST_REQUIRE_MESSAGE(false, "Error for " <<
										  hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op1_marto) ) <<
										  " and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op2_marto) ) <<
										  "\nexpecting\t" << hint::to_string(VivadoWrapper<FORMAT_SIZE, false>{sum_repr} ) <<
										  "\ngot  \t" << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(sum_marto) ));
				}
			}
		}
		counter++;
		if((counter % PRINT_EVERY) == 0){
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(FORMAT_LIMIT)*100, counter,FORMAT_LIMIT);
		}
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)\n", static_cast<double>(FORMAT_LIMIT)/static_cast<double>(FORMAT_LIMIT)*100, FORMAT_LIMIT,FORMAT_LIMIT);
}

BOOST_AUTO_TEST_CASE(TestIEEEAdd_5_10_SP_RNTA, *utf::disabled() * utf::label("long"))
{
	constexpr unsigned int WE = 5;
	constexpr unsigned int WF = 10;

	using MartoIEEE = IEEENumber<WE, WF, hint::VivadoWrapper>;
	constexpr uint64_t FORMAT_SIZE = 1 + WE + WF;
	constexpr uint64_t FORMAT_LIMIT = 1 << FORMAT_SIZE;
	uint64_t counter = 0;
	softfloat_roundingMode = softfloat_round_near_maxMag;
	auto marto_roundingmode = IEEERoundingMode::RoundNearestTieAway;
	for (uint64_t count1 = 0 ; count1 < FORMAT_LIMIT ; ++count1) {
		uint16_t op1_repr = count1;
		float16_t op1_sf{op1_repr};
		MartoIEEE op1_marto{{op1_repr}};

		for (uint32_t count2=0 ; count2 < FORMAT_LIMIT ; count2++ ) {
			uint16_t op2_repr = count2;
			float16_t op2_sf{op2_repr};
			MartoIEEE op2_marto{{op2_repr}};
			auto sum_sf = f16_add(op1_sf, op2_sf);
			auto sum_repr = sum_sf.v;
			auto sum_marto = ieee_add_sub_impl(op1_marto, op2_marto, marto_roundingmode);
			if (isNan(sum_sf)) {//result is NaN
				bool marto_is_nan = (sum_marto.isNaN().unravel() == 1);
				if(not marto_is_nan){
					BOOST_REQUIRE_MESSAGE(false, "NAN CASE Error for \t" <<
											hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op1_marto) ) <<
											" and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op2_marto) ));
				}
			}
			else{
				bool must = (sum_sf.v == sum_marto.unravel());
				// bool must = ((marto_prod.unravel()-resgmp) <= 1) and ((resgmp-marto_prod.unravel()) <= 1);
				if(not must){
					BOOST_REQUIRE_MESSAGE(false, "Error for " <<
										  hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op1_marto) ) <<
										  " and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(op2_marto) ) <<
										  "\nexpecting\t" << hint::to_string(VivadoWrapper<FORMAT_SIZE, false>{sum_repr} ) <<
										  "\ngot  \t" << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(sum_marto) ));
				}
			}
		}
		counter++;
		if((counter % PRINT_EVERY) == 0){
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(FORMAT_LIMIT)*100, counter,FORMAT_LIMIT);
		}
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)\n", static_cast<double>(FORMAT_LIMIT)/static_cast<double>(FORMAT_LIMIT)*100, FORMAT_LIMIT,FORMAT_LIMIT);
}
#endif

BOOST_AUTO_TEST_CASE(TestIEEEAdd_4_7)
{
	constexpr unsigned int WE = 4;
	constexpr unsigned int WF = 7;

	using SIEEE = libnumform::SmallIEEENumber<WE, WF>;

	constexpr unsigned int FORMAT_SIZE = 1 + WE + WF;
	constexpr uint32_t FORMAT_LIMIT = 1 << FORMAT_SIZE;
	uint32_t counter = 0;
	for (uint32_t count1 = 0 ; count1 < FORMAT_LIMIT ; ++count1) {
		SIEEE op1{0};
		uint32_t op1_repr = count1;
		IEEENumber<WE, WF, hint::VivadoWrapper> var1ieee{{op1_repr}};
		op1 = op1_repr;

		for (uint32_t count2=0 ; count2 < FORMAT_LIMIT ; count2++ ) {
			SIEEE op2{0};
			// std::cerr << count1 << ", " << count2 << endl;
			uint32_t op2_repr = count2;
			IEEENumber<WE, WF, hint::VivadoWrapper> var2ieee{{op2_repr}};
			op2 = op2_repr;
			auto sum = op1 + op2;
			auto sum_repr = sum.getRepr();

			auto marto_sum = ieee_add_sub_impl(var1ieee, var2ieee);
			if (sum.isNaN()){//result is NaN
				bool marto_is_nan = (marto_sum.isNaN().unravel() == 1);
				if(not marto_is_nan){
					BOOST_REQUIRE_MESSAGE(false, "NAN CASE Error for \t" <<
											hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(var1ieee) ) <<
											" and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(var2ieee) ));
				}
			}
			else{
				bool must = (sum.getRepr()==marto_sum.unravel());
				// bool must = ((marto_prod.unravel()-resgmp) <= 1) and ((resgmp-marto_prod.unravel()) <= 1);
				if(not must){
					BOOST_REQUIRE_MESSAGE(false, "Error for " <<
										  hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(var1ieee) ) <<
										  " and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(var2ieee) ) <<
										  "\nexpecting\t" << hint::to_string(VivadoWrapper<FORMAT_SIZE, false>{sum_repr} ) <<
										  "\ngot  \t" << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(marto_sum) ));
				}
			}
		}
		counter++;
		if((counter % PRINT_EVERY) == 0){
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)", static_cast<double>(counter)/static_cast<double>(FORMAT_LIMIT)*100, counter,FORMAT_LIMIT);
		}
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)\n", static_cast<double>(FORMAT_LIMIT)/static_cast<double>(FORMAT_LIMIT)*100, FORMAT_LIMIT,FORMAT_LIMIT);

}//*/

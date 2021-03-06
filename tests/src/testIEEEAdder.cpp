#include <fstream>
#include <string.h>
#include <iostream>
#include <stdio.h>

#include <boost/test/unit_test.hpp>

#include "hint.hpp"
#include "tools/printing.hpp"
#include "ieeefloats/ieee_adder.hpp"

#ifdef HAS_OMP
#include <omp.h>
#endif 

#include "ieeefloats/ieeetype.hpp"

#include "numeric_formats/ieee_small.hpp"

using namespace std;
using namespace hint;
namespace utf = boost::unit_test;
constexpr int PRINT_EVERY = 10;
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

struct SumError {
		enum struct Code : uint8_t {
			OK = 0,
			WaitingNaN = 1,
			ResDiffer = 2
		};

		uint16_t op1;
		uint16_t op2;
		uint16_t expected_result;
		uint16_t result;
		Code err_code;
};

bool test_f16_sum(uint16_t op1_repr, uint16_t op2_repr, decltype(softfloat_roundingMode) sf_rnd_mode, IEEERoundingMode marto_round_mode)
{
	float16_t op1_sf{op1_repr}, op2_sf{op2_repr};
	IEEENumber<5, 10, hint::VivadoWrapper> op1_marto{{op1_repr}}, op2_marto{{op2_repr}};
	softfloat_roundingMode = sf_rnd_mode;
	auto sum_sf = f16_add(op1_sf, op2_sf);
	auto sum_marto = ieee_add_sub_impl(op1_marto, op2_marto, marto_round_mode);
	uint16_t res_repr_sf = sum_sf.v;
	uint16_t res_repr_marto = sum_marto.unravel();
	//cout << "SF: " << res_repr_sf << endl << "Marto :" << res_repr_marto << endl;
	return (isNan(sum_sf) and sum_marto.isNaN().unravel()) or (res_repr_sf == res_repr_marto);
}

BOOST_AUTO_TEST_CASE(TestCaseRndUp) {
	BOOST_REQUIRE(test_f16_sum(448, 1601, softfloat_round_max, IEEERoundingMode::RoundUp));
	BOOST_REQUIRE(test_f16_sum(448, 40964, softfloat_round_max, IEEERoundingMode::RoundUp));
	BOOST_REQUIRE(test_f16_sum(53696, 64511, softfloat_round_max, IEEERoundingMode::RoundUp));
	BOOST_REQUIRE(test_f16_sum(448, 64513, softfloat_round_max, IEEERoundingMode::RoundUp));
}

BOOST_AUTO_TEST_CASE(TestCaseRndDown) {
	BOOST_REQUIRE(test_f16_sum(448, 33216, softfloat_round_min, IEEERoundingMode::RoundDown));
	BOOST_REQUIRE(test_f16_sum(0, 32768, softfloat_round_min, IEEERoundingMode::RoundDown));
	BOOST_REQUIRE(test_f16_sum(256, 31745, softfloat_round_min, IEEERoundingMode::RoundDown));
}

BOOST_AUTO_TEST_CASE(TestCaseRndZero) {
	BOOST_REQUIRE(test_f16_sum(448, 40965, softfloat_round_minMag, IEEERoundingMode::RoundTowardZero));
	BOOST_REQUIRE(test_f16_sum(20928, 31743, softfloat_round_minMag, IEEERoundingMode::RoundTowardZero));
	BOOST_REQUIRE(test_f16_sum(320, 31745, softfloat_round_minMag, IEEERoundingMode::RoundTowardZero));
}

BOOST_AUTO_TEST_CASE(TestCaseRndTieAway) {
	BOOST_REQUIRE(test_f16_sum(0, 1, softfloat_round_near_maxMag, IEEERoundingMode::RoundNearestTieAway));
	BOOST_REQUIRE(test_f16_sum(448, 40964, softfloat_round_near_maxMag, IEEERoundingMode::RoundNearestTieAway));
	BOOST_REQUIRE(test_f16_sum(448, 1601, softfloat_round_near_maxMag, IEEERoundingMode::RoundNearestTieAway));
}

BOOST_AUTO_TEST_CASE(TestCaseRndTieEven) {
	BOOST_REQUIRE(test_f16_sum(128, 32768, softfloat_round_near_even, IEEERoundingMode::RoundNearestTieEven));
}

void compute_ieee_sum(decltype(softfloat_roundingMode) sf_rnd_mode, IEEERoundingMode marto_round_mode)
{
	constexpr unsigned int WE = 5;
	constexpr unsigned int WF = 10;

	softfloat_roundingMode = sf_rnd_mode;

	using MartoIEEE = IEEENumber<WE, WF, hint::VivadoWrapper>;
	constexpr uint64_t FORMAT_SIZE = 1 + WE + WF;
	constexpr uint64_t FORMAT_LIMIT = 1 << FORMAT_SIZE;
	uint64_t global_counter = 0;
	uint64_t counter = 0;

	int keep_going = -1;
	SumError retval = {0, 0, 0, 0, SumError::Code::OK};

	#pragma omp parallel for private(counter) shared(global_counter) schedule(static, 64)
	for (uint64_t count1 = 0 ; count1 < FORMAT_LIMIT; ++count1) {
		uint16_t op1_repr = count1;
		float16_t op1_sf{op1_repr};
		MartoIEEE op1_marto{{op1_repr}};

		for (uint32_t count2=count1 ; count2 < FORMAT_LIMIT && keep_going < 0; count2++ ) {
			uint16_t op2_repr = count2;
			float16_t op2_sf{op2_repr};
			MartoIEEE op2_marto{{op2_repr}};
			auto sum_sf = f16_add(op1_sf, op2_sf);
			auto sum_marto = ieee_add_sub_impl(op1_marto, op2_marto, marto_round_mode);
			if (isNan(sum_sf)) {//result is NaN
				bool marto_is_nan = (sum_marto.isNaN().unravel() == 1);
				if(not marto_is_nan){
					if (keep_going < 0) {
						#pragma omp critical (res)
						{
							retval.op1 = op1_repr;
							retval.op2 = op2_repr;
							retval.expected_result = sum_sf.v;
							retval.result = sum_marto.unravel();
							retval.err_code = SumError::Code::WaitingNaN;
							keep_going = 1;
						}
					}
				}
			}
			else{
				bool must = (sum_sf.v == sum_marto.unravel());
				// bool must = ((marto_prod.unravel()-resgmp) <= 1) and ((resgmp-marto_prod.unravel()) <= 1);

				if(not must){
					if (keep_going < 0) {
						#pragma omp critical (res)
						{
							retval.op1 = op1_repr;
							retval.op2 = op2_repr;
							retval.expected_result = sum_sf.v;
							retval.result = sum_marto.unravel();
							retval.err_code = SumError::Code::ResDiffer;
							keep_going = 1;
						}
					}
				}
			}
		}
		counter++;
		if((counter % PRINT_EVERY) == 0){
			#pragma omp critical (print)
			{
				global_counter += PRINT_EVERY;
				fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)", static_cast<double>(global_counter)/static_cast<double>(FORMAT_LIMIT)*100, global_counter,FORMAT_LIMIT);
			}
		}
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%lu\t/%lu)\n", static_cast<double>(FORMAT_LIMIT)/static_cast<double>(FORMAT_LIMIT)*100, FORMAT_LIMIT,FORMAT_LIMIT);
	BOOST_REQUIRE_MESSAGE(retval.err_code == SumError::Code::OK, "Sum error on operands " << retval.op1 << " and " << retval.op2);
};

BOOST_AUTO_TEST_CASE(TestIEEEAdd_5_10_SP_RNTE, *utf::disabled() * utf::label("long"))
{
	compute_ieee_sum(softfloat_round_near_even, IEEERoundingMode::RoundNearestTieEven);
}

BOOST_AUTO_TEST_CASE(TestIEEEAdd_5_10_SP_RNUp, *utf::disabled() * utf::label("long"))
{
	compute_ieee_sum(softfloat_round_max, IEEERoundingMode::RoundUp);
}

BOOST_AUTO_TEST_CASE(TestIEEEAdd_5_10_SP_RNDown, *utf::disabled() * utf::label("long"))
{
	compute_ieee_sum(softfloat_round_min, IEEERoundingMode::RoundDown);
}

BOOST_AUTO_TEST_CASE(TestIEEEAdd_5_10_SP_RNDZero, *utf::disabled() * utf::label("long"))
{
	compute_ieee_sum(softfloat_round_minMag, IEEERoundingMode::RoundTowardZero);
}

BOOST_AUTO_TEST_CASE(TestIEEEAdd_5_10_SP_RNTA, *utf::disabled() * utf::label("long"))
{
	compute_ieee_sum(softfloat_round_near_maxMag, IEEERoundingMode::RoundNearestTieAway);
}
#endif

bool test_ieee_4_7_val(u_int32_t in0, uint32_t in1)
{
	constexpr unsigned int WE = 4;
	constexpr unsigned int WF = 7;
	using SIEEE = libnumform::SmallIEEENumber<WE, WF>;

	IEEENumber<WE, WF, hint::VivadoWrapper> var0ieee{{in0}};
	IEEENumber<WE, WF, hint::VivadoWrapper> var1ieee{{in1}};

	SIEEE op0{in0}, op1{in1};
	auto sum = op0 + op1;
	auto sum_repr = sum.getRepr();

	auto marto_sum = ieee_add_sub_impl(var0ieee, var1ieee);
	bool res;
	if (marto_sum.isNaN().unravel() == 1){
		res = sum.isNaN();
	} else {
		res = marto_sum.unravel() == sum_repr;
	}
	if (!res)
		cerr << "Expecting " << sum_repr << " got " << marto_sum.unravel() << endl;
	return res;
}

BOOST_AUTO_TEST_CASE(TestSomeAdd4_7)
{
	BOOST_REQUIRE(test_ieee_4_7_val(255, 2304));
}

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
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%u\t/%u)", static_cast<double>(counter)/static_cast<double>(FORMAT_LIMIT)*100, counter,FORMAT_LIMIT);
		}
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%u\t/%u)\n", static_cast<double>(FORMAT_LIMIT)/static_cast<double>(FORMAT_LIMIT)*100, FORMAT_LIMIT,FORMAT_LIMIT);

}//*/

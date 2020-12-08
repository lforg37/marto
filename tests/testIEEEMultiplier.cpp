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
namespace utf = boost::unit_test;
#include <bitset>

#ifdef SOFTFLOAT
extern "C" {
#include "softfloat.h"
}

struct ProdError {
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

inline bool isNan(float16_t val)
{
    constexpr uint16_t nan_inf_flag = ((1 << 5) - 1) << 10;
    constexpr uint16_t absval_flag = (1 << 15) - 1;
    uint16_t repr = val.v & absval_flag;
    return ((repr & nan_inf_flag) == nan_inf_flag) && (repr ^ nan_inf_flag);
}

void compute_ieee_mult()
{
    constexpr unsigned int WE = 5;
    constexpr unsigned int WF = 10;

    auto softfloat_roundingMode = softfloat_round_near_even;

    using MartoIEEE = IEEENumber<WE, WF, hint::VivadoWrapper>;
    constexpr uint64_t FORMAT_SIZE = 1 + WE + WF;
    constexpr uint64_t FORMAT_LIMIT = 1 << FORMAT_SIZE;
    uint64_t global_counter = 0;
    uint64_t counter = 0;
    constexpr uint64_t PRINT_EVERY = 10;

    int keep_going = -1;
    ProdError retval = {0, 0, 0, 0, ProdError::Code::OK};

    #pragma omp parallel for private(counter) shared(global_counter) schedule(static, 64)
    for (uint64_t count1 = 0 ; count1 < FORMAT_LIMIT; ++count1) {
        uint16_t op1_repr = count1;
        float16_t op1_sf{op1_repr};
        MartoIEEE op1_marto{{op1_repr}};

        for (uint32_t count2=count1 ; count2 < FORMAT_LIMIT && keep_going < 0; count2++ ) {
            uint16_t op2_repr = count2;
            float16_t op2_sf{op2_repr};
            MartoIEEE op2_marto{{op2_repr}};
            auto prod_sf = f16_mul(op1_sf, op2_sf);
            auto prod_repr = prod_sf.v;
            auto prod_marto = ieee_product(op1_marto, op2_marto);
            if (isNan(prod_sf)) {//result is NaN
                bool marto_is_nan = (prod_marto.isNaN().unravel() == 1);
                if(not marto_is_nan){
                    if (keep_going < 0) {
                        #pragma omp critical (res)
                        {
                            retval.op1 = op1_repr;
                            retval.op2 = op2_repr;
                            retval.expected_result = prod_sf.v;
                            retval.result = prod_marto.unravel();
                            retval.err_code = ProdError::Code::WaitingNaN;
                            keep_going = 1;
                        }
                    }
                }
            }
            else{
                bool must = (prod_sf.v == prod_marto.unravel());
                // bool must = ((marto_prod.unravel()-resgmp) <= 1) and ((resgmp-marto_prod.unravel()) <= 1);

                if(not must){
                    if (keep_going < 0) {
                        #pragma omp critical (res)
                        {
                            retval.op1 = op1_repr;
                            retval.op2 = op2_repr;
                            retval.expected_result = prod_sf.v;
                            retval.result = prod_marto.unravel();
                            retval.err_code = ProdError::Code::ResDiffer;
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
    BOOST_REQUIRE_MESSAGE(retval.err_code == ProdError::Code::OK, "Sum error on operands " << retval.op1 << " and " << retval.op2);
};

BOOST_AUTO_TEST_CASE(TestIEEEMult_5_10_SF, *utf::disabled() * utf::label("long"))
{
    compute_ieee_mult();
}


#endif

BOOST_AUTO_TEST_CASE(TestIEEEMul0)
{
	constexpr unsigned int WE = 3;
	constexpr unsigned int WF = 4;

	using SIEEE = libnumform::SmallIEEENumber<WE, WF>;

	SIEEE a{0b00000011}, b{0b00001011};

	auto c = a*b;

	BOOST_REQUIRE_MESSAGE(c.getRepr()==1, "Error");
}


BOOST_AUTO_TEST_CASE(TestIEEMul_4_7)
{
	constexpr unsigned int WE = 4;
	constexpr unsigned int WF = 7;

	using SIEEE = libnumform::SmallIEEENumber<WE, WF>;


	// uint16_t var1_ui16 = 0b000001001;
	// uint16_t var2_ui16 = 0b001010111;
	// SIEEE op1{var1_ui16}, op2{var2_ui16};
	// auto prod = op1 * op2;
	// auto prod_repr = prod.getRepr();
	// std::bitset<9> x1(var1_ui16);
	// std::bitset<9> x2(var2_ui16);

	// hint::VivadoWrapper<9, false> resgmp{prod.getRepr()};
	// IEEENumber<WE, WF, hint::VivadoWrapper> ieee_res_gmp{resgmp};

	// hint::VivadoWrapper<9, false> var1gmp{var1_ui16};
	// hint::VivadoWrapper<9, false> var2gmp{var2_ui16};
	// IEEENumber<WE, WF, hint::VivadoWrapper> var1ieee{var1gmp};
	// IEEENumber<WE, WF, hint::VivadoWrapper> var2ieee{var2gmp};
	// auto marto_prod = ieee_product(var1ieee, var2ieee);

	// cerr << "SIEEE \t" << hint::to_string(static_cast<hint::VivadoWrapper<9, false> >(ieee_res_gmp) ) << endl;
	// cerr << "marto \t" << hint::to_string(static_cast<hint::VivadoWrapper<9, false> >(marto_prod) ) << endl;

	// exit(0);


	// omp_set_num_threads(16);
	constexpr unsigned int FORMAT_SIZE = 1 + WE + WF;
	constexpr uint32_t FORMAT_LIMIT = 1 << FORMAT_SIZE;
	uint32_t global_counter = 0;
	uint32_t counter = 0;
	#pragma omp parallel for private(counter)
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
			auto prod = op1 * op2;
			auto prod_repr = prod.getRepr();

			auto marto_prod = ieee_product(var1ieee, var2ieee);
			if (prod.isNaN()){//result is NaN
				bool marto_is_nan = (marto_prod.isNaN().unravel() == 1);
				if(not marto_is_nan){
					#pragma omp critical
					BOOST_REQUIRE_MESSAGE(false, "NAN CASE Error for \t" <<
											hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(var1ieee) ) <<
											" and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(var2ieee) ));
				}
			}
			else{
				bool must = (prod.getRepr()==marto_prod.unravel());
				// bool must = ((marto_prod.unravel()-resgmp) <= 1) and ((resgmp-marto_prod.unravel()) <= 1);
				if(not must){
					#pragma omp critical
					BOOST_REQUIRE_MESSAGE(false, "Error for " <<
										  hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(var1ieee) ) <<
										  " and " << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(var2ieee) ) <<
										  "\nexpecting\t" << hint::to_string(VivadoWrapper<FORMAT_SIZE, false>{prod_repr} ) <<
										  "\ngot  \t" << hint::to_string(static_cast<VivadoWrapper<FORMAT_SIZE, false> >(marto_prod) ));
				}
			}
		}
		counter++;
		if(counter == 10){
			#pragma omp critical
			global_counter += 10;
			counter = 0;

			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%u\t/%u)", static_cast<double>(global_counter)/static_cast<double>(FORMAT_LIMIT)*100, global_counter,FORMAT_LIMIT);
		}
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%u\t/%u)\n", static_cast<double>(FORMAT_LIMIT)/static_cast<double>(FORMAT_LIMIT)*100, FORMAT_LIMIT,FORMAT_LIMIT);

}//*/

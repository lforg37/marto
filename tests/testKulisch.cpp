#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE KulischAcc

#include <boost/test/unit_test.hpp>
#define AP_INT_MAX_W 4500

#include <array>
#include <stdlib.h>
#include <mpfr.h>

#include "kulisch/kulisch_dim.hpp"
#include "kulisch/fp_exact_prod.hpp"
#include "kulisch/acc_rounding.hpp"
#include "kulisch/K1.hpp"
#include "kulisch/K3.hpp"

#include "hint.hpp"
#include "tools/printing.hpp"

namespace utf = boost::unit_test;

static constexpr unsigned int N = 32;
static constexpr unsigned int SEGMENT = 32;

using std::array;

template <unsigned int W>
using Wrapper = hint::VivadoWrapper<N, false>;

using hint::VivadoWrapper;
using fp_dim = StandardIEEEDim<N>;

using ieee_t = IEEENumber<fp_dim::WE, fp_dim::WF, VivadoWrapper>;


typedef union {
  uint32_t i;
  float f;
} float_to_fix;

template<int low, int high>
inline float rand_float()
{
	return static_cast<float>(low) + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX)/static_cast<float>(high - low));
}


BOOST_AUTO_TEST_CASE(TestRand2CK1_32, *utf::disabled() * utf::label("long"))
{
	float_to_fix conv1, conv2;

	constexpr size_t SIZE = 1000000;

	mpfr_t mpfr_acc;
	mpfr_init2(mpfr_acc,1000);
	mpfr_t current_mpfr_value1, current_mpfr_value2;
	mpfr_init2(current_mpfr_value1,1000);
	mpfr_init2(current_mpfr_value2,1000);
	mpfr_set_d(mpfr_acc, 0.0, GMP_RNDZ);

	float& f_op1 = conv1.f;
	float& f_op2 = conv2.f;
	uint32_t& i_op1 = conv1.i;
	uint32_t& i_op2 = conv2.i;

	KulischAcc<fp_dim::WE, fp_dim::WF, VivadoWrapper> acc{{0}};
	float mpfr_sum;

	for(unsigned int i = 0 ; i < SIZE ; ++i){
		f_op1 = rand_float<-1, 1>();
		f_op2 = rand_float<-1, 1>();

		ieee_t m_op1{{i_op1}}, m_op2{{i_op2}};

		auto prod = exact_prod(m_op1, m_op2);
		acc = add_2CK1(acc, prod);

		mpfr_set_flt(current_mpfr_value1, f_op1, MPFR_RNDN);
		mpfr_set_flt(current_mpfr_value2, f_op2, MPFR_RNDN);
		mpfr_mul(current_mpfr_value1, current_mpfr_value1, current_mpfr_value2, MPFR_RNDN);
		mpfr_add(mpfr_acc, mpfr_acc, current_mpfr_value1, MPFR_RNDN);

		auto res = acc_IEEE_rounding(acc);
		float_to_fix conv_r;
		conv_r.i = res.unravel();
		float res_f = conv_r.f;
		mpfr_sum = mpfr_get_flt(mpfr_acc, MPFR_RNDN);

		if(mpfr_sum != res_f){
			fprintf(stderr, "Iteration %d\n", i);
			fprintf(stderr, "op_1 : %1.15f\n", f_op1);
			fprintf(stderr, "op_2 : %1.15f\n", f_op2);
			fprintf(stderr, "Result : %1.15f\n", res_f);
			fprintf(stderr, "mpfr : %1.15f\n", mpfr_sum);
			conv_r.f = mpfr_sum;
			VivadoWrapper<32, false> mpfr_bits{{conv_r.i}};
			cerr << hint::to_string(mpfr_bits) << endl;
			BOOST_REQUIRE_MESSAGE(false, "Accumulator holding incorrect result");
		}
		// fprintf(stderr, "%d\n", i);
	}
}

BOOST_AUTO_TEST_CASE(TestRand2CK3_32, *utf::disabled() * utf::label("long"))
{
	constexpr unsigned int SIZE = 1000000;

	mpfr_t mpfr_acc;
	mpfr_init2(mpfr_acc,1000);
	mpfr_t current_mpfr_value1, current_mpfr_value2;
	mpfr_init2(current_mpfr_value1,1000);
	mpfr_init2(current_mpfr_value2,1000);
	mpfr_set_d(mpfr_acc, 0.0, MPFR_RNDN);

	float_to_fix conv1, conv2;
	float& f_op1 = conv1.f;
	float& f_op2 = conv2.f;
	uint32_t& i_op1 = conv1.i;
	uint32_t& i_op2 = conv2.i;

	acc_2CK3<fp_dim::WE, fp_dim::WF, SEGMENT, VivadoWrapper> acc{{{0}}};
	KulischAcc<fp_dim::WE, fp_dim::WF, VivadoWrapper> acc2{{0}};
	float mpfr_sum;
	for (unsigned int i = 0 ; i < SIZE ; ++i) {
		f_op1 = rand_float<-1, 1>();
		f_op2 = rand_float<-1, 1>();

		ieee_t m_op1{{i_op1}}, m_op2{{i_op2}};

		auto prod = exact_prod(m_op1, m_op2);
		acc = acc.add_2CK3(prod);

		acc2 = add_2CK1(acc2, prod);
		auto accprop = acc.propagate_carries_2CK3();


		mpfr_set_flt(current_mpfr_value1, f_op1, MPFR_RNDN);
		mpfr_set_flt(current_mpfr_value2, f_op2, MPFR_RNDN);
		mpfr_mul(current_mpfr_value1, current_mpfr_value1, current_mpfr_value2, MPFR_RNDN);
		mpfr_add(mpfr_acc, mpfr_acc, current_mpfr_value1, MPFR_RNDN);

		auto resbis = acc_IEEE_rounding(acc2);//TODO remove;
		auto res = acc_IEEE_rounding(accprop);
		float_to_fix conv_r;
		conv_r.i = res.unravel();
		float res_f = conv_r.f;
		mpfr_sum = mpfr_get_flt(mpfr_acc, MPFR_RNDN);

		if(mpfr_sum != res_f){
			fprintf(stderr, "Result : %1.15f\n", res_f);
			fprintf(stderr, "mpfr : %1.15f\n", mpfr_sum);

			BOOST_REQUIRE_MESSAGE(false, "Accumulator holding incorrect result");
		}

		if(((i%1000) == 0) and (i != 0)){
			fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%u\t/%u)", static_cast<double>(i)/static_cast<double>(SIZE)*100, i,SIZE);
		}
		// fprintf(stderr, "%d\n", i);
	}
	fprintf(stderr, "\33[2K\rCompletion: \t%1.1f%% (%u\t/%u)\n", static_cast<double>(SIZE)/static_cast<double>(SIZE)*100, SIZE,SIZE);
}

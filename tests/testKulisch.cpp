#define BOOST_TEST_DYN_LINK   
#define BOOST_TEST_MODULE KulischAcc

#include <boost/test/unit_test.hpp>


#include "../includes/kulisch_acc.hpp"
#include "../includes/kulisch_dim.hpp"
#include "../includes/utils.hpp"
#include <stdlib.h>
#include <mpfr.h>

namespace utf = boost::unit_test;


#define N 32
#define SEGMENT 64

typedef union {
  unsigned int i;
  float f;
} float_to_fix;

BOOST_AUTO_TEST_CASE(TestRandKulischAcc32, *utf::disabled() * utf::label("long")) 
{
	int SIZE = 10000000;
	float *tab1, *tab2, sign1, sign2;
	ap_uint<N> *tab1_ap, *tab2_ap;

	tab1 = (float*) malloc(sizeof(float)*SIZE);
	tab2 = (float*) malloc(sizeof(float)*SIZE);
	tab1_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);
	tab2_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);

	for(int i=0; i<SIZE; i++){
		tab1[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		tab2[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		sign1 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		sign2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		if(((int) sign1) == 0){
			tab1[i] = -tab1[i]; 
		}
		if(((int) sign2) == 0){
			tab2[i] = -tab2[i]; 
		}

	}

	// tab1[0] = 293.992; 
	// tab2[0] = 0.19882; 


	for(int i=0; i<SIZE; i++){
		float_to_fix conv;
		conv.f = tab1[i] ;
		tab1_ap[i] = conv.i;
		conv.f = tab2[i] ;
		tab2_ap[i] = conv.i;
	};


	mpfr_t mpfr_acc;
	mpfr_init2(mpfr_acc,1000);
	mpfr_t current_mpfr_value1, current_mpfr_value2;
	mpfr_init2(current_mpfr_value1,1000);
	mpfr_init2(current_mpfr_value2,1000);
 	mpfr_set_d(mpfr_acc, 0.0, GMP_RNDZ);
		




	KulischAcc<N> acc(0);
	ap_uint<32> res;
	float mpfr_sum;
	fprintf(stderr, "%d\n", SIZE);
	for(int i=0; i<SIZE; i++){
		// fprintf(stderr, "tab1 : %f, tab2: %f\n", tab1[i], tab2[i]);
		auto prod = exact_prod<N>(tab1_ap[i], tab2_ap[i]);
		acc = kulisch_accumulator<N>(acc, prod);
		
		mpfr_set_flt(current_mpfr_value1, tab1[i], GMP_RNDZ);
		mpfr_set_flt(current_mpfr_value2, tab2[i], GMP_RNDZ);
		mpfr_mul(current_mpfr_value1, current_mpfr_value1, current_mpfr_value2, GMP_RNDZ);
		mpfr_add(mpfr_acc, mpfr_acc, current_mpfr_value1, GMP_RNDZ);
		
		res = acc_to_fp<N>(acc);
		float_to_fix conv_r;
		conv_r.i = res;
		float res_f = conv_r.f;
		mpfr_sum = mpfr_get_flt(mpfr_acc, MPFR_RNDN);
		
		if(mpfr_sum != res_f){
			fprintf(stderr, "Result : %1.15f\n", res_f);
			printApUint(res);
			fprintf(stderr, "mpfr : %1.15f\n", mpfr_sum);
			conv_r.f = mpfr_sum;
			ap_uint<32> mpfr_bits = conv_r.i;
			printApUint(mpfr_bits);

			free(tab1);
			free(tab2);
			free(tab1_ap);
			free(tab2_ap);
			BOOST_REQUIRE_MESSAGE(false, "Accumulator holding incorrect result");
		}
		// fprintf(stderr, "%d\n", i);
	}	

	free(tab1);
	free(tab2);
	free(tab1_ap);
	free(tab2_ap);
}

BOOST_AUTO_TEST_CASE(TestRandSegmentedKulischAcc32, *utf::disabled() * utf::label("long")) 
{
	int SIZE = 1000000000;
	float *tab1, *tab2, sign1, sign2;
	ap_uint<N> *tab1_ap, *tab2_ap;

	tab1 = (float*) malloc(sizeof(float)*SIZE);
	tab2 = (float*) malloc(sizeof(float)*SIZE);
	tab1_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);
	tab2_ap = (ap_uint<32>*) malloc(sizeof(ap_uint<32>)*SIZE);

	for(int i=0; i<SIZE; i++){
		tab1[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		tab2[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
		sign1 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		sign2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/2));
		if(((int) sign1) == 0){
			tab1[i] = -tab1[i]; 
		}
		if(((int) sign2) == 0){
			tab2[i] = -tab2[i]; 
		}

	}

	// tab1[0] = 293.992; 
	// tab2[0] = 0.19882; 


	for(int i=0; i<SIZE; i++){
		float_to_fix conv;
		conv.f = tab1[i] ;
		tab1_ap[i] = conv.i;
		conv.f = tab2[i] ;
		tab2_ap[i] = conv.i;
	};


	mpfr_t mpfr_acc;
	mpfr_init2(mpfr_acc,1000);
	mpfr_t current_mpfr_value1, current_mpfr_value2;
	mpfr_init2(current_mpfr_value1,1000);
	mpfr_init2(current_mpfr_value2,1000);
 	mpfr_set_d(mpfr_acc, 0.0, GMP_RNDZ);
		




	SegmentedKulischAcc<N, SEGMENT> acc(KulischAcc<N>(0));
	KulischAcc<N> acc2(0); 
	ap_uint<32> res;
	float mpfr_sum;
	fprintf(stderr, "%d\n", SIZE);
	for(int i=0; i<SIZE; i++){
		// fprintf(stderr, "tab1 : %f, tab2: %f\n", tab1[i], tab2[i]);
		auto prod = exact_prod<N>(tab1_ap[i], tab2_ap[i]);
		acc = segmented_kulisch_accumulator<N, SEGMENT>(acc, prod);
		acc2 = kulisch_accumulator<N>(acc2, prod);
		// acc.printContent();
		KulischAcc<N> accprop = Kulisch_propagate_carries(acc);
		
		// printApUint(accprop);
		// printApUint(acc2);

		mpfr_set_flt(current_mpfr_value1, tab1[i], GMP_RNDZ);
		mpfr_set_flt(current_mpfr_value2, tab2[i], GMP_RNDZ);
		mpfr_mul(current_mpfr_value1, current_mpfr_value1, current_mpfr_value2, GMP_RNDZ);
		mpfr_add(mpfr_acc, mpfr_acc, current_mpfr_value1, GMP_RNDZ);
		
		res = acc_to_fp<N>(accprop);
		float_to_fix conv_r;
		conv_r.i = res;
		float res_f = conv_r.f;
		mpfr_sum = mpfr_get_flt(mpfr_acc, MPFR_RNDN);
		
		if(mpfr_sum != res_f){
			fprintf(stderr, "Result : %1.15f\n", res_f);
			printApUint(res);
			fprintf(stderr, "mpfr : %1.15f\n", mpfr_sum);
			conv_r.f = mpfr_sum;
			ap_uint<32> mpfr_bits = conv_r.i;
			printApUint(mpfr_bits);

			free(tab1);
			free(tab2);
			free(tab1_ap);
			free(tab2_ap);
			BOOST_REQUIRE_MESSAGE(false, "Accumulator holding incorrect result");
		}
		// fprintf(stderr, "%d\n", i);
	}	

	free(tab1);
	free(tab2);
	free(tab1_ap);
	free(tab2_ap);
}
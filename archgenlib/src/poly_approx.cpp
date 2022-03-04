#include "runtime/poly_approx.hpp"
#include "runtime/sollya_fix_function.hpp"
#include "runtime/sollya_handler.hpp"
#include <sollya.h>
#include <vector>


namespace archgenlib{

	BasicPolyApprox::BasicPolyApprox(SollyaFunction sf, int degree_, int width)
	{
        sollya_obj_t fS = sf.function; // no need to free this one
		sollya_obj_t inputRangeS = sf.domain; // no need to free this one
		SollyaHandler degreeS {sollya_lib_constant_from_int(degree)};

        std::vector<SollyaHandler> coeffsHandler;
        std::vector<sollya_obj_t> coeffsSollya;
        coeffsHandler.reserve(degree_ + 1);
        coeffsSollya.reserve(degree_ + 1);
        for(int i=0; i<=degree ; i++) {
            coeffsSollya.emplace_back(sollya_lib_constant_from_int(width));
            coeffsHandler.emplace_back(coeffsSollya[i]);
		}

        SollyaHandler coeffSizeListS {sollya_lib_list(coeffsSollya.data(), degree_ + 1)};
		SollyaHandler fixed {sollya_lib_fixed()};
        SollyaHandler absolute{sollya_lib_absolute()};
        
        auto fixedS = static_cast<sollya_obj_t>(fixed);
        auto absoluteS = static_cast<sollya_obj_t>(absolute);
        
        // Tadaaa! After all this we may launch fpminimax
		polynomialS = sollya_lib_fpminimax(fS, degreeS, coeffSizeListS, inputRangeS, fixedS, absoluteS, NULL);
		/*/ TODO
		// Checking its approximation error;
		sollya_obj_t supNormS; // it will end up there
		sollya_obj_t supNormAccS = sollya_lib_parse_string("1b-10"); // This is the size of the returned interval... 10^-3 should be enough for anybody
		sollya_obj_t supNormInputRangeS = sollya_lib_supnorm(polynomialS, fS, inputRangeS, absoluteS, supNormAccS);
		if(sollya_lib_obj_is_error(supNormInputRangeS)) {
			cout <<  ">   Sollya infnorm failed, but do not loose all hope yet: launching dirtyinfnorm:" << endl;
			sollya_obj_t pminusfS = sollya_lib_sub(polynomialS, fS);
			supNormS = sollya_lib_dirtyinfnorm(pminusfS, inputRangeS);
			sollya_lib_clear_obj(pminusfS);
		}
		else{ // supnorm succeeded, we are mostly interested in the sup of the interval
			supNormS = sollya_lib_sup(supNormInputRangeS);
		}
		sollya_lib_clear_obj(supNormAccS);
		sollya_lib_clear_obj(supNormInputRangeS);

		sollya_lib_get_constant_as_double(& approxErrorBound, supNormS);
		sollya_lib_clear_obj(supNormS);

		REPORT(DETAILED, "Polynomial accuracy is " << approxErrorBound);
		// Please leave the memory in the state you would like to find it when entering
		sollya_lib_clear_obj(degreeS);



		buildFixFormatVector();*/
	}


	// This is a static (class) method.
	void BasicPolyApprox::guessDegree(sollya_obj_t fS, sollya_obj_t inputRangeS, sollya_obj_t targetAccuracy, int& degreeInf, int& degreeSup) {
		// initial evaluation of the required degree
		// guessdegree function prototype
		// 	guessdegree(f,I,eps,w,bound ) : (function, range, constant, function, constant) → range
		// guessdegree function parameters
		// 	f is the function to be approximated.
	    // 	I is the interval where the function must be approximated.
	    // 	eps is the maximal acceptable error.
	    // 	w (optional) is a weight function. Default is 1.
	    // 	bound (optional) is a bound on the degree. Default is currently 128.
		// guessdegree function reyurn value
		// 	returns an interval: for common cases, this interval is reduced to a single number (i.e. the minimal degree).
		// 	But in certain cases, guessdegree does not succeed in finding the minimal degree.
		// 	In such cases the returned interval is of the form [n,p] such that:
		SollyaHandler degreeIntervalS{sollya_lib_guessdegree(fS, inputRangeS, targetAccuracy, NULL)};
        SollyaHandler degreeInfS{sollya_lib_inf(degreeIntervalS)};
		SollyaHandler degreeSupS{sollya_lib_sup(degreeIntervalS)};
		sollya_lib_get_constant_as_int(&degreeInf, degreeInfS);
		sollya_lib_get_constant_as_int(&degreeSup, degreeSupS);
    }

	void BasicPolyApprox::buildFixFormatVector()
	{
        #if 0
        //TODO
		// compute the MSBs
		int msb,lsb;
		for (int i=0; i<=degree; i++){
			sollya_obj_t iS = sollya_lib_constant_from_int(i);
			//sollya_lib_printf("i = %d  = %b\n", i, iS);
			sollya_obj_t coeffS = sollya_lib_coeff(polynomialS, iS);
			//sollya_lib_printf(">  c%d = %b \n", i, coeffS);
			sollya_lib_clear_obj(iS);

			mpfr_t mpcoeff, mptmp;
			// First a tentative conversion to double to get an estimate of the MSB and zeroness
			// FIXME do everything in MPFR before it kills somebody
			double dcoeff;
			sollya_lib_get_constant_as_double(&dcoeff, coeffS);
			if(0.0==dcoeff) {
				msb=f->msbOut; // for the sign
				lsb=f->msbOut-1;
				mpfr_init2(mpcoeff, 32);
				mpfr_set_d(mpcoeff, 0.0, GMP_RNDN);
			}
			else{
				msb = floor(log2(fabs(dcoeff)));  // ex. 2.01
				msb++; // For the sign
				lsb = LSB;
				// now we may safely allocate the proper size for the mpfr_t. Add two bits for sign + rounding.
				mpfr_init2(mpcoeff, msb-lsb+3);
				mpfr_init2(mptmp, msb-lsb+3);
				sollya_lib_get_constant(mpcoeff, coeffS);
				// Now recompute the MSB explicitely.
				mpfr_abs(mptmp, mpcoeff, GMP_RNDN); // exact
				mpfr_log2(mptmp, mptmp, GMP_RNDU);
				mpfr_floor(mptmp, mptmp);
				msb = mpfr_get_si(mptmp, GMP_RNDU);
				mpfr_clear(mptmp);
				msb++;
			}
			FixConstant* fixcoeff =	new FixConstant(msb, lsb, true/*signed*/, mpcoeff);
			coeff.push_back(fixcoeff);

			sollya_lib_clear_obj(coeffS);
			mpfr_clear(mpcoeff);
		}
        #endif
	}	
} //namespace


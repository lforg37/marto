#include <cstddef>
#include <cstdint>
#include <gmpxx.h>
#include <sollya.h>

#include "fixedpoint/fixedpoint.hpp"
#include "runtime/mpfr_handler.hpp"
#include "runtime/sollya_fix_function.hpp"
#include "runtime/sollya_handler.hpp"
#include "runtime/utility.hpp"

namespace{
sollya_obj_t compute_output_domain(sollya_obj_t function, sollya_obj_t input_domain, int prec) {
 sollya_lib_set_prec(sollya_lib_constant_from_int(prec));
 return sollya_lib_evaluate(function, input_domain);
};
}


namespace archgenlib {
SollyaFunction::SollyaFunction(SollyaHandler& f, FPDimRTRepr input_f):function{f}, input_format{input_f},domain{sollya_interval_from_rtfpdim(input_f)}, output_domain{compute_output_domain(f, domain, 100*input_f.width)}{
  SollyaHandler sup{sollya_lib_sup(output_domain)}, inf{sollya_lib_inf(output_domain)};
  MPFRHandler supMPFR{100*input_f.width}, infMPFR{100*input_f.width};
  sollya_lib_get_constant(supMPFR, sup);
  sollya_lib_get_constant(infMPFR, inf);
  signed_output = mpfr_sgn(static_cast<mpfr_t&>(infMPFR)) < 0;
  mpfr_abs(supMPFR, supMPFR, MPFR_RNDU);
  mpfr_abs(infMPFR, infMPFR, MPFR_RNDU);
  mpfr_max(supMPFR, infMPFR, supMPFR, MPFR_RNDU);
  mpfr_log2(supMPFR, supMPFR, MPFR_RNDU);
  mpfr_floor(supMPFR, supMPFR);
  msb_output = mpfr_get_si(supMPFR, MPFR_RNDU);
}

mpz_class SollyaFunction::getmpz(mpfr_t evaluationpoint, bitweight_t LSBOut, mpfr_t intermResHolder) {
    mpz_class ret;
    sollya_lib_evaluate_function_at_point(intermResHolder, function, evaluationpoint, nullptr);
    mpfr_mul_2si(intermResHolder,intermResHolder, -LSBOut, MPFR_RNDN);
    mpfr_get_z(ret.get_mpz_t(), intermResHolder, MPFR_RNDN);
    return ret;
}

std::vector<mpz_class> SollyaFunction::faithful_at_weight(bitweight_t LSBOut) {
    assert(LSBOut <= msb_output);
    std::vector<mpz_class> ret;
    std::uint64_t size{1};
    size <<= input_format.width;
    ret.reserve(size);

    // Handle positive values 
    vecwidth_t precision = msb_output - LSBOut + 1;
    MPFRHandler curval{input_format.width}, increment{input_format.width}, output{precision};

    mpfr_set_uj(curval, 0, MPFR_RNDN);
    mpfr_set_sj_2exp(increment, 1, input_format.lsb_weight, MPFR_RNDN);
    uint64_t limit = (input_format.is_signed) ? size >> 1 : size;
    sollya_lib_set_prec(sollya_lib_constant_from_int(precision + 1));
    for (uint64_t i = 0 ; i < limit ; ++i) {
      ret.emplace_back(getmpz(curval, LSBOut, output));
      mpfr_add(curval, curval, increment, MPFR_RNDN);
    }
    if (input_format.is_signed) { // Should now handle the negative inputs
      mpfr_set_sj_2exp(curval, -1, input_format.msb_weight, MPFR_RNDN);
      for (uint64_t i = limit ; i < size ; ++i) {
        ret.emplace_back(getmpz(curval, LSBOut, output));
        mpfr_add(curval, curval, increment, MPFR_RNDN);
      }
    }
    return ret;
};
}

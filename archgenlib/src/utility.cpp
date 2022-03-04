#include <gmpxx.h>
#include <mpfr.h>
#include <sollya.h>

#include "runtime/sollya_handler.hpp"
#include "runtime/utility.hpp"

namespace archgenlib {
sollya_obj_t sollya_interval_from_rtfpdim(FPDimRTRepr dim) {
  mpfr_t min, max;
  mpfr_init2(min, dim.width);
  mpfr_init2(max, dim.width);
  if (dim.is_signed) {
    mpfr_set_sj_2exp(min, -1, dim.msb_weight-1, MPFR_RNDNA);
    mpfr_set_uj(max, 0, MPFR_RNDN);
    mpfr_sub(max, max, min, MPFR_RNDNA);
    mpfr_t lsb;
    mpfr_init2(lsb, dim.width);
    mpfr_set_uj_2exp(lsb, 1, dim.lsb_weight, MPFR_RNDZ);
    mpfr_sub(max, max, lsb, MPFR_RNDNA);
    mpfr_clear(lsb);
  } else {
    mpfr_set_uj_2exp(min, 1, dim.lsb_weight, MPFR_RNDZ);
    mpfr_set_uj_2exp(max, 1, dim.msb_weight, MPFR_RNDU);
    mpfr_sub(max, max, min, MPFR_RNDNA);
    mpfr_set_uj(min, 0, MPFR_RNDZ);
  }
  auto retval = sollya_lib_range_from_bounds(min, max);
  mpfr_clears(max, min, NULL);
  return retval;
}
}

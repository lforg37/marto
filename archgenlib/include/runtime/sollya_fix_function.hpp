#ifndef RUNTIME_SOLLYA_FIX_FUNCTION_HPP
#define RUNTIME_SOLLYA_FIX_FUNCTION_HPP

#include <vector>

#include <gmpxx.h>
#include <sollya.h>

#include "fixedpoint/fixedpoint.hpp"
#include "runtime/expression_tree.hpp"
#include "runtime/sollya_handler.hpp"

namespace archgenlib {
///
///@brief Non owning view of a sollya function associated with the input interval
///
struct SollyaFunction {
  SollyaHandler& function;
  SollyaHandler domain;
  const FPDimRTRepr input_format;
  SollyaHandler output_domain;
  bool signed_output;
  bitweight_t msb_output;

  ///
  ///@brief Compute the table of values of function, faithfully rounded 
  ///       to 2^LSBout       
  ///
  ///@param LSBOut 
  ///@return std::vector<mpz_class> 
  std::vector<mpz_class> faithful_at_weight(bitweight_t LSBOut);

  SollyaFunction(SollyaHandler& f, FPDimRTRepr input_f);
  
  private:
  mpz_class getmpz(mpfr_t evalPoint, bitweight_t LSBOut, mpfr_t intermResHolder);

};
}
#endif

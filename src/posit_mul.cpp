#include "posit_mul.hpp"

//Instantiation for standard sizes
template PositProd<8> posit_mul(PositValue<8> in1, PositValue<8> in2);
template PositProd<16> posit_mul(PositValue<16> in1, PositValue<16> in2);
template PositProd<32> posit_mul(PositValue<32> in1, PositValue<32> in2);
template PositProd<64> posit_mul(PositValue<64> in1, PositValue<64> in2);

#include "posit_decoder.hpp"

//Instantiation for standard posit values
template PositValue<8> posit_decoder(PositEncoding<8> positN);
template PositValue<16> posit_decoder(PositEncoding<16> positN);
template PositValue<32> posit_decoder(PositEncoding<32> positN);
template PositValue<64> posit_decoder(PositEncoding<64> positN);

#include "posit_encoder.hpp"

//Instanciation for standard posits
template PositEncoding<8> posit_encoder(PositValue<8> positValue);
template PositEncoding<16> posit_encoder(PositValue<16> positValue);
template PositEncoding<32> posit_encoder(PositValue<32> positValue);
template PositEncoding<64> posit_encoder(PositValue<64> positValue);

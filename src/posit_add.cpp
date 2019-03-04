 #include "posit_add.hpp"

//Instanciation for standard posit sizes
template PositValue<8> posit_add(PositValue<8> in1, PositValue<8> in2);
template PositValue<16> posit_add(PositValue<16> in1, PositValue<16> in2);
template PositValue<32> posit_add(PositValue<32> in1, PositValue<32> in2);
template PositValue<64> posit_add(PositValue<64> in1, PositValue<64> in2);

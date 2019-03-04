#include "quire_to_posit.hpp"

//Instantiation for standard values
template PositValue<8> quire_to_posit(Quire<8> quire);
template PositValue<16> quire_to_posit(Quire<16> quire);
template PositValue<32> quire_to_posit(Quire<32> quire);
template PositValue<64> quire_to_posit(Quire<64> quire);

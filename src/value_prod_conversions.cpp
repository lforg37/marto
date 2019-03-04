#include "value_prod_conversions.hpp"

//Instanciation for standard posit sizes
template PositProd<8> PositValue_to_PositProd(PositValue<8> val); 
template PositProd<16> PositValue_to_PositProd(PositValue<16> val); 
template PositProd<32> PositValue_to_PositProd(PositValue<32> val); 
template PositProd<64> PositValue_to_PositProd(PositValue<64> val); 

template PositValue<8> PositProd_to_PositValue(PositProd<8> val);
template PositValue<16> PositProd_to_PositValue(PositProd<16> val);
template PositValue<32> PositProd_to_PositValue(PositProd<32> val);
template PositValue<64> PositProd_to_PositValue(PositProd<64> val);

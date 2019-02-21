#ifndef STATIC_MATH_HPP
#define STATIC_MATH_HPP

template <int num, int denom>
struct Static_Ceil_Div{
	static constexpr int val = (num % denom) ? (num / denom) + 1 : num / denom;
};

#endif

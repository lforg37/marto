#ifndef STATIC_MATH_HPP
#define STATIC_MATH_HPP

template <int num, int denom>
struct Static_Ceil_Div{
	static constexpr int val = (num % denom) ? (num / denom) + 1 : num / denom;
};

constexpr int ceilLog2(int N, uint8_t remains = 0)
{
	return (N <= 1) ? remains : 1 + ceilLog2(N>>1, remains | (N%2));
}

constexpr int ceil2Power(int N)
{
	return 1 << ceilLog2(N);
}

constexpr int log2(int N)
{
	 return ((N<2) ? 1 : 1+log2(N>>1));
}

constexpr int r2pow(int N)
{
	 return 1 << log2(N);
}



template<int N>
class Static_Val
{
	public:
		static constexpr int _rlog2 = log2(N);
		static constexpr int _r2pow = r2pow(N);
		static constexpr int _log2 = ceilLog2(N);
		static constexpr int _2pow = ceil2Power(N);
};

#endif

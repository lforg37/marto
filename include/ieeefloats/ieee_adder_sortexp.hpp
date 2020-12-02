#ifndef IEEE_ADDER_SORTEXP_HPP
#define IEEE_ADDER_SORTEXP_HPP

#include "ieeefloats/ieeetype.hpp"
#include "ieeefloats/ieee_rounding.hpp"
#include "tools/static_math.hpp"
#include "primitives/lzoc.hpp"
#include "primitives/shifter_sticky.hpp"

#ifdef IEEE_ADDER_SORTEXP_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
inline IEEENumber<WE, WF, Wrapper> ieee_add_sortexp(
	   IEEENumber<WE, WF, Wrapper> in0,
	   IEEENumber<WE, WF, Wrapper> in1,
	   IEEERoundingMode roundingMode = IEEERoundingMode::RoundNearestTieEven
	)
{
	auto exp0 = in0.getExponent();
	auto exp1 = in1.getExponent();

	auto exp0IsZero = (exp0 == Wrapper<WE, false>{0});
	auto exp1IsZero = (exp1 == Wrapper<WE, false>{0});
	auto exp0IsAllOne = (exp0 == Wrapper<WE, false>{(1<<WE)-1});
	auto exp1IsAllOne = (exp1 == Wrapper<WE, false>{(1<<WE)-1});

	auto sign0 = in0.getSign();
	auto sign1 = in1.getSign();

	auto frac0 = in0.getFractionnalPart();
	auto frac1 = in1.getFractionnalPart();

	auto frac0IsZero = (frac0 == Wrapper<WF, false>{0});
	auto frac1IsZero = (frac1 == Wrapper<WF, false>{0});

	auto diff0 = exp0.modularSub(exp1);
	auto diff1 = exp1.modularSub(exp0);

	auto aligned = (exp0 == exp1);
	auto unaligned = aligned.invert();

	// Sorting according to exp
	auto swap = exp1 > exp0;
	auto keep = swap.invert();

	auto maxExp = Wrapper<WE, false>::mux(swap, exp1, exp0);
	auto minExp = Wrapper<WE, false>::mux(swap, exp0, exp1);
	auto expdiff = Wrapper<WE, false>::mux(swap, diff1, diff0);

	auto effsub = sign0 ^ sign1;

	auto tmpMaxSign = Wrapper<1, false>::mux(swap, sign1, sign0);
	auto tmpMinSign = Wrapper<1, false>::mux(swap, sign0, sign1);

	auto maxFrac = Wrapper<WF, false>::mux(swap, frac1, frac0);
	auto minFrac = Wrapper<WF, false>::mux(swap, frac0, frac1);

	auto maxExpIsZero = (swap & exp1IsZero) | (keep & exp0IsZero);
	auto minExpIsZero = (swap & exp0IsZero) | (keep & exp1IsZero);
	auto maxExpAllOne = (swap & exp1IsAllOne) | (keep & exp0IsAllOne);
	auto minExpAllOne = (swap & exp0IsAllOne) | (keep & exp1IsAllOne);
	auto maxFracIsZero = (swap & frac1IsZero) | (keep & frac0IsZero);
	auto minFracIsZero = (swap & frac0IsZero) | (keep & frac1IsZero);
	auto maxExpIsOne = (maxExp == Wrapper<WE, false>{1});

	//Special case logic
	auto maxIsInfinity = maxExpAllOne & maxFracIsZero;
	auto maxIsNaN = maxExpAllOne & maxFracIsZero.invert();
	auto maxIsZero = maxExpIsZero & maxFracIsZero;
	auto minIsInfinity = minExpAllOne & minFracIsZero;
	auto minIsNaN = minExpAllOne & minFracIsZero.invert();
	auto minIsZero = minExpIsZero & minFracIsZero;
	auto oneIsZero = maxIsZero | minIsZero;

	auto bothSubNormals = exp0IsZero & exp1IsZero;
	auto maxIsNormal = maxExpIsZero.invert();
	auto minIsNormal = minExpIsZero.invert();

	auto infinitySub = maxIsInfinity & minIsInfinity & effsub;
	auto resultIsNan = maxIsNaN | minIsNaN |infinitySub;
	auto onlyOneSubnormal = minExpIsZero & maxIsNormal;

	// Reconstruct full fraction
	auto explicitedMaxFrac = maxIsNormal.concatenate(maxFrac);
	auto explicitedMinFrac = minIsNormal.concatenate(minFrac);

	//alignment
	auto maxShiftVal = Wrapper<WE, false>{WF+3};
	auto allShiftedOut = expdiff > maxShiftVal;

	auto shiftValue = Wrapper<WE, false>::mux(
					allShiftedOut,
					maxShiftVal,
					expdiff
				).modularSub(onlyOneSubnormal.template leftpad<WE>());

	//////// WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP ///////////////////////
	auto extendedMinFrac = explicitedMinFrac.concatenate(Wrapper<2, false>{0});
	auto shiftedMinFracSticky = shifter_sticky(extendedMinFrac, shiftValue);

	////////////////////////////////////////////////////////////////////////////
	/// Branche min +- max (max soustrait seulement si eff sub et alignement)///
	///   01.xxxxxxxxxx  OR    01.xxxxxxxxxx  OR   1x.xxxxxxxxxxxx           ///
	/// + 01.xxxxxxxxxx  --  + 0.000xxxxxxxx  -- + 01.xxxxxxxxxxxx           ///
	////////////////////////////////////////////////////////////////////////////
	constexpr unsigned int AdderWidth = WF+3;
	auto shiftedMinPBranch = shiftedMinFracSticky.template slice<WF+3, 3>().template leftpad<AdderWidth>();
	auto adderShoulSub = effsub & aligned & oneIsZero.invert();
	auto adderNegMask = Wrapper<AdderWidth, false>::generateSequence(adderShoulSub);
	auto adderMaxIn = explicitedMaxFrac.template leftpad<AdderWidth>().bitwise_xor(adderNegMask);
	auto adderRes = adderMaxIn.addWithCarry(shiftedMinPBranch, adderShoulSub);
	auto adderResNeg = adderRes.template get<AdderWidth-1>();
	auto fracAdd = adderRes.template slice<WF+1, 0>().concatenate(shiftedMinFracSticky.template slice<2, 1>());
	////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	/// Branche max - min                                                    ///
	///   1.xxxxxxxxxxx00  OR    1.xxxxxxxxxxx00                             ///
	/// - 1.xxxxxxxxxxx00  --  - 0.0000xxxxxxxgs                             ///
	////////////////////////////////////////////////////////////////////////////
	constexpr unsigned int SubWidth = WF+4;
	auto shiftedMinSBranch = shiftedMinFracSticky;
	auto SubMaxIn = explicitedMaxFrac.concatenate(Wrapper<3, false>{0});
	auto subRes = SubMaxIn.modularSub(shiftedMinSBranch);
	auto fracSub = subRes.template slice<SubWidth - 1, 1>().template leftpad<WF+4>();
	////////////////////////////////////////////////////////////////////////////

	auto resIsSub = effsub & (adderResNeg | unaligned);
	auto resIsAdd = resIsSub.invert();
	auto stickyRes = shiftedMinFracSticky.template get<0>();
	auto fracRes = Wrapper<WF+4, false>::mux(resIsSub, fracSub, fracAdd);
	auto maxSign = Wrapper<1, false>::mux(effsub & aligned & adderResNeg.invert() & minIsZero.invert(), tmpMinSign, tmpMaxSign);

	auto z1 = fracRes.template get<WF+3>();
	auto z0 = fracRes.template get<WF+2>();

	auto lzc = lzoc_wrapper(fracRes, {0});

	constexpr unsigned int lzcsize = hint::Static_Val<WF+2>::_storage;

	static_assert (lzcsize<=WE, "The adder works only for wE > log2(WF).\nAre you sure you need subnormals ?\nIf yes, contact us with your use case, we will be happy to make it work for you.");
	auto subnormalLimitVal = Wrapper<lzcsize, false>{WF+3};


	auto lzcGreaterEqExp = (lzc.template leftpad<WE>() >= maxExp);
	auto lzcSmallerEqExp = (lzc.template leftpad<WE>() <= maxExp);
	auto lzcSmallerMaxVal = lzc < subnormalLimitVal;
	auto fullCancellation = lzcSmallerMaxVal.invert();

	auto normalOverflow = z1;
	auto lzcOne = z1.invert() & z0;
	auto subnormalOverflow = lzcOne & maxExpIsZero;
	auto cancellation = z1.invert() & z0.invert();

	auto overflow = normalOverflow | subnormalOverflow;

	//////// WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP WIP ///////////////////

	auto isLeftShiftLZC = overflow |
			(lzcSmallerMaxVal.invert() & bothSubNormals) |
			(cancellation & maxIsNormal & lzcSmallerEqExp) |
			(maxIsNormal & lzcSmallerMaxVal.invert());
	auto isLeftShiftExp = lzcSmallerMaxVal & lzcGreaterEqExp & maxIsNormal;

	auto shiftFirstStage = Wrapper<lzcsize, false>::mux(isLeftShiftLZC, lzc, Wrapper<lzcsize, false>{1});
	auto normalisationShiftVal = Wrapper<lzcsize, false>::mux(isLeftShiftExp, maxExp.template slice<lzcsize-1, 0>(), shiftFirstStage);

	auto normalisedSignif = fracRes << normalisationShiftVal;
	auto significandPreRound = normalisedSignif.template slice<WF+2, 3>();
	auto lsb = normalisedSignif.template get<3>();
	auto roundBit = normalisedSignif.template get<2>();
	auto sticky = stickyRes| normalisedSignif.template slice <1, 0>().or_reduction();

	auto deltaExpIsZero = z1.invert() & (z0 ^ bothSubNormals);
	auto deltaExpIsMinusOne = z1 | (z0 & bothSubNormals);
	auto deltaExpIsLZC = (z1 | z0 | bothSubNormals).invert() & lzcSmallerEqExp & lzcSmallerMaxVal;
	auto deltaExpExp = (deltaExpIsLZC | deltaExpIsZero | deltaExpIsMinusOne).invert();

	auto deltaExpCin = deltaExpExp | deltaExpIsMinusOne | deltaExpIsLZC;
	auto deltaBigPartIsZero = deltaExpIsZero | deltaExpIsMinusOne;

	auto deltaExpUnmasked = Wrapper<WE, false>::mux(
				deltaExpIsLZC,
				(lzc.modularSub(Wrapper<lzcsize, false>{1})).template leftpad<WE>().invert(),
				maxExp.invert()
			);

	auto maskSequence = Wrapper<WE, false>::generateSequence(deltaBigPartIsZero.invert());
	auto deltaExpBeforeCorrection = deltaExpUnmasked & maskSequence;
	auto expPreRound = maxExp.addWithCarry(deltaExpBeforeCorrection, deltaExpCin).template slice<WE-1, 0>();
	auto expSigPreRound = expPreRound.concatenate(significandPreRound);

	auto roundUpBit = ieee_getRoundBit(maxSign, lsb, roundBit, sticky, roundingMode);

	auto unroundedInf = expPreRound.and_reduction();

	Wrapper<3, false> roundingCode{static_cast<uint8_t>(roundingMode)};
	auto b0 = roundingCode.template get<0>();
	auto b1 = roundingCode.template get<1>();
	auto b2 = roundingCode.template get<2>();
	auto forbiddent_inf = ((b2 & b0 & (b1 == maxSign)) | (roundingCode.or_reduction().invert())) & unroundedInf & maxIsInfinity.invert() & resultIsNan.invert();

	auto expSigRounded = expSigPreRound.modularAdd(roundUpBit.template leftpad<WE+WF>());
	auto finalExp = expSigRounded.template slice<WF+WE-1, WF>();

	auto resultIsZero = fullCancellation & finalExp.or_reduction().invert();
	auto resultIsInf = resultIsNan.invert() & (
					(maxIsInfinity & minIsInfinity & effsub.invert()) |
					(maxIsInfinity ^ minIsInfinity) |
					finalExp.and_reduction()
				);

	auto constInfNanExp = Wrapper<WE-1, false>::generateSequence({1}).concatenate(resultIsNan | forbiddent_inf.invert());
	auto constInfNanSignif = Wrapper<WF, false>::generateSequence(resultIsNan | forbiddent_inf);

	auto constInfNan = constInfNanExp.concatenate(constInfNanSignif);
	auto finalRes = Wrapper<WE+WF, false>::mux(resultIsNan | resultIsInf, constInfNan, expSigRounded);

	auto bothZeros = maxIsZero & minIsZero;
	auto signBothZero = tmpMinSign & tmpMaxSign;

	Wrapper<1, false> isRoundDown{roundingMode == IEEERoundingMode::RoundDown};
	auto negZeroOp = resultIsZero & isRoundDown & (effsub | bothZeros.invert());
	auto signR = resultIsNan.invert() & // NaN forces sign to be zero
			((resultIsZero & (signBothZero | negZeroOp)) | // If summing two zeros of opposite sign, set result to zero
			(resultIsZero.invert() & maxSign));

#ifdef IEEE_ADDER_SORTEXP_DEBUG
	cerr << "===== IEEE ADD =====" << endl;
	cerr << "in0: " << to_string(in0) << endl;
	cerr << "in1: " << to_string(in1) << endl;
	cerr << "exp0: " << to_string(exp0) << endl;
	cerr << "exp1: " << to_string(exp1) << endl;
	cerr << "exp0IsZero: " << to_string(exp0IsZero) << endl;
	cerr << "exp1IsZero: " << to_string(exp1IsZero) << endl;
	cerr << "exp0IsAllOne: " << to_string(exp0IsAllOne) << endl;
	cerr << "exp1IsAllOne: " << to_string(exp1IsAllOne) << endl;
	cerr << "sign0: " << to_string(sign0) << endl;
	cerr << "sign1: " << to_string(sign1) << endl;
	cerr << "frac0: " << to_string(frac0) << endl;
	cerr << "frac1: " << to_string(frac1) << endl;
	cerr << "frac0IsZero: " << to_string(frac0IsZero) << endl;
	cerr << "frac1IsZero: " << to_string(frac1IsZero) << endl;
	cerr << "diff0: " << to_string(diff0) << endl;
	cerr << "diff1: " << to_string(diff1) << endl;
	cerr << "aligned: " << to_string(aligned) << endl;
	cerr << "swap: " << to_string(swap) << endl;
	cerr << "keep: " << to_string(keep) << endl;
	cerr << "maxExp: " << to_string(maxExp) << endl;
	cerr << "minExp: " << to_string(minExp) << endl;
	cerr << "expdiff: " << to_string(expdiff) << endl;
	cerr << "effsub: " << to_string(effsub) << endl;
	cerr << "tmpMinSign: " << to_string(tmpMinSign) << endl;
	cerr << "tmpMaxSign: " << to_string(tmpMaxSign) << endl;
	cerr << "maxSign: " << to_string(maxSign) << endl;
	cerr << "maxFrac: " << to_string(maxFrac) << endl;
	cerr << "minFrac: " << to_string(minFrac) << endl;
	cerr << "maxExpIsZero: " << to_string(maxExpIsZero) << endl;
	cerr << "minExpIsZero: " << to_string(minExpIsZero) << endl;
	cerr << "maxExpAllOne: " << to_string(maxExpAllOne) << endl;
	cerr << "minExpAllOne: " << to_string(minExpAllOne) << endl;
	cerr << "maxFracIsZero: " << to_string(maxFracIsZero) << endl;
	cerr << "minFracIsZero: " << to_string(minFracIsZero) << endl;
	cerr << "maxExpIsOne: " << to_string(maxExpIsOne) << endl;
	cerr << "maxIsInfinity: " << to_string(maxIsInfinity) << endl;
	cerr << "maxIsNaN: " << to_string(maxIsNaN) << endl;
	cerr << "maxIsZero: " << to_string(maxIsZero) << endl;
	cerr << "minIsInfinity: " << to_string(minIsInfinity) << endl;
	cerr << "minIsNaN: " << to_string(minIsNaN) << endl;
	cerr << "minIsZero: " << to_string(minIsZero) << endl;
	cerr << "bothSubNormals: " << to_string(bothSubNormals) << endl;
	cerr << "maxIsNormal: " << to_string(maxIsNormal) << endl;
	cerr << "minIsNormal: " << to_string(minIsNormal) << endl;
	cerr << "infinitySub: " << to_string(infinitySub) << endl;
	cerr << "resultIsNan: " << to_string(resultIsNan) << endl;
	cerr << "onlyOneSubnormal: " << to_string(onlyOneSubnormal) << endl;
	cerr << "explicitedMaxFrac: " << to_string(explicitedMaxFrac) << endl;
	cerr << "explicitedMinFrac: " << to_string(explicitedMinFrac) << endl;
	cerr << "maxShiftVal: " << to_string(maxShiftVal) << endl;
	cerr << "allShiftedOut: " << to_string(allShiftedOut) << endl;
	cerr << "shiftValue: " << to_string(shiftValue) << endl;
	cerr << "extendedMinFrac: " << to_string(extendedMinFrac) << endl;
	cerr << "shiftedMinFracSticky: " << to_string(shiftedMinFracSticky) << endl;
	cerr << "shiftedMinPBranch: " << to_string(shiftedMinPBranch) << endl;
	cerr << "adderShoulSub: " << to_string(adderShoulSub) << endl;
	cerr << "adderNegMask: " << to_string(adderNegMask) << endl;
	cerr << "adderMaxIn: " << to_string(adderMaxIn) << endl;
	cerr << "adderRes: " << to_string(adderRes) << endl;
	cerr << "adderResNeg: " << to_string(adderResNeg) << endl;
	cerr << "fracAdd: " << to_string(fracAdd) << endl;
	cerr << "shiftedMinSBranch: " << to_string(shiftedMinSBranch) << endl;
	cerr << "SubMaxIn: " << to_string(SubMaxIn) << endl;
	cerr << "subRes: " << to_string(subRes) << endl;
	cerr << "fracSub: " << to_string(fracSub) << endl;
	cerr << "resIsSub: " << to_string(resIsSub) << endl;
	cerr << "resIsAdd: " << to_string(resIsAdd) << endl;
	cerr << "fracRes: " << to_string(fracRes) << endl;
	cerr << "stickyRes: " << to_string(stickyRes) << endl;
	cerr << "z1: " << to_string(z1) << endl;
	cerr << "z0: " << to_string(z0) << endl;
	cerr << "lzc: " << to_string(lzc) << endl;
	cerr << "subnormalLimitVal: " << to_string(subnormalLimitVal) << endl;
	cerr << "lzcGreaterEqExp: " << to_string(lzcGreaterEqExp) << endl;
	cerr << "lzcSmallerEqExp: " << to_string(lzcSmallerEqExp) << endl;
	cerr << "lzcSmallerMaxVal: " << to_string(lzcSmallerMaxVal) << endl;
	cerr << "fullCancellation: " << to_string(fullCancellation) << endl;
	cerr << "normalOverflow: " << to_string(normalOverflow) << endl;
	cerr << "lzcOne: " << to_string(lzcOne) << endl;
	cerr << "subnormalOverflow: " << to_string(subnormalOverflow) << endl;
	cerr << "cancellation: " << to_string(cancellation) << endl;
	cerr << "overflow: " << to_string(overflow) << endl;
	cerr << "isLeftShiftLZC: " << to_string(isLeftShiftLZC) << endl;
	cerr << "isLeftShiftExp: " << to_string(isLeftShiftExp) << endl;
	cerr << "shiftFirstStage: " << to_string(shiftFirstStage) << endl;
	cerr << "normalisationShiftVal: " << to_string(normalisationShiftVal) << endl;
	cerr << "normalisedSignif: " << to_string(normalisedSignif) << endl;
	cerr << "significandPreRound: " << to_string(significandPreRound) << endl;
	cerr << "lsb: " << to_string(lsb) << endl;
	cerr << "roundBit: " << to_string(roundBit) << endl;
	cerr << "sticky: " << to_string(sticky) << endl;
	cerr << "deltaExpIsZero: " << to_string(deltaExpIsZero) << endl;
	cerr << "deltaExpIsMinusOne: " << to_string(deltaExpIsMinusOne) << endl;
	cerr << "deltaExpIsLZC: " << to_string(deltaExpIsLZC) << endl;
	cerr << "deltaExpExp: " << to_string(deltaExpExp) << endl;
	cerr << "deltaExpCin: " << to_string(deltaExpCin) << endl;
	cerr << "deltaBigPartIsZero: " << to_string(deltaBigPartIsZero) << endl;
	cerr << "deltaExpUnmasked: " << to_string(deltaExpUnmasked) << endl;
	cerr << "maskSequence: " << to_string(maskSequence) << endl;
	cerr << "deltaExpBeforeCorrection: " << to_string(deltaExpBeforeCorrection) << endl;
	cerr << "expPreRound: " << to_string(expPreRound) << endl;
	cerr << "expSigPreRound: " << to_string(expSigPreRound) << endl;
	cerr << "roundUpBit: " << to_string(roundUpBit) << endl;
	cerr << "unroundedInf: " << to_string(unroundedInf) << endl;
	cerr << "b0: " << to_string(b0) << endl;
	cerr << "b1: " << to_string(b1) << endl;
	cerr << "b2: " << to_string(b2) << endl;
	cerr << "forbiddent_inf: " << to_string(forbiddent_inf) << endl;
	cerr << "expSigRounded: " << to_string(expSigRounded) << endl;
	cerr << "finalExp: " << to_string(finalExp) << endl;
	cerr << "resultIsZero: " << to_string(resultIsZero) << endl;
	cerr << "resultIsInf: " << to_string(resultIsInf) << endl;
	cerr << "constInfNanExp: " << to_string(constInfNanExp) << endl;
	cerr << "constInfNanSignif: " << to_string(constInfNanSignif) << endl;
	cerr << "constInfNan: " << to_string(constInfNan) << endl;
	cerr << "finalRes: " << to_string(finalRes) << endl;
	cerr << "bothZeros: " << to_string(bothZeros) << endl;
	cerr << "signBothZero: " << to_string(signBothZero) << endl;
	cerr << "negZeroOp: " << to_string(negZeroOp) << endl;
	cerr << "signR: " << to_string(signR) << endl;
#endif

	return {signR.concatenate(finalRes)};
}
#endif // IEEE_ADDER_HPP

#endif // IEEE_ADDER_SORTEXP_HPP

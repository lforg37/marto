#ifndef IEEE_ADDER_HPP
#define IEEE_ADDER_HPP

#include "ieeefloats/ieeetype.hpp"
#include "ieeefloats/ieee_rounding.hpp"
#include "tools/static_math.hpp"
#include "primitives/lzoc_shifter.hpp"
#include "primitives/lzoc.hpp"
#include "primitives/shifter_sticky.hpp"

#ifdef IEEE_ADDER_DEBUG
#include <iostream>
#include "tools/printing.hpp"
using hint::to_string;
using std::cerr;
#endif

template<unsigned int WE, unsigned int WF, template<unsigned int, bool> class Wrapper>
inline IEEENumber<WE, WF, Wrapper> ieee_add_sub_impl(
	   IEEENumber<WE, WF, Wrapper> const & in0,
	   IEEENumber<WE, WF, Wrapper> const & in1,
	   IEEERoundingMode const roundingMode = IEEERoundingMode::RoundNearestTieEven
	)
{
	auto exp0 = in0.getExponent();
	auto exp1 = in1.getExponent();


	auto sign0 = in0.getSign();
	auto sign1 = in1.getSign();

	auto frac0 = in0.getFractionnalPart();
	auto frac1 = in1.getFractionnalPart();


	auto exp0IsZero = exp0.nor_reduction();
	auto exp1IsZero = exp1.nor_reduction();
	auto exp0AllOne = exp0.and_reduction();
	auto exp1AllOne = exp1.and_reduction();
	auto frac0IsZero = frac0.nor_reduction();
	auto frac1IsZero = frac1.nor_reduction();

	auto expfrac0 = exp0.concatenate(frac0);
	//cerr << "expfrac0 : " << expfrac0.unravel().get_str(2) << endl;
	auto expfrac1 = exp1.concatenate(frac1);
	//cerr << "expfrac1 : " << expfrac1.unravel().get_str(2) << endl;

	auto diff0 = exp0.modularSub(exp1);
	auto diff1 = exp1.modularSub(exp0);

	// Sorting according to exp
	auto swap = expfrac1 > expfrac0;
	//cerr << "swap : " << swap.unravel().get_str(2) << endl;

	auto maxExp = Wrapper<WE, false>::mux(swap, exp1, exp0);
	auto minExp = Wrapper<WE, false>::mux(swap, exp0, exp1);
	auto expdiff = Wrapper<WE, false>::mux(swap, diff1, diff0);

	//cerr << "Expdiff : " << expdiff.unravel().get_str(2) << endl;

	auto effsub = sign0 ^ sign1;

	auto maxSign = Wrapper<1, false>::mux(swap, sign1, sign0);
	auto minSign = Wrapper<1, false>::mux(swap, sign0, sign1);

	auto maxFrac = Wrapper<WF, false>::mux(swap, frac1, frac0);
	auto minFrac = Wrapper<WF, false>::mux(swap, frac0, frac1);

	// Special case detection
	auto maxExpIsZero = Wrapper<1, false>::mux(swap, exp1IsZero, exp0IsZero);
	auto minExpIsZero = Wrapper<1, false>::mux(swap, exp0IsZero, exp1IsZero);
	auto maxExpAllOne = Wrapper<1, false>::mux(swap, exp1AllOne, exp0AllOne);
	auto minExpAllOne = Wrapper<1, false>::mux(swap, exp0AllOne, exp1AllOne);
	auto maxFracIsZero = Wrapper<1, false>::mux(swap, frac1IsZero, frac0IsZero);
	auto minFracIsZero = Wrapper<1, false>::mux(swap, frac0IsZero, frac1IsZero);

	//Special case logic
	auto maxIsInfinity = maxExpAllOne & maxFracIsZero;
	auto maxIsNaN = maxExpAllOne & maxFracIsZero.invert();
	auto maxIsZero = maxExpIsZero & maxFracIsZero;
	auto minIsInfinity = minExpAllOne & minFracIsZero;
	auto minIsNaN = minExpAllOne & minFracIsZero.invert();
	auto minIsZero = minExpIsZero & minFracIsZero;

	auto bothSubNormals = maxExpIsZero;
	//cerr << "Both subnormals : " << bothSubNormals.unravel().get_str(2) << endl;
	auto maxIsNormal = maxExpIsZero.invert();
	auto minIsNormal = minExpIsZero.invert();

	auto infinitySub = maxIsInfinity & minIsInfinity & effsub;
	auto resultIsNan = maxIsNaN | minIsNaN |infinitySub;
	auto onlyOneSubnormal = minExpIsZero & maxExpIsZero.invert();
	//cerr << "Only One  subnormal :" << onlyOneSubnormal.unravel().get_str(2) << endl;

	// Reconstruct full fraction
	auto explicitedMaxFrac = maxIsNormal.concatenate(maxFrac);
	auto explicitedMinFrac = minIsNormal.concatenate(minFrac);

	//cerr << "explicit fracs : \nMax : " << explicitedMaxFrac.unravel().get_str(2) << endl <<
	//		"Min : " << explicitedMinFrac.unravel().get_str(2) << endl;

	//alignment
	auto maxShiftVal = Wrapper<WE, false>{WF+3};
	auto allShiftedOut = expdiff > maxShiftVal;
	//cerr << "All shifted out : " << allShiftedOut.unravel().get_str(2) << endl;

	auto shiftValue = Wrapper<WE, false>::mux(
					allShiftedOut,
					maxShiftVal,
					expdiff
				).modularSub(onlyOneSubnormal.template leftpad<WE>());
	//cerr << "Shift value : " << shiftValue.unravel().get_str(2) << endl;

	Wrapper<WF+3, false> extendedMinFrac = explicitedMinFrac.concatenate(Wrapper<2, false>{0});
	//cerr << "extminfrac : " << extendedMinFrac.unravel().get_str(2) << endl;

	auto shiftedMinFracSticky = shifter_sticky(extendedMinFrac, shiftValue);
	auto beforeComp = Wrapper<1, false>{0}.concatenate(shiftedMinFracSticky.template slice<WF + 3, 1>());
	//cerr << "beforeComp : " << beforeComp.unravel().get_str(2) << endl;
	auto shiftedMinFrac = beforeComp ^ Wrapper<WF+4, false>::generateSequence(effsub);
	//cerr << "shiftedMinFrac : " << shiftedMinFrac.unravel().get_str(2) << endl;

	auto stickyMinFrac = shiftedMinFracSticky.template get<0>();
	//cerr << "stickyMinFrac : " << stickyMinFrac.unravel().get_str(2) << endl;

	// Addition
	auto carryIn = effsub & stickyMinFrac.invert();
	//cerr << "carryIn : " << carryIn.unravel().get_str(2) << endl;
	auto extendedMaxFrac = explicitedMaxFrac.concatenate(Wrapper<1, false>{0}).concatenate(carryIn).template leftpad<WF+4>();
	//cerr << "extMaxfrac : " << extendedMaxFrac.unravel().get_str(2) << endl;

	auto signifcandResult = extendedMaxFrac + shiftedMinFrac;
	//cerr << "Signif result : " << signifcandResult.unravel().get_str(2) << endl;
	// Renormalization

	auto isNeg = signifcandResult.template get<WF + 4>();
	auto z1 = signifcandResult.template get<WF+3>();
	auto z0 = signifcandResult.template get<WF+2>();
	//cerr << "Z1 :" << z1.unravel().get_str(2) << endl
	//	<< "Z0  :" << z0.unravel().get_str(2) << endl;

	auto lzcInput = signifcandResult.template slice<WF+3, 1>();
	auto lzc = lzoc_wrapper(lzcInput, {0});

	//cerr << "LZC Input :" << lzcInput.unravel().get_str(2) << endl
	//	<< "LZC  :" << lzc.unravel().get_str(2) << endl;
	constexpr unsigned int lzcsize = hint::Static_Val<WF+3>::_storage;

	static_assert (lzcsize<=WE, "The adder works only for wE > log2(WF).\nAre you sure you need subnormals ?\nIf yes, contact us with your use case, we will be happy to make it work for you.");
	auto subnormalLimitVal = Wrapper<lzcsize, false>{WF+3};

	auto maxExpIsOne = (maxExp == Wrapper<WE, false>{1});

	auto lzcGreaterEqExp = (lzc.template leftpad<WE>() >= maxExp);
	auto lzcSmallerEqExp = (lzc.template leftpad<WE>() <= maxExp);
	auto lzcSmallerMaxVal = lzc < subnormalLimitVal;
	auto fullCancellation = lzcSmallerMaxVal.invert();
	//cerr << "Fulll cancel : " << fullCancellation.unravel().get_str(2) << endl;

	//cerr << "lzcsmaller" << lzcSmallerMaxVal.unravel().get_str(2) << endl;

	auto normalOverflow = z1;
	//cerr << "Normal overflow : " << normalOverflow.unravel().get_str(2) << endl;
	auto lzcOne = z1.invert() & z0;
	auto subnormalOverflow = lzcOne & maxExpIsZero;
	//cerr << "Subnormal overflow : " << subnormalOverflow.unravel().get_str(2) << endl;
	auto cancellation = z1.invert() & z0.invert();

	auto overflow = normalOverflow | subnormalOverflow;

	auto isLeftShiftLZC = overflow |
			(lzcSmallerMaxVal.invert() & bothSubNormals) |
			(cancellation & maxIsNormal & lzcSmallerEqExp) |
			(maxIsNormal & lzcSmallerMaxVal.invert());
	auto isLeftShiftExp = lzcSmallerMaxVal & lzcGreaterEqExp & maxIsNormal;

	auto shiftFirstStage = Wrapper<lzcsize, false>::mux(isLeftShiftLZC, lzc, Wrapper<lzcsize, false>{1});
	auto normalisationShiftVal = Wrapper<lzcsize, false>::mux(isLeftShiftExp, maxExp.template slice<lzcsize-1, 0>(), shiftFirstStage);

	//cerr << "Normalisation shift : " << normalisationShiftVal.unravel().get_str(2) << endl;
	auto normalisedSignif = signifcandResult << normalisationShiftVal;
	//cerr << "Normalised signif : " << normalisedSignif.unravel().get_str(2) << endl;
	auto significandPreRound = normalisedSignif.template slice<WF+2, 3>();
	//cerr << "Preround signif : " << significandPreRound.unravel().get_str(2) << endl;
	auto lsb = normalisedSignif.template get<3>();
	//cerr << "lsb : " << lsb.unravel().get_str(2) << endl;
	auto roundBit = normalisedSignif.template get<2>();
	//cerr << "roundBit : " << roundBit.unravel().get_str(2) << endl;
	auto sticky = stickyMinFrac | normalisedSignif.template slice <1, 0>().or_reduction();
	//cerr << "sticky : " << sticky.unravel().get_str(2) << endl;

	auto deltaExpIsZero = z1.invert() & (z0 ^ bothSubNormals);
	//cerr << "DeltaExpIsZero : " << deltaExpIsZero.unravel().get_str(2) << endl;
	auto deltaExpIsMinusOne = z1 | (z0 & bothSubNormals);
	//cerr << "deltaExpIsMinusOne : " << deltaExpIsMinusOne.unravel().get_str(2) << endl;
	auto deltaExpIsLZC = (z1 | z0 | bothSubNormals).invert() & lzcSmallerEqExp & lzcSmallerMaxVal;
	//cerr << "deltaExpIsLZC : " << deltaExpIsLZC.unravel().get_str(2) << endl;
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

	//cerr << "DeltaExp : " << deltaExpBeforeCorrection.unravel().get_str(2) << endl;
	//cerr << "DeltaEXP cin : " << deltaExpCin.unravel().get_str(2) <<endl;

	auto expPreRound = maxExp.addWithCarry(deltaExpBeforeCorrection, deltaExpCin).template slice<WE-1, 0>();
	//cerr << "expPreRound : " << expPreRound.unravel().get_str(2) << endl;
	auto expSigPreRound = expPreRound.concatenate(significandPreRound);
	//cerr << "expSigPreRound : " << to_string(expPreRound) << endl;

	auto roundUpBit = ieee_getRoundBit(maxSign, lsb, roundBit, sticky, roundingMode);
	//cerr << "rub : " << to_string(roundUpBit) << endl;

	auto unroundedInf = expPreRound.and_reduction();
	//cerr << "Infinity pre-round : " << to_string(unroundedInf) << endl;

	Wrapper<3, false> roundingCode{static_cast<uint8_t>(roundingMode)};
	auto b0 = roundingCode.template get<0>();
	auto b1 = roundingCode.template get<1>();
	auto b2 = roundingCode.template get<2>();
	auto forbiddent_inf = ((b2 & b0 & (b1 == maxSign)) | (roundingCode.or_reduction().invert())) & unroundedInf & maxIsInfinity.invert() & resultIsNan.invert();

	//auto roundUpBit = roundBit & (sticky | lsb);
	auto expSigRounded = expSigPreRound.modularAdd(roundUpBit.template leftpad<WE+WF>());
	auto finalExp = expSigRounded.template slice<WF+WE-1, WF>();

	auto resultIsZero = fullCancellation & finalExp.or_reduction().invert();
	//cerr << "ResIsZero : " << resultIsZero.unravel().get_str(2) << endl;
	auto resultIsInf = resultIsNan.invert() & (
					(maxIsInfinity & minIsInfinity & effsub.invert()) |
					(maxIsInfinity ^ minIsInfinity) |
					finalExp.and_reduction()
				);

	auto constInfNanExp = Wrapper<WE-1, false>::generateSequence({1}).concatenate(resultIsNan | forbiddent_inf.invert());
	auto constInfNanSignif = Wrapper<WF, false>::generateSequence(resultIsNan | forbiddent_inf);

	auto constInfNan = constInfNanExp.concatenate(constInfNanSignif);
	//cerr << "Constinfnan : " << to_string(constInfNan) << endl;
	auto finalRes = Wrapper<WE+WF, false>::mux(resultIsNan | resultIsInf, constInfNan, expSigRounded);

	auto bothZeros = maxIsZero & minIsZero;
	auto signBothZero = minSign & maxSign;

	Wrapper<1, false> isRoundDown{roundingMode == IEEERoundingMode::RoundDown};
	auto negZeroOp = resultIsZero & isRoundDown & (effsub | bothZeros.invert());
	auto signR = resultIsNan.invert() & // NaN forces sign to be zero
			((resultIsZero & (signBothZero | negZeroOp)) | // If summing two zeros of opposite sign, set result to zero
			(resultIsZero.invert() & maxSign));

#ifdef IEEE_ADDER_DEBUG
	cerr << "====== IEEE_ADD (cmp frac) ======" << endl;
	cerr << "exp0: " << to_string(exp0) << endl;
	cerr << "exp1: " << to_string(exp1) << endl;
	cerr << "sign0: " << to_string(sign0) << endl;
	cerr << "sign1: " << to_string(sign1) << endl;
	cerr << "frac0: " << to_string(frac0) << endl;
	cerr << "frac1: " << to_string(frac1) << endl;
	cerr << "exp0IsZero: " << to_string(exp0IsZero) << endl;
	cerr << "exp1IsZero: " << to_string(exp1IsZero) << endl;
	cerr << "exp0AllOne: " << to_string(exp0AllOne) << endl;
	cerr << "exp1AllOne: " << to_string(exp1AllOne) << endl;
	cerr << "frac0IsZero: " << to_string(frac0IsZero) << endl;
	cerr << "frac1IsZero: " << to_string(frac1IsZero) << endl;
	cerr << "expfrac0: " << to_string(expfrac0) << endl;
	cerr << "expfrac1: " << to_string(expfrac1) << endl;
	cerr << "diff0: " << to_string(diff0) << endl;
	cerr << "diff1: " << to_string(diff1) << endl;
	cerr << "swap: " << to_string(swap) << endl;
	cerr << "maxExp: " << to_string(maxExp) << endl;
	cerr << "minExp: " << to_string(minExp) << endl;
	cerr << "expdiff: " << to_string(expdiff) << endl;
	cerr << "effsub: " << to_string(effsub) << endl;
	cerr << "maxSign: " << to_string(maxSign) << endl;
	cerr << "minSign: " << to_string(minSign) << endl;
	cerr << "maxFrac: " << to_string(maxFrac) << endl;
	cerr << "minFrac: " << to_string(minFrac) << endl;
	cerr << "maxExpIsZero: " << to_string(maxExpIsZero) << endl;
	cerr << "minExpIsZero: " << to_string(minExpIsZero) << endl;
	cerr << "maxExpAllOne: " << to_string(maxExpAllOne) << endl;
	cerr << "minExpAllOne: " << to_string(minExpAllOne) << endl;
	cerr << "maxFracIsZero: " << to_string(maxFracIsZero) << endl;
	cerr << "minFracIsZero: " << to_string(minFracIsZero) << endl;
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
	cerr << "shiftedMinFracSticky: " << to_string(shiftedMinFracSticky) << endl;
	cerr << "beforeComp: " << to_string(beforeComp) << endl;
	cerr << "shiftedMinFrac: " << to_string(shiftedMinFrac) << endl;
	cerr << "stickyMinFrac: " << to_string(stickyMinFrac) << endl;
	cerr << "carryIn: " << to_string(carryIn) << endl;
	cerr << "extendedMaxFrac: " << to_string(extendedMaxFrac) << endl;
	cerr << "signifcandResult: " << to_string(signifcandResult) << endl;
	cerr << "isNeg: " << to_string(isNeg) << endl;
	cerr << "z1: " << to_string(z1) << endl;
	cerr << "z0: " << to_string(z0) << endl;
	cerr << "lzcInput: " << to_string(lzcInput) << endl;
	cerr << "lzc: " << to_string(lzc) << endl;
	cerr << "subnormalLimitVal: " << to_string(subnormalLimitVal) << endl;
	cerr << "maxExpIsOne: " << to_string(maxExpIsOne) << endl;
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
	cerr << "roundUpBit: " << to_string(roundUpBit) << endl;
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
	cerr << "=================================" << endl;
#endif
	return {signR.concatenate(finalRes)};
}
#endif // IEEE_ADDER_HPP

#include <iostream>
#include "utils.hpp"
#include "QAM_ConstellationGenerator.hpp"
#include "ModulateDemodulate.hpp"

constexpr uint nSym = 131072; // 2^17
constexpr uint nBits = nSym * 4;

int main()
{
	//region QAM16 constellation
	auto qam16 = renormed<4>(generateQAMConstellation<4>(), 1. / 0.5270462766947298 / 6);//taskQAM16;
	//endregion

	//region Reading reference symbols
	iqPoints2Pol refSym;
	refSym = readIqFromFile("refSym.txt", 2, nSym);
	iqPoints refSymX(nSym), refSymY(nSym);
	for (uint idx = 0; idx < nSym; ++idx)
	{
		refSymX[idx] = refSym[0][idx];
		refSymY[idx] = refSym[1][idx];
	}
	//endregion

	//region Reading distorted symbols
	iqPoints2Pol distSym;
	distSym = readIqFromFile("distSym.txt", 2, nSym);
	iqPoints distSymX(nSym), distSymY(nSym);
	for (uint idx = 0; idx < nSym; ++idx)
	{
		distSymX[idx] = distSym[0][idx];
		distSymY[idx] = distSym[1][idx];
	}
	//endregion

	//region Reference bits
//    std::vector<std::vector<bool>> refBitsVec = readMulticolFile<bool>("srcPermBitDataMod.txt", 2, nBits);
	std::bitset<nBits> refBitsX, refBitsY;
	refBitsX = demodulate<nSym, 4>(refSymX, qam16);
	refBitsY = demodulate<nSym, 4>(refSymY, qam16);
	//endregion

	//region Check mod-demod
	iqPoints tmp = modulate<nBits, 4>(refBitsX, qam16);
	for (uint i = 0; i < nSym; ++i)
		if (std::norm(tmp[i] - refSymX[i]) > 1e-16)
			throw std::runtime_error("incorrect mod demod at " + std::to_string(i) + ": " +
									 std::to_string(std::log10(std::norm(tmp[i] - refSymX[i]))));
	std::cout << "Mod demod test passed" << std::endl;

	//endregion

	//region Pre-restoration metrics
	std::cout << "Pre-restoration mean NMSE: " << NMSE_metrics(refSym, distSym) << " dB" << std::endl;
	//endregion


	//region Post-restoration metrics
	int M = 5;
	iqPoints2Pol restoredSym = calcCoeffVector<nSym, 5>(refSym, distSym);
	std::cout << "M=5, Post-restoration mean NMSE: " << NMSE_metrics(refSym, restoredSym, (2 * M + 1) * (2 * M + 1) * 2)
			  << " dB" << std::endl;

	M = 10;
	restoredSym = calcCoeffVector<nSym, 10>(refSym, distSym);
	std::cout << "M=10, Post-restoration mean NMSE: "
			  << NMSE_metrics(refSym, restoredSym, (2 * M + 1) * (2 * M + 1) * 2) << " dB" << std::endl;

	M = 20;
	restoredSym = calcCoeffVector<nSym, 20>(refSym, distSym);
	std::cout << "M=10, Post-restoration mean NMSE: "
			  << NMSE_metrics(refSym, restoredSym, (2 * M + 1) * (2 * M + 1) * 2) << " dB" << std::endl;
	//endregion
//    for (uint i = 0; i < nSym; ++i)
//    {
//        std::cout << restoredSym[0][i] << std::endl;
//    }
//    std::cout << U << std::endl;

	return 0;
}

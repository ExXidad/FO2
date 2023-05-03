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
	bits refBitsX, refBitsY;
	refBitsX = demodulate<4>(refSymX, qam16);
	refBitsY = demodulate<4>(refSymY, qam16);
	//endregion

	//region Check mod-demod
	iqPoints tmp = modulate<4>(refBitsX, qam16);
	for (uint i = 0; i < nSym; ++i)
		if (std::norm(tmp[i] - refSymX[i]) > 1e-16)
			throw std::runtime_error("incorrect mod demod at " + std::to_string(i) + ": " +
									 std::to_string(std::log10(std::norm(tmp[i] - refSymX[i]))));
	std::cout << "Mod demod test passed" << std::endl;

	//endregion

	//region Pre-restoration metrics
	std::cout << "Pre-restoration mean NMSE: " << NMSE_metrics(refSym, distSym) << " dB" << std::endl;
	//endregion

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	//region Post-restoration metrics
	std::fstream resultFile;
	resultFile.open(PROJECT_DIR + "/result.txt", std::ios::out);
	resultFile << "M\tNMSEX\tNMSEY\tNMSEAV\tBERX\tBERY\tBERAV" << std::endl;

	int Mmax = 20;
	for (int M = 20; M <= Mmax; ++M)
	{
		iqPoints2Pol restoredSym = restorationModel(refSym, distSym, M);
		bits restoredBitsX = demodulate<4>(restoredSym[0], qam16);
		bits restoredBitsY = demodulate<4>(restoredSym[1], qam16);
		double BERX = estimateBER(refBitsX, restoredBitsX);
		double BERY = estimateBER(refBitsY, restoredBitsY);

		resultFile
				<< M << "\t"
				<< NMSE_metrics(refSym[0], restoredSym[0]) << "\t"
				<< NMSE_metrics(refSym[1], restoredSym[1]) << "\t"
				<< NMSE_metrics(refSym, restoredSym) << "\t"
				<< BERX << "\t"
				<< BERY << "\t"
				<< 0.5 * (BERX + BERY) << "\t"
				<< std::endl;
		std::cout << "Progress: " << static_cast<double>(M+1)/Mmax*100. << "%" << std::endl;
	}

	resultFile.close();

	//endregion

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count()
			  << " minutes"
			  << std::endl;

	return 0;
}

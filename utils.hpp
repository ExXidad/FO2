//
// Created by Ivan Kalesnikau on 18.04.2023.
//

#ifndef FO2_UTILS_HPP
#define FO2_UTILS_HPP

#include <iostream>
#include <complex>
#include <bitset>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <ios>
#include <fstream>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <chrono>

#define PROJECT_DIR std::string(PROJECT_DIR_RAW)

using complex = std::complex<double>;
using iqPoint = complex;
using iqPoints = std::vector<iqPoint>;
using iqPoints2Pol = std::vector<iqPoints>;
using inDataType = std::vector<std::tuple<iqPoint, iqPoint, bool>>;
using string = std::string;
using DoubleVector = std::vector<double>;
template<int size>
using MatrixSqD = Eigen::Matrix<double, size, size>;
template<int size>
using MatrixSqC = Eigen::Matrix<complex, size, size>;
using bits = std::vector<bool>;

using namespace std::complex_literals;

template<uint constellationOrder>
using constellation = std::unordered_map<std::bitset<constellationOrder>, iqPoint>;

double estimateBER(const bits &refBits, const bits &distBits)
{
	if (((refBits.size() - distBits.size()) % 2) != 0)
		throw std::runtime_error("Bad dimensions");
	const uint n = distBits.size();
	const uint referenceShift = (refBits.size() - distBits.size()) / 2;
	uint errSum = 0;
	for (uint idx = 0; idx < n; ++idx)
		errSum += (refBits[idx + referenceShift] == distBits[idx]) ? 0 : 1;
	return static_cast<double>(errSum) / static_cast<double>(n);

}

double NMSE_metrics(const iqPoints &referenceSignal, const iqPoints &distortedSignal, const uint nExclude = 0)
{
	if (((referenceSignal.size() - distortedSignal.size()) % 2) != 0)
		throw std::runtime_error("Bad signal dimensions");
	if (referenceSignal.size() < distortedSignal.size())
		throw std::runtime_error("Bad signal dimensions");
	if (referenceSignal.empty() || distortedSignal.empty())
		throw std::runtime_error("Vectors should be non-empty");
	const uint n = distortedSignal.size();
	const uint referenceShift = (referenceSignal.size() - distortedSignal.size()) / 2;

	double norm = 0;
	double diff = 0;
	for (uint idx = nExclude; idx < n - nExclude; ++idx)
	{
		diff += std::pow(std::abs(distortedSignal[idx] - referenceSignal[idx + referenceShift]), 2);
		norm += std::pow(std::abs(referenceSignal[idx]), 2);
	}
	return 10. * std::log10(diff / norm);
}

double
NMSE_metrics(const iqPoints2Pol &referenceSignal, const iqPoints2Pol &distortedSignal,
			 const uint nExclude = 0)
{
	const uint N = referenceSignal.size();
	DoubleVector results(N);
	for (uint idx = 0; idx < N; ++idx)
		results[idx] = NMSE_metrics(referenceSignal[idx], distortedSignal[idx], nExclude);
	double tmp = 0;
	for (uint idx = 0; idx < N; ++idx)
		tmp += std::pow(10., results[idx] * 0.1);
	return 10. * std::log10(tmp / N);
}

std::vector<iqPoints> readIqFromFile(const string fname, const uint iqInRow, const uint nRows)
{
	std::vector<iqPoints> data(iqInRow, iqPoints(nRows));
	std::fstream file;
	file.open(PROJECT_DIR + "/" + fname, std::ios::in);
	if (!file.is_open())
		throw std::runtime_error("Failed to open file " + fname);
	for (uint idx = 0; idx < nRows; ++idx)
	{
		for (uint j = 0; j < iqInRow; ++j)
		{
			double x, y;
			file >> x;
			file >> y;
			data[j][idx] = x + 1.i * y;
		}
	}
	file.close();
	return data;
}

template<typename T>
std::vector<std::vector<T>> readMulticolFile(const string fname, const uint nCols, const uint nRows)
{
	std::vector<std::vector<T>> data(nCols, std::vector<T>(nRows));
	std::fstream file;
	file.open(PROJECT_DIR + "/" + fname, std::ios::in);
	if (!file.is_open())
		throw std::runtime_error("Failed to open file " + fname);
	for (uint idx = 0; idx < nRows; ++idx)
	{
		for (uint j = 0; j < nCols; ++j)
		{
			T x;
			file >> x;
			data[j][idx] = x;
		}
	}
	file.close();
	return data;
}

template<std::size_t N>
std::bitset<N> reversed(const std::bitset<N> &b)
{
	std::bitset<N> result = b;
	for (uint i = 0; i < N / 2; ++i)
	{
		result[i] = b[N - 1 - i];
		result[N - 1 - i] = b[i];
	}
	return result;
}

iqPoint getFromIdx(const iqPoints &vec, const int idx)
{
	const uint max = vec.size();
	if (idx >= 0 && idx < max)
		return vec[idx];
//	else if (idx < 0)
//		return vec[0];
//	else if (idx >= max)
//		return vec[max - 1];
	else throw std::runtime_error("wtf just happened");
}

template<typename T>
Eigen::MatrixX<T> regularized(Eigen::MatrixX<T> matrix)
{
	int N = matrix.rows();
	double power = 0;
	for (uint i = 0; i < N; ++i)
		power += std::pow(std::abs(matrix(i, i)), 2);
	power /= N;
	power = std::sqrt(power);
	for (uint i = 0; i < N; ++i)
		matrix(i, i) += 1e-10 * power;
	return matrix;
}

iqPoints2Pol restorationModel(const iqPoints2Pol &refSym, const iqPoints2Pol &distSym, const int M)
{
	if (distSym[0].size() % distSym[0].size() != 0) throw std::runtime_error("Incorrect dimension");
	const int blockSize = 2 * M + 1;
	const int matrixSize = 2 * blockSize * blockSize;
	const int exclusionN = 2 * M;
	const int N = distSym[0].size() - 2 * exclusionN;
	Eigen::MatrixXcd Ux(N, matrixSize), Uy(N, matrixSize);
	Eigen::VectorXcd cx(matrixSize), cy(matrixSize);
	Eigen::VectorXcd dx(N), dy(N);
	iqPoints yx(N), yy(N);

	//region Create Ux, Uy
	for (int k = 0; k < N; ++k)
	{
		for (int m = -M; m <= M; ++m)
			for (int n = -M; n <= M; ++n)
			{
				Ux(k, M + n + (M + m) * blockSize) =
						distSym[0][k + exclusionN + m] *
						distSym[0][k + exclusionN + n] *
						std::conj(distSym[0][k + exclusionN + n + m]);
				Ux(k, blockSize * blockSize + M + n + (M + m) * blockSize) =
						distSym[0][k + exclusionN + m] *
						distSym[1][k + exclusionN + n] *
						std::conj(distSym[1][k + exclusionN + n + m]);

				Uy(k, M + n + (M + m) * blockSize) =
						distSym[1][k + exclusionN + m] *
						distSym[1][k + exclusionN + n] *
						std::conj(distSym[1][k + exclusionN + n + m]);
				Uy(k, blockSize * blockSize + M + n + (M + m) * blockSize) =
						distSym[1][k + exclusionN + m] *
						distSym[0][k + exclusionN + n] *
						std::conj(distSym[0][k + exclusionN + n + m]);
			}
		dx[k] = refSym[0][k + exclusionN] - distSym[0][k + exclusionN];
		dy[k] = refSym[1][k + exclusionN] - distSym[1][k + exclusionN];
	}
	//endregion

	//region Least squares
//	Eigen::MatrixXcd tmpcx =
//			regularized<complex>(Ux.adjoint() * Ux).inverse() * Ux.adjoint() * dx;
//	Eigen::MatrixXcd tmpcy =
//			regularized<complex>(Uy.adjoint() * Uy).inverse() * Uy.adjoint() * dy;
//	for (uint idx = 0; idx < matrixSize; ++idx)
//	{
//		cx(idx) = tmpcx(idx, 0);
//		cy(idx) = tmpcy(idx, 0);
//	}
	cx = regularized<complex>(Ux.adjoint() * Ux).inverse() * Ux.adjoint() * dx;
	cy = regularized<complex>(Uy.adjoint() * Uy).inverse() * Uy.adjoint() * dy;
//	Eigen::LeastSquaresConjugateGradient<Eigen::MatrixXcd> lscgx, lscgy;
//	lscgx.compute(Ux);
//	lscgy.compute(Uy);
//	cx = lscgx.solve(dx);
//	cy = lscgy.solve(dy);
	//endregion

	//region Model output
	Eigen::VectorXcd yxVec = Ux * cx, yyVec = Uy * cy;
	for (uint k = 0; k < N; ++k)
	{
		yx[k] = yxVec(k) + distSym[0][k + exclusionN];
		yy[k] = yyVec(k) + distSym[1][k + exclusionN];
	}
	//endregion

	iqPoints2Pol result{yx, yy};
	return result;
}

#endif //FO2_UTILS_HPP


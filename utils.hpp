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

using namespace std::complex_literals;

template<uint constellationOrder>
using constellation = std::unordered_map<std::bitset<constellationOrder>, iqPoint>;

template<uint n>
double estimateBER(const std::bitset<n> &b1, const std::bitset<n> &b2)
{
    uint errSum = 0;
    for (uint i = 0; i < n; ++i)
        errSum += b1[i] != b2[i];
    return static_cast<double>(errSum) / static_cast<double>(n);
}

double NMSE_metrics(const iqPoints &referenceSignal, const iqPoints &distortedSignal)
{
    if (referenceSignal.size() != distortedSignal.size())
        throw std::runtime_error("Bad signal dimensions");
    if (referenceSignal.empty() || distortedSignal.empty())
        throw std::runtime_error("Vectors should be non-empty");

    const uint n = referenceSignal.size();

    double norm = 0;
    double diff = 0;
    for (uint idx = 0; idx < n; ++idx)
    {
        diff += std::pow(std::abs(distortedSignal[idx] - referenceSignal[idx]), 2);
        norm += std::pow(std::abs(distortedSignal[idx]), 2);
    }
    return 10. * std::log10(diff / norm);
}

double
NMSE_metrics(const std::vector<iqPoints> &referenceSignal, const std::vector<iqPoints> &distortedSignal)
{
    const uint N = referenceSignal.size();
    DoubleVector results(N);
    for (uint idx = 0; idx < N; ++idx)
        results[idx] = NMSE_metrics(referenceSignal[idx], distortedSignal[idx]);
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
    else if (idx < 0)
        return vec[0];
    else if (idx >= max)
        return vec[max - 1];
    else throw std::runtime_error("wtf just happened");
}

template<typename T, uint N>
Eigen::Matrix<T, N, N> regularized(Eigen::Matrix<T, N, N> matrix)
{
    double tmp = 0;
    for (uint i = 0; i < N; ++i)
        tmp += std::abs(matrix(i, i));
    tmp /= N;
    tmp = std::sqrt(tmp);
    matrix += Eigen::Matrix<T, N, N>::Identity() * 1e-10 * tmp;
    return matrix;
}

template<int averagingBlockSize, int M>
iqPoints2Pol calcCoeffVector(const iqPoints2Pol &sig)
{
    if (sig[0].size() % averagingBlockSize != 0) throw std::runtime_error("Incorrect dimension");
    constexpr int blockSize = 2 * M + 1;
    const uint N = sig[0].size();
    using UType = Eigen::Matrix<complex, averagingBlockSize, 2 * blockSize * blockSize>;
    UType Ux, Uy;
    Eigen::Vector<complex, 2 * blockSize * blockSize> cx, cy;
    Eigen::Vector<complex, averagingBlockSize> dx, dy;
    iqPoints yx(sig[0].size()), yy(sig[0].size());
    for (uint averagingBlockIdx = 0; averagingBlockIdx < N / averagingBlockSize; ++averagingBlockIdx)
    {
        //region Create Ux, Uy
        for (int k = 0; k < averagingBlockSize; ++k)
        {
            for (int m = -M; m <= M; ++m)
                for (int n = -M; n <= M; ++n)
                {
                    Ux(k, M + n + (M + m) * blockSize) =
                            getFromIdx(sig[0], k + averagingBlockIdx * averagingBlockSize + m) *
                            getFromIdx(sig[0], k + averagingBlockIdx * averagingBlockSize + n) *
                            std::conj(getFromIdx(sig[0],
                                                 k + averagingBlockIdx * averagingBlockSize + n + m));
                    Ux(k, blockSize * blockSize + M + n + (M + m) * blockSize) =
                            getFromIdx(sig[0], k + averagingBlockIdx * averagingBlockSize + m) *
                            getFromIdx(sig[1], k + averagingBlockIdx * averagingBlockSize + n) *
                            std::conj(getFromIdx(sig[1],
                                                 k + averagingBlockIdx * averagingBlockSize + n + m));

                    Uy(k, M + n + (M + m) * blockSize) =
                            getFromIdx(sig[1], k + averagingBlockIdx * averagingBlockSize + m) *
                            getFromIdx(sig[1], k + averagingBlockIdx * averagingBlockSize + n) *
                            std::conj(getFromIdx(sig[1],
                                                 k + averagingBlockIdx * averagingBlockSize + n + m));
                    Uy(k, blockSize * blockSize + M + n + (M + m) * blockSize) =
                            getFromIdx(sig[1], k + averagingBlockIdx * averagingBlockSize + m) *
                            getFromIdx(sig[0], k + averagingBlockIdx * averagingBlockSize + n) *
                            std::conj(getFromIdx(sig[0],
                                                 k + averagingBlockIdx * averagingBlockSize + n + m));
                }
            dx[k] = sig[0][k + averagingBlockIdx * averagingBlockSize];
            dy[k] = sig[1][k + averagingBlockIdx * averagingBlockSize];
        }
        //endregion

        //region Least squares
        Eigen::Matrix<complex, 2 * blockSize * blockSize, 1> tmpcx =
                regularized<complex, 2 * blockSize * blockSize>(Ux.adjoint() * Ux).inverse() * Ux.adjoint() * dx;
        Eigen::Matrix<complex, 2 * blockSize * blockSize, 1> tmpcy =
                regularized<complex, 2 * blockSize * blockSize>(Uy.adjoint() * Uy).inverse() * Uy.adjoint() * dy;
        for (uint idx = 0; idx < 2 * blockSize * blockSize; ++idx)
        {
            cx(idx) = tmpcx(idx, 0);
            cy(idx) = tmpcy(idx, 0);
        }
        //endregion

        //region Model output
        Eigen::Vector<complex, averagingBlockSize> yxVec = Ux * cx, yyVec = Uy * cy;
        for (uint k = 0; k < averagingBlockSize; ++k)
        {
            yx[k + averagingBlockIdx * averagingBlockSize] = yxVec(k);
            yy[k + averagingBlockIdx * averagingBlockSize] = yyVec(k);
        }
        //endregion
    }
    iqPoints2Pol result{yx,yy};
    return result;
}

#endif //FO2_UTILS_HPP


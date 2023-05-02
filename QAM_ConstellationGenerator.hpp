//
// Created by Ivan Kalesnikau on 21.03.2023.
//

#ifndef FO1_QAM_CONSTELLATIONGENERATOR_HPP
#define FO1_QAM_CONSTELLATIONGENERATOR_HPP

#include "utils.hpp"
#include "GrayCodeGenerator.hpp"

template<uint constellationOrder>
constellation<constellationOrder> generateQAMConstellation(const bool normalize = true)
{
    using namespace std::complex_literals;
    uint numberOfIQPoints = std::pow(2, constellationOrder);
    uint sideSize = std::pow(2, constellationOrder / 2);

    // Расчет iq точек
    std::vector<complex> iqPoints(numberOfIQPoints);
    const double h = 1. / (sideSize - 1);
    for (uint idx = 0; idx < numberOfIQPoints; ++idx)
        iqPoints[idx] = -0.5 + 0.5i - 1.i * h * static_cast<double>(idx % sideSize) +
                        h * static_cast<uint>(idx / sideSize);

    for (uint row = 0; row < sideSize; row += 2)
        std::reverse(std::next(iqPoints.begin(), row * sideSize),
                     std::next(iqPoints.begin(), (1 + row) * sideSize));

    // Подсчет нормы
    double constellationNorm = 0.;
    for (uint idx = 0; idx < numberOfIQPoints; ++idx)
        constellationNorm += std::norm(iqPoints[idx]);
    constellationNorm = normalize ? std::sqrt(constellationNorm / numberOfIQPoints) : 1.;


    // Код Грея
    std::vector<std::bitset<constellationOrder>> grayCode = generateGrayCode<constellationOrder>();

    // Итог
    constellation<constellationOrder> resultingConstellation;
    resultingConstellation.reserve(numberOfIQPoints);
    for (uint idx = 0; idx < numberOfIQPoints; ++idx)
        resultingConstellation.insert({grayCode[idx], iqPoints[idx] / constellationNorm});

    return resultingConstellation;
}

constellation<4> qam16SA()
{
    using namespace std::complex_literals;
    constellation<4> qam16{
            {0b0000, 1. + 1.i},
            {0b0010, 3. + 1.i},
            {0b0001, 1. + 3.i},
            {0b0011, 3. + 3.i},
            {0b0100, 1. - 1.i},
            {0b0110, 3. - 1.i},
            {0b0101, 1. - 3.i},
            {0b0111, 3. - 3.i},
            {0b1000, -1. + 1.i},
            {0b1010, -3. + 1.i},
            {0b1001, -1. + 3.i},
            {0b1011, -3. + 3.i},
            {0b1100, -1. - 1.i},
            {0b1110, -3. - 1.i},
            {0b1101, -1. - 3.i},
            {0b1111, -3. - 3.i}
    };
    // Подсчет нормы
    double constellationNorm = 0.;
    for (auto &[code, iqPoint]: qam16)
        constellationNorm += std::norm(iqPoint);
    constellationNorm = std::sqrt(constellationNorm / 16);
    for (auto &[code, iqPoint]: qam16)
        iqPoint /= constellationNorm;
    return qam16;
}

constellation<4> qam16SAn()
{
    using namespace std::complex_literals;
    constellation<4> qam16{
            {0b1101, -3. + 3.i},
            {0b0101, -1. + 3.i},
            {0b0111, 1. + 3.i},
            {0b1111, 3. + 3.i},
            {0b1001, -3. + 1.i},
            {0b0001, -1. + 1.i},
            {0b0011, 1. + 1.i},
            {0b1011, 3. + 1.i},
            {0b1000, -3. - 1.i},
            {0b0000, -1. - 1.i},
            {0b0010, 1. - 1.i},
            {0b1010, 3. - 1.i},
            {0b1100, -3. - 3.i},
            {0b0100, -1. - 3.i},
            {0b0110, 1. - 3.i},
            {0b1110, 3. - 3.i}
    };
    // Подсчет нормы
    double constellationNorm = 0.;
    for (auto &[code, iqPoint]: qam16)
        constellationNorm += std::norm(iqPoint);
    constellationNorm = std::sqrt(constellationNorm / 16);
    for (auto &[code, iqPoint]: qam16)
        iqPoint /= constellationNorm;
    std::cout << constellationNorm << std::endl;

    return qam16;
}

constellation<4> taskQAM16{{0b1111,1.-1.i},
                           {0b1011,1.+1.i},
                           {0b1110,-1.-1.i},
                           {0b0101,3.-3.i},
                           {0b1010,-1.+1.i},
                           {0b1100,-3.-1.i},
                           {0b1000,-3.+1.i},
                           {0b0111,1.-3.i},
                           {0b1001,3.+1.i},
                           {0b0000,-3.+3.i},
                           {0b1101,3.-1.i},
                           {0b0100,-3.-3.i},
                           {0b0011,1.+3.i},
                           {0b0010,-1.+3.i},
                           {0b0001,3.+3.i},
                           {0b0110,-1.-3.i}};


template<uint constellationOrder>
constellation<constellationOrder> renormed(const constellation<constellationOrder> &c, const double norm)
{
    constellation<constellationOrder> r = c;
    for (auto &[code, iqPoint]:r)
        iqPoint /= norm;
    return r;
}

#endif //FO1_QAM_CONSTELLATIONGENERATOR_HPP

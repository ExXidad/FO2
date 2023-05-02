//
// Created by Ivan Kalesnikau on 21.03.2023.
//

#ifndef FO1_MODULATEDEMODULATE_HPP
#define FO1_MODULATEDEMODULATE_HPP

#include "utils.hpp"

template<uint dataSize, uint constellationOrder>
std::vector<iqPoint> modulate(const std::bitset<dataSize> &data, constellation<constellationOrder> &constellation_)
{
    static_assert(dataSize % constellationOrder == 0);
    uint numberOfIQPoints = dataSize / constellationOrder;
    iqPoints modulated(numberOfIQPoints);
    for (uint idx = 0; idx < numberOfIQPoints; ++idx)
    {
        std::bitset<constellationOrder> dataBits;
        for (uint j = 0; j < constellationOrder; ++j)
            dataBits[j] = data[idx * constellationOrder + j];
        modulated[idx] = constellation_[dataBits];
    }
    return modulated;
}

template<uint iqSize, uint constellationOrder>
std::bitset<iqSize * constellationOrder>
demodulate(const std::vector<iqPoint> &iqPoints, constellation<constellationOrder> &constellation_)
{
    assert(iqSize == iqPoints.size());
    constexpr uint numberOfDataPoints = iqSize * constellationOrder;
    const uint numberOfConstellationPoints = std::pow(2, constellationOrder);
    std::bitset<numberOfDataPoints> data;
    for (uint iqPointIdx = 0; iqPointIdx < iqSize; iqPointIdx++)
    {
        std::vector<double> distances(numberOfConstellationPoints);
        for (uint j = 0; j < numberOfConstellationPoints; ++j)
            distances[j] = std::norm(iqPoints[iqPointIdx] - (std::next(constellation_.begin(), j)->second));
        uint minElIdx = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));
        std::bitset<constellationOrder> demodulatedBits = std::next(constellation_.begin(), minElIdx)->first;

        for (uint k = 0; k < constellationOrder; ++k)
            data[constellationOrder * iqPointIdx + k] = demodulatedBits[k];

    }
    return data;
}

#endif //FO1_MODULATEDEMODULATE_HPP

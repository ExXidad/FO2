//
// Created by Ivan Kalesnikau on 21.03.2023.
//

#ifndef FO1_GRAYCODEGENERATOR_HPP
#define FO1_GRAYCODEGENERATOR_HPP

#include "utils.hpp"

template <uint constellationOrder>
std::vector<std::bitset<constellationOrder>> generateGrayCode()
{
    using bitset = std::bitset<constellationOrder>;
    std::vector<bitset> result;
    result.reserve(std::pow(2, constellationOrder));
    for (int i = 0; i < (1 << constellationOrder); i++)
        result.emplace_back((i ^ (i >> 1)));
    return result;
}


#endif //FO1_GRAYCODEGENERATOR_HPP

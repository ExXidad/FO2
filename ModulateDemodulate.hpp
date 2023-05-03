//
// Created by Ivan Kalesnikau on 21.03.2023.
//

#ifndef FO1_MODULATEDEMODULATE_HPP
#define FO1_MODULATEDEMODULATE_HPP

#include "utils.hpp"

template<uint constellationOrder>
std::vector<iqPoint> modulate(const bits &data, constellation<constellationOrder> &constellation_)
{
	if(data.size() % constellationOrder != 0)
		throw std::runtime_error("incorrect dimensions");
	uint numberOfIQPoints = data.size() / constellationOrder;
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

template<uint constellationOrder>
bits
demodulate(const std::vector<iqPoint> &iqPoints, constellation<constellationOrder> &constellation_)
{
	const uint numberOfDataPoints = iqPoints.size() * constellationOrder;
	const uint numberOfConstellationPoints = std::pow(2, constellationOrder);
	bits data(iqPoints.size()*constellationOrder);
	for (uint iqPointIdx = 0; iqPointIdx < iqPoints.size(); iqPointIdx++)
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

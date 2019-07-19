#ifndef _FILL_HOLE_H_
#define _FILL_HOLE_H_
#include "stdafx.h"
#include "graph_typedefs.h"

template <typename TIterator>
std::vector<std::pair<int, int>> fillHole(TIterator holeBegin, TIterator holeEnd, G& g, const std::vector<std::pair<int, int>>& prohibitedEdges)
{
	std::vector<std::pair<int, int>> result;
	std::map<size_t, std::set<int>> myCurves;
	for (TIterator it = holeBegin; it != holeEnd; ++it)
	{
		for (const auto& p : g[*it].clusterPoints)
			myCurves[*it].insert(p.curve);
	}

	typedef std::pair<std::pair<int, int>, size_t> UglyType;
	std::vector<UglyType> nSharedCurves;

	for (TIterator it = holeBegin; it != holeEnd; ++it)
	{
		for (TIterator kIt = holeBegin; kIt != it; ++kIt)
		{
			std::pair<int, int> p = { (int)*it,(int)*kIt }, p1 = {(int)*kIt, (int)*it};
			if ((std::find(prohibitedEdges.begin(), prohibitedEdges.end(), p) != prohibitedEdges.end())
				|| (std::find(prohibitedEdges.begin(), prohibitedEdges.end(), p1) != prohibitedEdges.end()))
				continue;

			std::vector<int> intersection;
			std::set_intersection(myCurves[*it].begin(), myCurves[*it].end(), myCurves[*kIt].begin(), myCurves[*kIt].end(), std::back_inserter(intersection));
			nSharedCurves.push_back({ { *it, *kIt }, intersection.size() });
		}
	}

	std::map<size_t, bool> connected;
	std::sort(nSharedCurves.begin(), nSharedCurves.end(), [](const UglyType& one, const UglyType& two) {return one.second > two.second; });
	int nEdgesAdded = 0;
	for (auto it : nSharedCurves)
	{
		if ((std::distance(holeBegin,holeEnd) == 4) && (connected[it.first.first] || connected[it.first.second]))
			continue;
		connected[it.first.first] = true;
		connected[it.first.second] = true;

		std::cout << "CONNECTING vertex " << it.first.first << " to " << it.first.second << "(" << it.second << " shared curves)" << std::endl;

		result.push_back(it.first);
		++nEdgesAdded;
		if (nEdgesAdded == 2)
			break;
	}
	return result;
}

#endif

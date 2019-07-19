#include "stdafx.h"
#include "ContractDeg2.h"
#include "ChainDecomposition.h"
std::vector<std::pair<size_t, size_t>> topoGraph(const G & g)
{
	std::map<edge_descriptor, size_t> myChain;
	auto chains = chainDecomposition(g, myChain);
	std::vector<std::pair<size_t, size_t>> result(chains.size());
	for (int i = 0; i < chains.size(); ++i)
		result[i] = { chains[i].front().m_source,chains[i].back().m_target };
	return result;
}

std::vector<std::pair<size_t, size_t>> topoGraphHighValenceSeparated(const G & g, std::vector<std::vector<edge_descriptor>>& chainsSeparated, bool onlyLoops)
{
	std::map<edge_descriptor, size_t> myChain;
	auto chains = chainDecomposition(g, myChain);
	std::vector<std::pair<size_t, size_t>> result;
	for (int i = 0; i < chains.size(); ++i)
	{
		if (((boost::degree(chains[i].front().m_source, g) > 2 && boost::degree(chains[i].back().m_target, g) > 2) && !onlyLoops) ||
			(chains[i].front().m_source==chains[i].back().m_target))
		{
			if (chains[i].size() == 1)
			{
				//two high valence vertices connected by a single edge
				result.push_back({ chains[i].front().m_source,chains[i].back().m_target });
				chainsSeparated.push_back(chains[i]);
			}
			else
			{
				size_t sepVertex = chains[i][chains[i].size() / 2].m_source;
				if (g[sepVertex].split)
				{
					//find another one!
					double d = std::numeric_limits<double>::max();
					int bestSepVertex = -1;
					for (size_t vv = 1; vv < chains[i].size(); ++vv)
					{
						size_t vtx = chains[i][vv].m_source;
						if (!g[vtx].split)
						{
							double dist = fabs(chains[i].size() / 2 - vv);
							if (dist < d)
							{
								d = dist;
								bestSepVertex = vtx;
							}
						}
					}
					if (bestSepVertex >= 0)
						sepVertex = bestSepVertex;
				}
				result.push_back({ chains[i].front().m_source, sepVertex});
				result.push_back({ sepVertex, chains[i].back().m_target });
				
				for (int j = 0; j < chains[i].size(); ++j)
				{
					if (chains[i][j].m_source == sepVertex)
					{
						std::vector<edge_descriptor> chain1, chain2;
						chain1.insert(chain1.begin(), chains[i].begin(), chains[i].begin() + j);
						chain2.insert(chain2.begin(), chains[i].begin() + j, chains[i].end());
						chainsSeparated.push_back(chain1);
						chainsSeparated.push_back(chain2);
						break;
					}
				}
			}

		}
		else
		{
			result.push_back({ chains[i].front().m_source,chains[i].back().m_target });
			chainsSeparated.push_back(chains[i]);
		}
	}
	return result;
}

void contractDeg2(G & g)
{
	auto edges = topoGraph(g);
	for (size_t v = 0; v < boost::num_vertices(g); ++v)
	{
		boost::clear_vertex(v, g);
	}
	for (auto e : edges)
		boost::add_edge(e.first, e.second, g);
}

#include "stdafx.h"
#include "ChainDecomposition.h"

std::vector<std::vector<edge_descriptor>> chainDecomposition(const G & g, std::map<edge_descriptor, size_t>& myChain)
{
	std::vector<std::vector<edge_descriptor>> chains;
	//1. Break into chains, record their adjacencies into another graph
	edge_iter eit, eend;
	for (std::tie(eit, eend) = boost::edges(g); eit != eend; ++eit)
	{
		//find a seed edge for a new chain
		if (myChain.find(*eit) != myChain.end())
			continue;

		//grow the chain into both direction until it hits a high-valence vertex or reaches an end
		std::vector<edge_descriptor> newChain = { *eit };
		myChain[*eit] = chains.size();


		for (int dir : { -1, 1 })
		{
			size_t curVtx = dir == -1 ? newChain.front().m_source : newChain.back().m_target;
			bool continuationFound = true;
			//try extending starting from this vertex
			while (continuationFound && (boost::degree(curVtx, g) == 2))
			{
				continuationFound = false;
				oedge_iter oeit, oeend;
				for (std::tie(oeit, oeend) = boost::out_edges(curVtx, g); oeit != oeend; ++oeit)
				{
					if (myChain.find(*oeit) == myChain.end())
					{
						curVtx = oeit->m_target;
						if (dir == -1)
							newChain.insert(newChain.begin(), boost::edge(oeit->m_target, oeit->m_source, g).first);
						else
							newChain.push_back(*oeit);

						myChain[*oeit] = chains.size();
						continuationFound = true;
						break;
					}
				}
			}
		}

		chains.push_back(newChain);
	}
	return chains;
}

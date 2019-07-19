#include "stdafx.h"
#include "RemoveShortBranches.h"
#include "boost/graph/depth_first_search.hpp"
#include "ChainDecomposition.h"
#include <boost/graph/connected_components.hpp>
#include "Params.h"

double clusterWidth(const Cluster& cl)
{
	if (cl.clusterPoints.empty())
		return 0;
	return (cl.clusterPoints.back().p - cl.clusterPoints.front().p).norm();
}

double chainLength(const G& g, const std::vector<edge_descriptor>& chain, bool includeLast, bool projectOntoRoot)
{
	double length = 0;
	int limit = includeLast ? chain.size() : (int)chain.size() - 1;
	for (int j = 0; j < limit; ++j)
	{
		auto v_source = g[chain[j].m_source].location;
		auto v_target = g[chain[j].m_target].location;
		if (projectOntoRoot)
		{
			if (g[chain[j].m_source].root.norm() < 100)
				length += fabs((v_source - v_target).dot(g[chain[j].m_source].root));
		}
		else
			length += (v_source - v_target).norm();
	}
	return length;
}

void removeBranchesFilter1(G& g, bool onlyAtIntersections, const std::map<edge_descriptor,size_t>& myChainsBeforeSharpCornersSplit)
{
	std::cout << "Removing short branches...";
	int n = boost::num_vertices(g);
	bool somethingDeleted = true;

	while (somethingDeleted)
	{
		somethingDeleted = false;
		//let's see what that does
		std::map<edge_descriptor, size_t> myChain;
		auto chains = chainDecomposition(g, myChain);
		std::vector<double> chainLengths(chains.size());
		for (int i = 0; i < chains.size(); ++i)
			chainLengths[i] = chainLength(g, chains[i], true, false);

		auto orient = [&chains, &g](int chainIdx, int root)
		{
			if (chains[chainIdx].front().m_source == root)
				return chains[chainIdx];
			else
			{
				std::vector<edge_descriptor> reverseChain(chains[chainIdx].size());
				for (int i = 0; i < chains[chainIdx].size(); ++i)
				{
					reverseChain[chains[chainIdx].size() - 1 - i] = boost::edge(chains[chainIdx][i].m_target, chains[chainIdx][i].m_source, g).first;
				}
				return reverseChain;
			}
		};


		std::set<size_t> branchesToDelete;
		for (size_t v = 0; v < boost::num_vertices(g); ++v)
		{
			if (boost::degree(v, g) > 2)
			{
				if (onlyAtIntersections && g[v].root.norm() < 100)
					continue;

				std::map<size_t, std::vector<edge_descriptor>> adjChains;
				std::set<size_t> fixedChains;
				std::map<size_t, std::set<int>> intersectionsToRemove; //if we remove some chains, invalidate their intersections

				//we have a junction. see which adjacent chains are actually branches
				oedge_iter eit, eend;
				for (std::tie(eit, eend) = boost::out_edges(v, g); eit != eend; ++eit)
				{
					int chainIdx = myChain[*eit];
					size_t valEnd1 = boost::degree(chains[chainIdx].front().m_source, g), valEnd2 = boost::degree(chains[chainIdx].back().m_target, g);
					auto orientedChain = orient(chainIdx, v);
					adjChains[chainIdx] = orientedChain;
					if (valEnd1 != 1 && valEnd2 != 1)
						fixedChains.insert(chainIdx);
				}

				while (fixedChains.size() < adjChains.size() - 1)
				{
					//choose the longest one and mark as fixed
					int bestBranch = -1;
					double maxLength = 0;
					for (const auto& br : adjChains)
					{
						if (fixedChains.find(br.first) == fixedChains.end())
						{
							double len = chainLengths[br.first];
							if (maxLength < len)
							{
								maxLength = len;
								bestBranch = br.first;
							}
						}
					}

					if (bestBranch != -1)
						fixedChains.insert(bestBranch);
				}

				auto distTest = [&fixedChains, &adjChains, &g, &myChainsBeforeSharpCornersSplit](const Eigen::Vector2d& myPt, double myWidth)
				{
					
					Eigen::Vector2d myLoc = myPt;//.cast<int>();

						double myR = std::max(myWidth / 2,0.5+1e-6);
						for (size_t fc : fixedChains)
						{
							for (const auto& e : adjChains[fc])
							{
								double r = std::max(g[e.m_source].width / 2, 0.5 + 1e-6);
								Eigen::Vector2d hisLoc = g[e.m_source].location;//.cast<int>();
								double dist = (hisLoc - myLoc).squaredNorm();
								if (dist < std::pow(r + myR, 2))
									return true;
							}
							double r = g[adjChains[fc].back().m_target].width / 2;
							size_t lastVertex = adjChains[fc].back().m_target;
							Eigen::Vector2d hisLoc = g[lastVertex].location;//.cast<int>();
							double dist = (hisLoc - myLoc).squaredNorm();
							if (dist < std::pow(r + myR, 2))
								return true;
						}
					return false;
				};


				bool shouldDelete = true;
				for (const auto& adjChain : adjChains)
				{
					if (fixedChains.find(adjChain.first) != fixedChains.end())
						continue;

					double coveredLength = 0;
					double totalLength = 0;
					for (const auto& e : adjChain.second)
					{
						double edgeLength = (g[e.m_source].location - g[e.m_target].location).norm();
						int nSamples = 10;//std::max((int)edgeLength,3);
						std::vector<std::pair<Eigen::Vector2d, double>> samples;
						for (int k = 0; k < nSamples; ++k)
						{
							double alpha = k / (nSamples - 1);
							Eigen::Vector2d pt = (1 - alpha)*g[e.m_source].location + alpha*g[e.m_target].location;
							double w = g[e.m_source].width*(1 - alpha) + g[e.m_target].width*alpha;
							samples.push_back(std::make_pair(pt, w));
						}
						
						int nCovered = 0;
						for (int k = 0; k < nSamples; ++k)
						{
							if (distTest(samples[k].first, samples[k].second))
								nCovered++;
						}
						coveredLength += (nCovered / nSamples)*edgeLength;
						totalLength += edgeLength;
					}

					/* the last condition is not necessary - I'm just lazy. fix the stupid bug*/
					if (((totalLength - coveredLength > 1) && (adjChain.second.size()>2) && (coveredLength / totalLength < PRUNE_SHORT_BRANCHES_RATIO)) || (totalLength > 10))
						shouldDelete = false;
				}

				if (shouldDelete)
				{
					//std::cout << "Deleting. " << std::endl;
					for (const auto& br : adjChains)
					{
						if (fixedChains.find(br.first) == fixedChains.end())
						{
							branchesToDelete.insert(br.first);
						}
					}
				}
				/*else
					std::cout << "Keeping. " << std::endl;*/
			}
		}

		for (size_t branch : branchesToDelete)
		{
			somethingDeleted = true;
			for (const auto& e : chains[branch])
			{
				boost::remove_edge(e, g);
			}
		}
	}
	std::cout << "done." << std::endl;
}
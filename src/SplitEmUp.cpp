#include "stdafx.h"
#include "SplitEmUp.h"
#include "FillHole.h"
#include "ChainDecomposition.h"
#include "ContractDeg2.h"
void splitEmUpCorrectly(G& g)
{
	std::cout << "Splitting stuff... ";
	//1. Break into chains, record their adjacencies into another graph
	//2. For each chain in that little graph, if its degree is 3, 'open up' the chain
	//3. If degree is 4, split the chain into two
	//4. PROFIT

	int n = boost::num_vertices(g);

	//auto chains = chainDecomposition(g, myChain);
	std::vector<std::vector<edge_descriptor>> chains;
	auto tGraph = topoGraphHighValenceSeparated(g, chains, true);
	std::map<edge_descriptor, size_t> myChain;
	for (int i = 0; i < chains.size(); ++i)
	{
		for (const auto& e : chains[i])
			myChain[e] = i;
	}

	std::cout << "chains done... ";
	//2. For each chain in that little graph, if its degree is 3, 'open up' the chain

	auto avgTangent = [&](int chainIdx, bool end)
	{
		Eigen::Vector2d avg(0, 0);
		int count = 0;
		for (int i = 0; i < std::min(chains[chainIdx].size(), (size_t)20); ++i)
		{
			auto e = end ? chains[chainIdx][chains[chainIdx].size() - i - 1] : chains[chainIdx][i];
			Eigen::Vector2d edgeVec = g[e.m_target].location - g[e.m_source].location;
			avg += edgeVec.normalized();
			++count;
		}

		if (end)
			avg = -avg;
		Eigen::Vector2d result = avg / count;
		return result;
	};

	std::vector<ChainSplitInfo> chainsToSplit;
	std::map<std::pair<size_t, int>, bool> shouldSplit; //for convenience, has the same information as chainsToSplit. {chainIdx, endIdx}
	std::map<size_t, bool> willBeSplit; //for convenience, per vertex
	for (int i = 0; i < chains.size(); ++i)
	{
		ChainSplitInfo c = { i,{false,false} };
		
		for (int end = 0; end < 2; ++end)
		{
			size_t vtx = (end == 0) ? chains[i].front().m_source : chains[i].back().m_target;
			if (boost::degree(vtx, g) == 3)
			{
				//we have a candidate
				//should have a different sign of dot product with the root than the other 2 branches
				edge_descriptor myEdge = (end == 0) ? chains[i].front() : boost::edge(chains[i].back().m_target, chains[i].back().m_source, g).first;
				Eigen::Vector2d myEdgeVec = avgTangent(i, (end == 1));
				double myDot = myEdgeVec.dot(g[vtx].root);

				oedge_iter oeit, oend;
				bool shouldSplit = true;
				for (std::tie(oeit, oend) = boost::out_edges(vtx, g); oeit != oend; ++oeit)
				{
					if (*oeit == myEdge)
						continue;

					bool neighEnd = (chains[myChain[*oeit]].back().m_target == vtx);
					Eigen::Vector2d edgeVec = avgTangent(myChain[*oeit], neighEnd);
					if (edgeVec.dot(g[vtx].root)*myDot > 0)
						shouldSplit = false;
				}

				c.splitEnd[end] = shouldSplit;
			}
		}

		//the last condition is only to fight a bug in the sheriff's chin. I don't have time to debug it :(
		if (c.splitEnd[0] && c.splitEnd[1] && chains[i].size()<100)
		{
			willBeSplit[chains[i].front().m_source] = true;
			willBeSplit[chains[i].back().m_target] = true;
			shouldSplit[{i, 0}] = c.splitEnd[0];
			shouldSplit[{i, 1}] = c.splitEnd[1];
			chainsToSplit.push_back(c);
		}
	}

	//3. Do the actual splitting
	std::vector<std::vector<size_t>> chainVertices(chains.size());
	for (int i = 0; i < chains.size(); ++i)
	{
		std::vector<size_t> myChainVerts;
		myChainVerts.reserve(chains[i].size() + 1);
		for (const auto& e : chains[i])
			myChainVerts.push_back(e.m_source);
		myChainVerts.push_back(chains[i].back().m_target);
		chainVertices[i] = myChainVerts;
	}

	auto duplicateVtx = [&](size_t v)
	{
		size_t newV = boost::add_vertex(g);
		g[v].split = true;
		g[newV] = g[v];
		g[newV].clusterIdx = boost::num_vertices(g)-1;
		g[newV].sharpCorner = false;
		return newV;
	};

	auto createEdge = [&](size_t u, size_t v)
	{
		auto e = boost::add_edge(u, v, g);
		g[e.first].edgeCurve = -1; //not used in final optimization
		g[e.first].weight = 1; //default weight, HOPEFULLY not used in the optimization
		//std::cout << "Creating edge between " << u << " and " << v << std::endl;
		return e.first;
	};
	
	std::map<size_t, std::map<size_t, size_t>> newAdjacencies; //[vertex][incomingChainIdx]
	std::map<std::pair<size_t,int>,std::vector<size_t>> newVertices; //per {chain,dir}
	for (const auto& c : chainsToSplit)
	{
		const auto& myChainVertices = chainVertices[c.chain];

		for (int dir = 0; dir < 2; ++dir)
		{
			size_t v = (dir == 0) ? myChainVertices.front() : myChainVertices.back();
			std::vector<size_t> myNewVerts = { v };
			//todo: fill in info for the new vertex
			if (c.splitEnd[dir])
			{
				size_t newVertex = duplicateVtx(v);
				myNewVerts.push_back(newVertex);
				newVertices[{c.chain, dir}] = myNewVerts;

				oedge_iter eit, eend;
				for (std::tie(eit, eend) = boost::out_edges(v, g); eit != eend; ++eit)
				{
					int otherChain = myChain[*eit];
					if (otherChain != c.chain)
					{
						newAdjacencies[v][otherChain] = myNewVerts[newAdjacencies[v].size() % myNewVerts.size()]; //if we're not splitting at this end, then all the adjacent chains should connect to that vertex
					}
				}
			}
		}
	}

	auto adjustYJunction = [&](size_t& v1, size_t& v2, size_t yJunctionVtx)
	{
		std::map<size_t, bool> traversed;
		traversed[yJunctionVtx] = true;
		size_t initialV1 = v1, initialV2 = v2;
		bool haveSharedCurves;
		bool fixed;

		do
		{
			fixed = false;
			haveSharedCurves = false;
			traversed[v1] = true;
			traversed[v2] = true;
			for (const auto& c1 : g[v1].clusterPoints)
			{
				for (const auto& c2 : g[v2].clusterPoints)
				{
					if (c1.curve == c2.curve)
					{
						haveSharedCurves = true;
						break;
					}
				}
			}

			if (haveSharedCurves)
			{
				oedge_iter eit, eend;
				if (boost::degree(v1, g) == 2)
				{
					for (std::tie(eit, eend) = boost::out_edges(v1, g); eit != eend; ++eit)
					{
						if (!traversed[eit->m_target] && (boost::degree(eit->m_target,g)==2))
						{
							v1 = eit->m_target;
							fixed = true;
						}
					}
				}
				if (boost::degree(v2, g) == 2)
				{
					for (std::tie(eit, eend) = boost::out_edges(v2, g); eit != eend; ++eit)
					{
						if (!traversed[eit->m_target] && (boost::degree(eit->m_target, g) == 2))
						{
							v2 = eit->m_target;
							fixed = true;
						}
					}
				}
			}

			if ((boost::degree(v1, g) != 2) && (boost::degree(v2, g) != 2))
				break;
		} while (haveSharedCurves && fixed);

		if (haveSharedCurves)
		{
			//failed, return everything
			std::cout << "FAILED to adjust" << std::endl;
			v1 = initialV1;
			v2 = initialV2;
		}
	};

	//now we know which vertex to connect to
	//for every chain, now connect beginning vertices to the ending vertices

	for (size_t i = 0; i < chains.size(); ++i)
	{
		const auto& myChainVertices = chainVertices[i];
		if (!willBeSplit[myChainVertices.front()] && !willBeSplit[myChainVertices.back()])
			continue;

		std::array<std::vector<size_t>,2> activeVerts;
		//std::vector<edge_descriptor> edgesToRemove;
		for (int dir = 0; dir < 2; ++dir)
		{
			size_t vtx = dir == 0 ? myChainVertices.front() : myChainVertices.back();
			size_t adjVtx = (dir == 0) ? myChainVertices[1] : myChainVertices[myChainVertices.size() - 2];

			boost::remove_edge(vtx, adjVtx, g);
			//edgesToRemove.push_back(boost::edge(vtx, adjVtx, g).first);

			if (shouldSplit[{ i, dir }])
				activeVerts[dir] = newVertices[{i, dir}];
			else
			{
				if (willBeSplit[vtx])
				{
					auto it = newAdjacencies.find(vtx);
					assert(it != newAdjacencies.end());
					if (it == newAdjacencies.end())
						std::cout << "ERROR 1" << std::endl;

					auto it2 = it->second.find(i);
					assert(it2 != it->second.end());

					if (it2 == it->second.end())
					{
						std::cout << "ERROR 2, size: " << it->second.size() << ", vtx = " << vtx << std::endl;
					}

					activeVerts[dir] = { it2->second };
				}
				else
					activeVerts[dir] = { vtx };
			}
		}

		std::vector<std::pair<int, int>> finalEdges;
		if ((activeVerts[0].size() == 2) && (activeVerts[1].size() == 2))
		{
			std::vector<size_t> pseudoHole;
			//look at all the adjacent clusters
			std::map<size_t, size_t> connectedTo; //for convenience to replace the vertices in the edges returned by fillHole to finalEdges
			for (int dir = 0; dir < 2; ++dir)
			{
				size_t vtx = dir == 0 ? myChainVertices.front() : myChainVertices.back();
				for (auto it : newAdjacencies[vtx])
				{
					int chainIdx = it.first;
					if (chainIdx != i)
					{
						size_t theirVertex = (chainVertices[chainIdx].back() == vtx) ? chainVertices[chainIdx][chainVertices[chainIdx].size() - 2] : chainVertices[chainIdx][1];
						pseudoHole.push_back(theirVertex);
						connectedTo[theirVertex] = it.second;
					}
				}
			}

			std::cout << "Hole: ";
			for (size_t v : pseudoHole)
				std::cout << v << " ";
			std::cout << std::endl;

			if (pseudoHole.size() != 4)
			{
				std::cout << "Skipping this chain: ";
				for (size_t v : myChainVertices)
					std::cout << v << " ";
				return;
				continue;
			}

			auto origHole = pseudoHole;

			//make sure vertices (0,1) and (2,3) don't share any curves in common
			adjustYJunction(pseudoHole[0], pseudoHole[1], myChainVertices.front());
			adjustYJunction(pseudoHole[2], pseudoHole[3], myChainVertices.back());
			std::vector<std::pair<int, int>> prohibitedEdges = { {pseudoHole[0], pseudoHole[1]},{ pseudoHole[2], pseudoHole[3] } };

			std::cout << "Adjusted to: ";
			for (size_t v : pseudoHole)
				std::cout << v << " ";
			std::cout << std::endl;

			for (int i = 0; i < 4; ++i)
				connectedTo[pseudoHole[i]] = connectedTo[origHole[i]];

			auto holeEdges = fillHole(pseudoHole.begin(), pseudoHole.end(), g, prohibitedEdges);



			//now we have connections,but those are connecting adjacent chains
			//convert the vertex indices into my chain indices
			for (auto it : holeEdges)
			{
				int a = connectedTo[it.first], b = connectedTo[it.second];
				if ((activeVerts[1][0] == a) || (activeVerts[1][1] == a))
					std::swap(a, b);
				finalEdges.push_back({ a , b });
			}
		}
		else
		{
			for (size_t v1 : activeVerts[0])
				for (size_t v2 : activeVerts[1])
					finalEdges.push_back({ (int)v1, (int)v2 });
		}

		assert(!finalEdges.empty());
		std::vector<std::vector<size_t>> midChain;
		midChain.resize(finalEdges.size());
		for (int j = 1; j + 1 < myChainVertices.size(); ++j)
			midChain[0].push_back(myChainVertices[j]);
		
		for (int k = 1; k < finalEdges.size(); ++k)
		{
			//for each new edge, add a midChain
			std::vector<size_t> newMidChain;
			for (size_t v : midChain[0])
			{
				size_t newV = duplicateVtx(v);
				newMidChain.push_back(newV);
			}
			for (int j = 1; j < newMidChain.size(); ++j)
				createEdge(newMidChain[j - 1], newMidChain[j]);

			midChain[k] = newMidChain;
		}
		for (int j=0; j<finalEdges.size(); ++j)
		{
			if (midChain[j].size() > 2)
			{
				createEdge(finalEdges[j].first, midChain[j].front());
				createEdge(finalEdges[j].second, midChain[j].back());
			}
			else
				createEdge(finalEdges[j].first, finalEdges[j].second);
		}
	}

	std::cout << "Processing deg 4 verts: ";
	//similar processing for valence for vertices
	//I could have not duplicated the code, but I don't want to mess with the logic above
	for (size_t v = 0; v < boost::num_vertices(g); ++v)
	{
		if (boost::degree(v, g) == 4)
		{
			std::vector<size_t> pseudoHole;
			oedge_iter eit, eend;
			for (std::tie(eit, eend) = boost::out_edges(v, g); eit != eend; ++eit)
				pseudoHole.push_back(eit->m_target);

			auto holeEdges = fillHole(pseudoHole.begin(), pseudoHole.end(), g, {});
			boost::clear_vertex(v, g);
			for (auto e : holeEdges)
			{
				createEdge(e.first, e.second);
			}

			std::cout << "Deg 4 vertex: " << v << ", filling the hole" << std::endl;
		}
	}

	std::cout << "done." << std::endl;
}



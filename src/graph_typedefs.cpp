#include "stdafx.h"
#include "graph_typedefs.h"
#include <queue>
#include "boost/graph/filtered_graph.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include "FillHole.h"

std::set<edge_descriptor> contract_edges(const std::vector<std::vector<boost::graph_traits<G>::edge_descriptor>>& edgeLoops, G& g)
{
	typedef boost::graph_traits<G>::edge_descriptor ed;
	typedef boost::graph_traits<G>::out_edge_iterator oei;
	std::map<ed, bool> willBeRemoved;
	std::set<size_t> anchors;
	//std::cout << "Contracting loop: ";
	for (int i = 0; i < edgeLoops.size(); ++i)
	{
		for (const auto& e : edgeLoops[i])
		{
			willBeRemoved[e] = true;
		}
	}

	for (size_t i = 0; i < boost::num_vertices(g); ++i)
	{
		if (g[i].split)
		{
			anchors.insert(i);
			continue;
		}
		
		oei eit, e_end;
		for (std::tie(eit, e_end) = boost::out_edges(i, g); eit != e_end; ++eit)
		{
			if (!willBeRemoved[*eit])
			{
				anchors.insert(i);
				break;
			}
		}
	}

	//split loops into connected components
	std::vector<std::set<int>> loopIdx(boost::num_vertices(g));
	for (int i = 0; i < edgeLoops.size(); ++i)
	{
		for (const auto& e : edgeLoops[i])
			loopIdx[e.m_source].insert(i);
	}

	std::vector<std::set<int>> loopIdxOptimized(boost::num_vertices(g));
	for (size_t i = 0; i < boost::num_vertices(g); ++i)
	{
		if (loopIdx[i].size() != 1)
			loopIdxOptimized[i] = loopIdx[i];
	}

	std::vector<bool> processed(edgeLoops.size());
	bool finished = false;
	std::vector<int> componentByLoop(edgeLoops.size());
	std::vector<std::set<int>> components;

	std::vector<std::set<int>> adjacentLoops(edgeLoops.size());
	for (int i = 0; i < edgeLoops.size(); ++i)
	{
		adjacentLoops[i].insert(i);
		for (const auto& e : edgeLoops[i])
		{
			for (int j : loopIdxOptimized[e.m_source])
				adjacentLoops[i].insert(j);
		}
	}

	for (int i = 0; i < edgeLoops.size(); ++i)
	{
		if (processed[i])
			continue;

		std::set<int> newComponent = { i };
		
		bool foundContinuation = true;
		while (foundContinuation)
		{
			foundContinuation = false;
			for (int loop : newComponent)
			{
				if (processed[loop])
					continue;

				foundContinuation = true;
				processed[loop] = true;
				componentByLoop[loop] = components.size();
				for (int j : adjacentLoops[loop])
				{
					bool notInSet = newComponent.find(j) == newComponent.end();
					if (!processed[j] && notInSet)
						newComponent.insert(j);
				}
				break;
			}
		}

		components.push_back(newComponent);
		/*std::cout << "COMPONENT " << components.size() - 1 << ": ";
		for (int jj : newComponent)
		{
			std::cout << jj << " ";
		}
		std::cout << std::endl;*/
	}

	std::vector<std::set<size_t>> myAdjacentVerts(components.size());
	std::map<edge_descriptor, bool> keepThisEdge;
	for (int i = 0; i < edgeLoops.size(); ++i)
	{
		for (const auto& e : edgeLoops[i])
		{
			if (anchors.find(e.m_source) != anchors.end())
				myAdjacentVerts[componentByLoop[i]].insert(e.m_source);
		}
	}

	for (int i = 0; i < components.size(); ++i)
	{
		if (myAdjacentVerts[i].empty())
			continue;

		Eigen::Vector2d avg(0, 0);
		Eigen::Vector2d avgRoot(0, 0);
		int cnt = 0;
		for (int j : components[i])
		{
			for (const auto& e : edgeLoops[j])
			{
				avg += g[e.m_source].location;
				avgRoot += g[e.m_source].root;
				++cnt;
			}

		}
		avg /= cnt;
		avgRoot /= cnt;


		G sg; //subgraph
		std::map<size_t, size_t> localToGlobal, globalToLocal;

		std::set<size_t> componentVertices;
		for (int j : components[i])
			for (const auto& e : edgeLoops[j])
			{
				componentVertices.insert(e.m_source);
			}

		size_t chosenVertex = *myAdjacentVerts[i].begin();

		double distToAvg = std::numeric_limits<double>::max();
		for (size_t v : componentVertices)
		{
			double d = (g[v].location - avg).norm();
			if (distToAvg > d)
			{
				distToAvg = d;
				chosenVertex = v;
			}
		}
		
		for (size_t v : componentVertices)
		{
			size_t local = boost::add_vertex(sg);
			sg[local].clusterIdx = local;
			localToGlobal[local] = v;
			globalToLocal[v] = local;
		}
	
		std::set<edge_descriptor> componentEdges;
		for (int j : components[i])
			for (const auto& e : edgeLoops[j])
			{
				auto ee = boost::add_edge(globalToLocal[e.m_source], globalToLocal[e.m_target], sg);
				sg[ee.first].weight = (g[e.m_source].location-g[e.m_target].location).norm();
			}

		std::vector<vertex_descriptor> pDij(num_vertices(sg));
		std::vector<double> dDij(num_vertices(sg));
		auto predMap = make_iterator_property_map(pDij.begin(), get(&Cluster::clusterIdx, sg));
		auto distMap = make_iterator_property_map(dDij.begin(), get(&Cluster::clusterIdx, sg));
		boost::dijkstra_shortest_paths(sg, globalToLocal[chosenVertex],
			predecessor_map(predMap).
			distance_map(distMap).weight_map(get(&Edge::weight, sg)));

		if (componentVertices.find(41) != componentVertices.end())
		{
			std::cout << "!!!!!*** ";
			for (size_t k : myAdjacentVerts[i])
				std::cout << k << " ";
			std::cout << std::endl << " $$$ ";
			for (size_t k : componentVertices)
				std::cout << k << " ";
			std::cout << std::endl;
		}
		

		for (size_t k : myAdjacentVerts[i])
		{
			if (k == chosenVertex)
				continue;
			
			size_t cur = globalToLocal[k]; //local
			while (cur != globalToLocal[chosenVertex])
			{
				size_t prev = predMap[cur]; 
				keepThisEdge[boost::edge(localToGlobal[predMap[cur]],localToGlobal[cur],g).first] = true;
				cur = prev;
			}
		}
	}

	std::map<edge_descriptor, bool> deleted; //I shouldn't need this, but somehow....
	std::set<edge_descriptor> result;
	for (int i = 0; i < edgeLoops.size(); ++i)
	{
		for (const auto& e : edgeLoops[i])
		{
			if (!deleted[e] && (boost::edge(e.m_source, e.m_target, g).second) && !keepThisEdge[e])
			{
				deleted[e] = true;
				result.insert(e);
				boost::remove_edge(e, g);
			}
		}

	}
	return result;
}

G listToVec(const listG& g)
{
	G result;
	listG::vertex_iterator vit, vend;
	std::map<listG::vertex_descriptor, size_t> index;
	for (std::tie(vit, vend) = boost::vertices(g); vit != vend; ++vit)
	{
		size_t v = boost::add_vertex(result);
		index[*vit] = v;
		result[v] = g[*vit];
		result[v].clusterIdx = v;
	}

	listG::edge_iterator eit, eend;
	for (std::tie(eit, eend) = boost::edges(g); eit != eend; ++eit)
	{
		auto e = boost::add_edge(index[eit->m_source], index[eit->m_target], result);
		result[e.first] = g[*eit];
	}
	return result;
}

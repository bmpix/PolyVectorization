#include "stdafx.h"
#include "ContractLoops.h"
#include "IsLoopContractible.h"
std::vector<edge_descriptor> contractLoops(G & g, const cv::Mat & origMask, const std::vector<MyPolyline>& polys)
{
	using namespace boost;
	int n = boost::num_vertices(g);
	std::vector<edge_descriptor> result; //removed edges

	std::vector<std::vector<cv::Point2f>> realIncontractibleLoops;

	std::cout << "Reeb graph: " << num_vertices(g) << " vertices, " << num_edges(g) << " edges." << std::endl;

	int c = 0;
	//do
	//{
		std::map<size_t, std::vector < graph_traits < G >::vertex_descriptor >> p;

		std::cout << "Computing min spaning trees...";
		size_t treeRoot = 0;
		//create a boolean map if an edge is in the sp tree

		std::map<graph_traits<G>::edge_descriptor, bool> isInSpTree;
		std::vector<bool> vertexCovered(num_vertices(g));
		while (treeRoot < num_vertices(g))
		{
			p[treeRoot] = std::vector< graph_traits < G >::vertex_descriptor>(num_vertices(g));
			prim_minimum_spanning_tree(g, &p[treeRoot][0], boost::root_vertex(treeRoot).weight_map(boost::get(&Edge::weight, g)));

			vertexCovered[treeRoot] = true;

			for (size_t i = 0; i < num_vertices(g); ++i)
			{
				auto ee = edge(i, p[treeRoot][i], g);
				if (ee.second)
				{
					isInSpTree[ee.first] = true;
					vertexCovered[i] = true;
					vertexCovered[p[treeRoot][i]] = true;
				}
			}

			for (; treeRoot < num_vertices(g); ++treeRoot)
				if (!vertexCovered[treeRoot])
					break;
		}
		std::cout << "done. " << std::endl;

		//now for every edge not in the tree, find the smallest loop containing it
		auto eii = edges(g);

		for (auto it = eii.first; it != eii.second; ++it)
		{
			if (!isInSpTree[*it])
				g[*it].weight = 1e10;
		}


		std::vector < std::vector<edge_descriptor>> loops, incontractibleLoops /*for debug only*/;
		std::vector<std::vector<cv::Point2f>> realLoops;
		std::vector<edge_descriptor> origEdge, origEdgeForIncontractibleLoops;


		std::cout << "Computing loops...";
		c = 0;
		int tmpLoopIdx = 0;
		for (auto it = eii.first; it != eii.second; ++it)
		{
			if (!isInSpTree[*it])
			{
				std::vector<vertex_descriptor> pDij(num_vertices(g));
				std::vector<double> dDij(num_vertices(g));
				auto predMap = make_iterator_property_map(pDij.begin(), get(&Cluster::clusterIdx, g));
				auto distMap = make_iterator_property_map(dDij.begin(), get(&Cluster::clusterIdx, g));
				size_t source = it->m_source;
				dijkstra_shortest_paths(g, it->m_source,
					predecessor_map(predMap).
					distance_map(distMap).weight_map(get(&Edge::weight,g)));

				//now record the path
				std::vector<edge_descriptor> loop;
				auto cur = it->m_target;
				while (cur != it->m_source)
				{
					edge_descriptor ed = edge(cur, pDij[cur], g).first;
					loop.push_back(ed);
					cur = pDij[cur];
				}
				loop.push_back(*it);
				realLoops.push_back({});
				
				/*if (loop.size() > 100)
				{
					std::cout << "BIG LOOP: (size = " << loop.size() << "), ";
					std::cout << "non-tree edge: " << it->m_source << "-" << it->m_target << std::endl;
					std::cout << "loop: ";
					for (auto tt : loop)
					{
						std::cout << tt.m_source << " ";
					}
				}*/

				if (isLoopContractible(loop, origMask, g, polys,realLoops.back()))
				{
					loops.push_back(loop);
					c += loop.size();
					origEdge.push_back(*it);
				}
				else
				{
					std::cout << "Incontractible loop: ";
					for (const auto& e : loop)
					{
						std::cout << e.m_source << " ";
					}
					std::cout << std::endl;
					incontractibleLoops.push_back(loop);
					origEdgeForIncontractibleLoops.push_back(*it);
					realIncontractibleLoops.push_back(realLoops.back());
				}
				++tmpLoopIdx;
			}
		}
		std::cout << "done, found " << c << " edges to remove" << std::endl;

		std::cout << "Contracting loops...";
		/*for (int i = 0; i < loops.size(); ++i)
		{
			std::cout << "Loop " << i << std::endl;
			for (const auto& e : loops[i])
			{
				std::cout << e.m_source << " ";
			}
			std::cout << loops[i].back().m_target << std::endl;
		}*/
		auto removedEdges = contract_edges(loops, g);
		if (!removedEdges.empty())
			result.insert(result.end(), removedEdges.begin(), removedEdges.end());
		std::cout << "done." << std::endl;
	//} while (c != 0);
	

	std::cout << "all done." << std::endl;
	return result;
}

#include "ContractLoops2.h"
#include <boost/graph/hawick_circuits.hpp>
#include "IsLoopContractible.h"
struct Visitor
{
	 Visitor(std::vector<std::vector<size_t>>& loops)
	 :loops(loops){
	 }

	 template <typename Graph>
	void cycle(const std::vector<size_t>& p, const Graph & g)
	{
		//20 can be changed to anything large, 
		if ((p.size() > 2) && (p.size() < 20))
		{
			loops.push_back(p);
		}
	}

	std::vector<std::vector<size_t>>& loops;
};

std::vector<edge_descriptor> contractLoops2(G& g, const cv::Mat & origMask, const std::vector<MyPolyline>& polys)
{
	using dirG = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS>;
	dirG gCopy(boost::num_vertices(g));
	edge_iter eit, eend;
	for (std::tie(eit, eend) = boost::edges(g); eit != eend; ++eit)
	{
		boost::add_edge(eit->m_source, eit->m_target, gCopy);
		boost::add_edge(eit->m_target, eit->m_source, gCopy);
	}
	std::vector<std::vector<size_t>> loops;
	Visitor visitor(loops);
	boost::hawick_circuits(gCopy, visitor);

	std::vector<std::vector<cv::Point2f>> realLoops;
	std::vector<std::vector<edge_descriptor>> contractibleLoops;
	std::cout << "TOTAL # loops: " << loops.size() << std::endl;
	for (int i = 0; i < loops.size(); ++i)
	{
		const auto& loop = loops[i];
		std::vector<edge_descriptor> edge_loop;
		for (int j = 0; j < loop.size(); ++j)
		{
			edge_loop.push_back(boost::edge(loop[j], loop[(j + 1)% loop.size()],g).first);
		}

		std::vector<cv::Point2f> realLoop;
		if (isLoopContractible(edge_loop, origMask, g, polys, realLoop))
		{
			contractibleLoops.push_back(edge_loop);

			std::cout << "FOUND A CONTRACTIBLE LOOP: ";
			for (size_t j : loop)
				std::cout << j << " ";
			std::cout << std::endl;
		}
	}

	auto removedEdges = contract_edges(contractibleLoops, g);
	std::vector<edge_descriptor> result;
	result.insert(result.end(), removedEdges.begin(), removedEdges.end());
	return result;
}

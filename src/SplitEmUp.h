#ifndef _SPLIT_EM_UP_H_
#define _SPLIT_EM_UP_H_

#include "graph_typedefs.h"

struct ChainSplitInfo
{
	int chain;
	std::array<bool,2> splitEnd;
};

void splitEmUpCorrectly(G & g); //returns one-degree verts that became part of chains

#endif

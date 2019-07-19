#pragma once
#include "graph_typedefs.h"
std::vector<std::pair<size_t, size_t>> topoGraph(const G& g);
std::vector<std::pair<size_t, size_t>> topoGraphHighValenceSeparated(const G & g, std::vector<std::vector<edge_descriptor>>& chainsSeparated, bool onlyLoops=false);
void contractDeg2(G& g);


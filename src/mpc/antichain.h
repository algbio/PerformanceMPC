#pragma once
#include <vector>
#include <mpc/graph.h>

typedef std::vector<int> antichain;

antichain maxantichain_from_minflow(Flowgraph<Edge::Minflow> &minflow);
bool is_antichain(antichain &antichain, Graph &g);

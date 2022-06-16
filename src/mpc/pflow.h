#pragma once
#include "graph.h"
#include "naive.h"

std::unique_ptr<Flowgraph<Edge::Minflow>> pflowk3(Graph &g);
std::unique_ptr<Flowgraph<Edge::Minflow>> pflowk2(Graph &g);

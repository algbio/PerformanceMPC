#pragma once
#include "graph.h"
#include <cassert>
#include <algorithm>
#include <array>
#include <memory>
#include <vector>
#include <functional>

typedef std::vector<std::vector<int>> path_cover;

void maxflow_solve_edmonds_karp(Flowgraph<Edge::Maxflow>& fg);
void minflow_maxflow_reduction(Flowgraph<Edge::Minflow>&, std::function<void(Flowgraph<Edge::Maxflow>&)> maxflow_solver);
bool is_valid_minflow(Flowgraph<Edge::Minflow>&);
void minflow_maxflow_reduction(Graph&);
bool is_valid_cover(std::vector<std::vector<int>>&, Graph&);
path_cover minflow_reduction_path_recover(Flowgraph<Edge::Minflow>&);
path_cover minflow_reduction_path_recover_fast(Flowgraph<Edge::Minflow>&);
path_cover minflow_reduction_path_recover_faster(Flowgraph<Edge::Minflow>&);
void naive_minflow_solve(Flowgraph<Edge::Minflow>& flowgraph);
void maxflow_solve_edmonds_karp_DMOD(Flowgraph<Edge::Maxflow> &fg);

std::unique_ptr<Flowgraph<Edge::Minflow>> naive_minflow_reduction(Graph &g, std::function<int(int)> node_weight=[](int i){return 1;});
std::unique_ptr<Flowgraph<Edge::Minflow>> greedy_minflow_reduction(Graph &g, std::function<int(int)> node_weight=[](int i){return 1;});
std::unique_ptr<Flowgraph<Edge::Minflow>> greedy_minflow_reduction_sparsified(Graph &g, std::function<int(int)> node_weight=[](int i){return 1;});

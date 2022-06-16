#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include <mpc/graph.h>
#include <mpc/pflow.h>
#include <mpc/reach.h>
#include <mpc/naive.h>

std::unique_ptr<Graph> transitive_closure(Graph &g) {
	auto reach_idx = graph_reachability(g);
	auto result = std::make_unique<Graph>(g.n);
	for(int i=1; i<=g.n; i++) {
		for(int j=1; j<=g.n; j++) {
			if(i != j && reach_idx->reaches(i, j))
				result->add_edge(i, j);
		}
	}
	return result;
}

std::unique_ptr<Graph> transitive_reduction(Graph &g) {
	auto reach_idx = graph_reachability(g);
	Eigen::MatrixX<bool> A = Eigen::MatrixX<bool>::Zero(g.n, g.n);
	for(int i=1; i<=g.n; i++) {
		for(auto u:g.edge_out[i])
			A(i-1, u-1) = 1;
	}
	Eigen::MatrixX<bool> B = Eigen::MatrixX<bool>::Zero(g.n, g.n);
	for(int i=1; i<=g.n; i++) {
		for(int j=1; j<=g.n; j++) {
			if(i != j && reach_idx->reaches(i,j))
				B(i-1, j-1) = 1;
		}
	}
	Eigen::MatrixX<bool> AB = A*B;
	auto result = std::make_unique<Graph>(g.n);
	for(int i=1; i<=g.n; i++) {
		for(auto u:g.edge_out[i])
			if(AB(i-1,u-1) == 0)
				result->add_edge(i,u);
	}
	return result;
}

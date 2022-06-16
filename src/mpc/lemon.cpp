#include "graph.h"
#include "lemon.h"
#include <algorithm>
#include <lemon/list_graph.h>
#include <lemon/concepts/digraph.h>
#include <lemon/edmonds_karp.h>
#include <lemon/preflow.h>
#include <lemon/cycle_canceling.h>
#include <lemon/capacity_scaling.h>
#include <lemon/network_simplex.h>
#include <lemon/cost_scaling.h>
#include <lemon/bfs.h>
#include <lemon/dfs.h>
#include <vector>
#include <utility>
#include <memory>
#include <iostream>
#include <queue>

using namespace lemon;

// Preflow
void lemon_preflow(Flowgraph<Edge::Maxflow> &input_graph) {
	ListDigraph g;
	std::vector<ListDigraph::Node> v;
	for(int i=1; i<=input_graph.n; i++) {
		v.push_back(g.addNode());
	}
	std::vector<std::pair<ListDigraph::Arc, Edge::Maxflow*>> v_e;

	for(int i=1; i<input_graph.n; i++) {
		for(auto &[u,e]:input_graph.edge_out[i]) {
			v_e.push_back({g.addArc(v[i-1], v[u-1]), e});
		}
	}

	auto cap = ListDigraph::ArcMap<int>(g);
	for(auto &[a,e]:v_e) {
		cap[a] = e->capacity;
	}

	Preflow<ListDigraph> ek(g, cap, v[input_graph.source-1], v[input_graph.sink-1]);
	ek.run();
	for(auto &[a,e]:v_e) {
		e->flow = ek.flow(a);
	}
}

struct lemon_minflow {
	std::unique_ptr<ListDigraph::ArcMap<int>> demand, cost;
	ListDigraph::Node source, sink;
	std::vector<std::pair<ListDigraph::Arc, Edge::Minflow*>> v_e;
	ListDigraph g;

	lemon_minflow(Flowgraph<Edge::Minflow> &input_graph, bool mincost) {
		std::vector<ListDigraph::Node> v;
		for(int i=1; i<=input_graph.n; i++) {
			v.push_back(g.addNode());
		}
		for(int i=1; i<input_graph.n; i++) {
			for(auto &[u,e]:input_graph.edge_out[i]) {
				v_e.push_back({g.addArc(v[i-1], v[u-1]), e});
			}
		}
		demand = std::make_unique<ListDigraph::ArcMap<int>>(g);
		cost = std::make_unique<ListDigraph::ArcMap<int>>(g);

		source = v[input_graph.source-1];
		sink = v[input_graph.sink-1];
		auto ea = g.addArc(v[input_graph.source-1], v[input_graph.sink-1]);
		(*cost)[ea] = 0;
		for(auto &[a,e]:v_e) {
			(*demand)[a] = e->demand;
			if(mincost) {
				if(g.source(a) == v[input_graph.source-1]) {
					(*cost)[a] = input_graph.n;
				} else {
					(*cost)[a] = 1;
				}
			} else {
				if(g.source(a) == v[input_graph.source-1]) {
					(*cost)[a] = 1;
				} else {
					(*cost)[a] = 0;
				}
			}
		}
	}
};

// NetworkSimplex
void lemon_ns(Flowgraph<Edge::Minflow> &input_graph) {
	auto mf = lemon_minflow(input_graph, false);
	NetworkSimplex<ListDigraph> ns(mf.g);
	ns.stSupply(mf.source, mf.sink, input_graph.n);
	ns.lowerMap(*mf.demand);
	ns.costMap(*mf.cost);
	ns.run();
	for(auto &[a,e]:mf.v_e) {
		e->flow = ns.flow(a);
	}
}

// NetworkSimplex minimum length
void lemon_ns_minlen(Flowgraph<Edge::Minflow> &input_graph) {
	auto mf = lemon_minflow(input_graph, true);
	NetworkSimplex<ListDigraph> ns(mf.g);
	ns.stSupply(mf.source, mf.sink, input_graph.n);
	ns.lowerMap(*mf.demand);
	ns.costMap(*mf.cost);
	ns.run();
	for(auto &[a,e]:mf.v_e) {
		e->flow = ns.flow(a);
	}
}

// CostScaling
void lemon_cs(Flowgraph<Edge::Minflow> &input_graph) {
	auto mf = lemon_minflow(input_graph, false);
	CostScaling<ListDigraph> cs(mf.g);
	cs.stSupply(mf.source, mf.sink, input_graph.n);
	cs.lowerMap(*mf.demand);
	cs.costMap(*mf.cost);
	cs.run();
	for(auto &[a,e]:mf.v_e) {
		e->flow = cs.flow(a);
	}
}

// CostScaling minimum length
void lemon_cs_minlen(Flowgraph<Edge::Minflow> &input_graph) {
	auto mf = lemon_minflow(input_graph, true);
	CostScaling<ListDigraph> cs(mf.g);
	cs.stSupply(mf.source, mf.sink, input_graph.n);
	cs.lowerMap(*mf.demand);
	cs.costMap(*mf.cost);
	cs.run();
	for(auto &[a,e]:mf.v_e) {
		e->flow = cs.flow(a);
	}
}

void lemon_caps(Flowgraph<Edge::Minflow> &input_graph) {
	auto mf = lemon_minflow(input_graph, false);
	CapacityScaling<ListDigraph> cs(mf.g);
	cs.stSupply(mf.source, mf.sink, input_graph.n);
	cs.lowerMap(*mf.demand);
	cs.costMap(*mf.cost);
	cs.run();
	for(auto &[a,e]:mf.v_e) {
		e->flow = cs.flow(a);
	}
}

void lemon_caps_minlen(Flowgraph<Edge::Minflow> &input_graph) {
	auto mf = lemon_minflow(input_graph, true);
	CapacityScaling<ListDigraph> cs(mf.g);
	cs.stSupply(mf.source, mf.sink, input_graph.n);
	cs.lowerMap(*mf.demand);
	cs.costMap(*mf.cost);
	cs.run();
	for(auto &[a,e]:mf.v_e) {
		e->flow = cs.flow(a);
	}
}

void lemon_cc(Flowgraph<Edge::Minflow> &input_graph) {
	auto mf = lemon_minflow(input_graph, false);
	CycleCanceling<ListDigraph> cc(mf.g);
	cc.stSupply(mf.source, mf.sink, input_graph.n);
	cc.lowerMap(*mf.demand);
	cc.costMap(*mf.cost);
	cc.run();
	for(auto &[a,e]:mf.v_e) {
		e->flow = cc.flow(a);
	}
}

void lemon_cc_minlen(Flowgraph<Edge::Minflow> &input_graph) {
	auto mf = lemon_minflow(input_graph, true);
	CycleCanceling<ListDigraph> cc(mf.g);
	cc.stSupply(mf.source, mf.sink, input_graph.n);
	cc.lowerMap(*mf.demand);
	cc.costMap(*mf.cost);
	cc.run();
	for(auto &[a,e]:mf.v_e) {
		e->flow = cc.flow(a);
	}
}

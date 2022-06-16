#include "graph.h"

void lemon_preflow(Flowgraph<Edge::Maxflow> &input_graph);
void lemon_ns(Flowgraph<Edge::Minflow> &input_graph);
void lemon_ns_minlen(Flowgraph<Edge::Minflow> &input_graph);
void lemon_cs(Flowgraph<Edge::Minflow> &input_graph);
void lemon_cs_minlen(Flowgraph<Edge::Minflow> &input_graph);
void lemon_cc(Flowgraph<Edge::Minflow> &input_graph);
void lemon_cc_minlen(Flowgraph<Edge::Minflow> &input_graph);
void lemon_caps(Flowgraph<Edge::Minflow> &input_graph);
void lemon_caps_minlen(Flowgraph<Edge::Minflow> &input_graph);

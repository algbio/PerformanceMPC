# Description

This repository contains **C++ implementations** of different fast algorithms for **Minimum Path Cover** (MPC).

> **Definition (MPC, MA and Dilworth's theorem)**
> 
> Given a Directed Acyclic Graph (DAG) $G = (V, E)$, a *minimum path cover* (MPC) $\mathcal{P}$ of $G$ is a minimum-sized set of paths of $G$ such that every vertex of $V$ is contained in at least one path of $\mathcal{P}$. A *maxium antichain* (MA) $A$ of $G$ is a maximum-sized set of pairwise non-reachable vertices of $G$. Dilworth's theorem states that the number of paths in an MPC equals the number of vertices in an MA, this number $k$ is known as the *width*.

<p align="center">
  <img src="https://github.com/algbio/PerformanceMPC/assets/2347684/ba659bce-03ed-43b6-a9b4-70eb8e87b570" alt="MPC and MA of DAG" width="50%"/>
</p>

## Algorithms

All algorithms implemented are based on a simple folklore reduction to **minimum flow**. The reduction builds a flow network $\mathcal{G}$ with a global source $s$ and a global sink $t$, and where every vertex $v$ of $V$ is **split** into an edge $(v^{in}, v^{out})$ with **demand** equal to $d(v^{in}, v^{out}) = 1$.

<p align="center">
  <img src="https://github.com/algbio/PerformanceMPC/assets/2347684/e2950619-a092-4242-8c80-af53b3434ef0" alt="Minimum Flow Reduction"width="50%"/>
</p>

A **decomposition** of a minimum flow $f^*$ of $\mathcal{G}$ corresponds to an MPC of $G$.

### Initial solution

The corresponding minimum flow problem is reduced to a **maximum flow** instance by also providing a feasible flow $f$ of $\mathcal{G}$, equivalently a path cover in $G$, not necessarily an MPC. We implemented the following initial solutions:

- Naive: A path cover of $|V|$ paths, each path covers exactly one vertex. `naive` in program options.
- Greedy: $O(\log{|V|})$-approximation based on greedy set cover [[1]](#references). `greedy` in program options.
- Greedy Sparsified: Same as Greedy, but transitive edges are removed while obtaining new greedy paths. `greedy_sparsified` in program options.

### Solvers

We divide the flow solvers included into four categories.

#### MaxFlow-based own implementations

- `naive_minflow_solve`: Implements a simple DFS-based Ford-Fulkerson approach [[2]](#references)
- `maxflow_solve_edmonds_karp`: Implements Edmonds-Karp algorithm [[3]](#references)
- `maxflow_solve_edmonds_karp_DMOD`: Implements Dinitz's algorithm [[4]](#references)

#### MaxFlow-based LEMON [[5]](#references)

- `lemon_preflow`: LEMON solver using Goldberg-Tarjan algorithm [[6]](#references)

#### MinCostFlow-based LEMON [[5]](#references)

These solvers ignore the initial solution.

- `lemon_ns[_minlen]`: LEMON solver using Network Simplex [[7]](#references)
- `lemon_caps[_minlen]`: LEMON solver using Capacity Scaling [[8]](#references)
- `lemon_cc[_minlen]`: LEMON solver using Cycle Canceling [[9]](#references)
- `lemon_cs[_minlen]`: LEMON solver using Cost Scaling [[10]](#references)

For these solvers, if used with suffix `_minlen` it outputs an MPC of minimum total length (sets cost one to edges).

#### Parameterized algorithms (own implementations)

These solvers ignore the initial solution.

- `pflowk3`: The first parameterized linear time solution running in time $O(k^3|V|+|E|)$ [[11]](#references).
- `pflowk2`: Later improvement over `pflowk3` running in time $O(k^2|V|+|E|)$ [[12]](#references).

### Decomposition

All implementations use the same fast decomposition algorithm to obtain the MPC $\mathcal{P}$ from the minimum flow $f^*$.

## Pre-processing

There are also two heuristics for pre-processing $G$ before computing the MPC.

- `sparsify_dfs`: Removes transitive edges based on a DFS-traversal of $G$
- `contract_graph`: Contracts root-to-leaf paths (in trees induced subgraphs) into vertices

## Compiling

```
mkdir build
cd build
cmake ..
make
```
Run cmake with `-DCMAKE_BUILD_TYPE=RELEASE` to disable debugs enable `O3` flag and `-DCMAKE_BUILD_TYPE=RELEASE2` to also enable `-march=native`.

## Running Experiments

`build/prog/exp` can be used to run experiments, for more information use `-h` argument.
`ulimit -s unlimited` should be used if stack size is limited.
```Usage: ./exp [OPTIONS]

Some options:
  -f TEXT REQUIRED            Either a path to a file, which is of format
                              	<# of nodes> <# of edges>
                              	a line for each edge a->b of the form <a b>
                              OR
                              random_dag (N and M must be provided)
                              OR
                              random_x_chain (N, M, K must be provided)

  -r TEXT:{naive,greedy,greedy_sparsified}
                              initial solution to use
  -s TEXT:{pflowk2,pflowk3,lemon_ns[_minlen],lemon_cs[_minlen],lemon_cc[_minlen],lemon_caps[_minlen],naive_minflow_solve,lemon_preflow,maxflow_solve_edmonds_karp,maxflow_solve_edmonds_karp_DMOD} REQUIRED
                              solver to use

  --contract_graph BOOLEAN [0] 
                              contract graph in out tree
  --sparsify_dfs BOOLEAN [0]  sparsify graph beforehand using the dfs sparsification

  --seed INT [1337]           seed to use when generating graphs
  --transitive_reduction BOOLEAN [0] 
                              Use transitive reduction of the graph instead
  --transitive_closure BOOLEAN [0] 
                              Use transitive closure of the graph instead
  -N INT                      N parameter for the graph generation (number of vertices)
  -M INT                      M parameter for the graph generation (number of edges)
  -K INT                      K parameter for the graph generation (initial width)
```


## Additional Features

Even though the main objective of the repo. is to be used to *conduct performance experiments* and *compare different algorithms* for MPC, at the moment it is also possible to:

> :warning: These features are not clearly exposed to a final user as they are still work in progress.

- Compute an MA (see [`./src/mpc/antichain.cpp`](src/mpc/antichain.cpp))
- Compute a *weigthed* MA (see [`./prog/wac/`](prog/wac))
- Compute a *minimum chain cover* (MCC) (see [`./src/mpc/cc.cpp`](src/mpc/cc.cpp))
- Compute MPC-based reachability index (see [`./src/mpc/reach.cpp`](src/mpc/reach.cpp))
- Compute transitive closure/reduction (see [`src/mpc/transitive.cpp`](src/mpc/transitive.cpp))
- MPC heuristic preprocessing (*transitive edge sparsification* and *graph contraction*, see [`./src/mpc/preprocessing.cpp`](src/mpc/preprocessing.cpp))
- DAG generators (see [`./src/mpc/graph.cpp`](src/mpc/graph.cpp))

 ## References

- [1] Mäkinen, V., Tomescu, A. I., Kuosmanen, A., Paavilainen, T., Gagie, T., & Chikhi, R. (2019). Sparse dynamic programming on DAGs with small width. ACM Transactions on Algorithms (TALG), 15(2), 1-21.
- [2] Ford, L. R., & Fulkerson, D. R. (1956). Maximal flow through a network. Canadian journal of Mathematics, 8, 399-404.
- [3] Edmonds, J., & Karp, R. M. (1972). Theoretical improvements in algorithmic efficiency for network flow problems. Journal of the ACM (JACM), 19(2), 248-264.
- [4] Dinitz, Y. (2006). Dinitz’algorithm: The original version and Even’s version. In Theoretical Computer Science: Essays in Memory of Shimon Even (pp. 218-240). Berlin, Heidelberg: Springer Berlin Heidelberg.
- [5] Dezső, B., Jüttner, A., & Kovács, P. (2011). LEMON–an open source C++ graph template library. Electronic notes in theoretical computer science, 264(5), 23-45.
- [6] Goldberg, A. V., & Tarjan, R. E. (1988). A new approach to the maximum-flow problem. Journal of the ACM (JACM), 35(4), 921-940.
- [7] Dantzig, G. (1963). Linear programming and extensions. Princeton university press.
- [8] Edmonds, J., & Karp, R. M. (1972). Theoretical improvements in algorithmic efficiency for network flow problems. Journal of the ACM (JACM), 19(2), 248-264.
- [9] Goldberg, A. V., & Tarjan, R. E. (1989). Finding minimum-cost circulations by canceling negative cycles. Journal of the ACM (JACM), 36(4), 873-886.
- [10] Goldberg, A. V., & Tarjan, R. E. (1990). Finding minimum-cost circulations by successive approximation. Mathematics of Operations Research, 15(3), 430-466.
- [11] Cáceres, M., Cairo, M., Mumey, B., Rizzi, R., & Tomescu, A. I. (2022). Sparsifying, shrinking and splicing for minimum path cover in parameterized linear time. In Proceedings of the 2022 Annual ACM-SIAM Symposium on Discrete Algorithms (SODA) (pp. 359-376). Society for Industrial and Applied Mathematics.
- [12] Caceres, M., Cairo, M., Mumey, B., Rizzi, R., & Tomescu, A. I. (2022). Minimum path cover in parameterized linear time. arXiv preprint arXiv:2211.09659.

# Contact
 Any error, improvement or suggestion please contact the authors or create an issue in the repo.

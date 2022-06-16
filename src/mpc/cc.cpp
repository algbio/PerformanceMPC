#include <algorithm>
#include <vector>
#include <stack>
#include <mpc/graph.h>
#include <mpc/naive.h>
#include <memory>

typedef std::vector<std::vector<int>> chain_cover;

// O(pathlen)
std::unique_ptr<chain_cover> naive_chaincover_from_pathcover(path_cover &pc, Graph &g) {
	std::vector<bool> visited(g.n+1);
	auto cc = std::make_unique<chain_cover>(pc.size());
	for(int i=0; i<pc.size(); i++) {
		auto &path = pc[i];
		for(auto &u:path) {
			if(!visited[u]) {
				(*cc)[i].push_back(u);
				visited[u] = true;
			}
		}
	}
	return cc;
}

struct mergeable_dict_trie {
	struct node {
		int subtrie_size=0;
		node *left=nullptr, *right=nullptr;
		node(int subtrie_size=0):subtrie_size(subtrie_size) {
		}
	};
	node *root;
	int height = 1;
	mergeable_dict_trie() {
		root = new node();
	}
	mergeable_dict_trie(node *rt, int height):root(rt),height(height){
	}
	mergeable_dict_trie(int k) {
		int tmp = k;
		while(tmp) {
			height++;
			tmp /= 2;
		}
		std::vector<node*> nodes((1<<height)+1);
		for(int i=0; i<k; i++) {
			nodes[(1<<(height-1))+i] = new node(1);
		}
		for(int i=(1<<(height-1))-1; i>=1; i--) {
			nodes[i] = new node();
			nodes[i]->left = nodes[i*2];
			nodes[i]->right = nodes[i*2+1];
			assert(i*2+1 < nodes.size());
			if(nodes[i]->left != nullptr)
				nodes[i]->subtrie_size += nodes[i]->left->subtrie_size;
			if(nodes[i]->right != nullptr)
				nodes[i]->subtrie_size += nodes[i]->right->subtrie_size;
		}
		root = nodes[1];
	}
	int some() {
		auto cur = root;
		int ret = 0;
		for(int i=0; i<height-1; i++) {
			assert(cur != nullptr);
			ret <<= 1;
			if(cur->left != nullptr) {
				cur = cur->left;
			} else {
				cur = cur->right;
				ret |= 1;
			}
		}
		return ret;
	}
	void merge(node *root_from) {
		merge_recursive(root_from, root);
	}
	// for internal use
	void merge_recursive(node *from, node *to) {
		assert(from != nullptr && to != nullptr);
		to->subtrie_size += from->subtrie_size;
		if(from->left != nullptr) {
			if(to->left == nullptr) {
				to->left = from->left;
			} else {
				merge_recursive(from->left, to->left);
			}
		}
		if(from->right != nullptr) {
			if(to->right == nullptr) {
				to->right = from->right;
			} else {
				merge_recursive(from->right, to->right);
			}
		}
		delete from;
	}
	node* size_split(int s) {
		assert(0 < s && s <= root->subtrie_size);
		auto cur_old = root;
		node *new_root = new node();
		auto cur_new = new_root;
		while(s > 0) {
			cur_new->subtrie_size = s;
			cur_old->subtrie_size -= s;
			assert(cur_old->subtrie_size >= 0);
			if(cur_old->left != nullptr && cur_old->left->subtrie_size <= s) {
				// remove left
				cur_new->left = cur_old->left;
				s -= cur_new->left->subtrie_size;
				cur_old->left = nullptr;
			}
			if(cur_old->right != nullptr && cur_old->right->subtrie_size <= s) {
				// remove right
				cur_new->right = cur_old->right;
				s -= cur_new->right->subtrie_size;
				cur_old->right = nullptr;
			}
			if(s == 0)
				break;
			if(cur_old->left != nullptr) {
				auto lol = new node();
				cur_new->left = lol;
				cur_new = lol;
				cur_old = cur_old->left;
			} else {
				auto lol = new node();
				cur_new->right = lol;
				cur_new = lol;
				assert(cur_old->right != nullptr);
				cur_old = cur_old->right;
			}
		}
		return new_root;
	}
	void rd(node* node) {
		if(node->left != nullptr)
			rd(node->left);
		if(node->right != nullptr)
			rd(node->right);
		delete node;
	}
	~mergeable_dict_trie() {
		if(root != nullptr)
			rd(root);
	}
};


// O((|V|+|E|) log k)
std::unique_ptr<chain_cover> minflow_reduction_cc_fast(Flowgraph<Edge::Minflow> &fg) {
	auto v_r = [](int v){return (v+1)/2;}; // fg -> original Graph
	auto v_in = [](int v){return v*2-1;};
	std::vector<bool> visited(fg.n+1);
	std::vector<int> topo;
	int k = 0;
	for(auto &[u, e]:fg.edge_out[fg.source]) {
		for(int i=0; i<e->flow; i++) {
			k++;
		}
	}
	auto cc = std::make_unique<chain_cover>(k);
	auto dfs = [&fg, &visited, &topo, &v_r](auto dfs, int s, bool skip) {
		if(visited[s] || s == fg.sink)
			return;
		visited[s] = 1;
		for(auto &[u, e]:fg.edge_out[s]) {
			dfs(dfs, u, !skip);
		}
		if(!skip)
			topo.push_back(v_r(s));
	};
	dfs(dfs, fg.source, true);
	std::reverse(topo.begin(), topo.end());

	auto src_trie = mergeable_dict_trie(k);
	std::vector<mergeable_dict_trie> mdt_v(fg.n/2+1);
	for(auto &u:mdt_v)
		u.height = src_trie.height;

	for(auto &v:topo) {
		for(auto &[u,e]:fg.edge_in[v_in(v)]) {
			if(e->flow == 0)
				continue;
			if(u == fg.source) {
				mdt_v[v].merge(src_trie.size_split(e->flow));
			} else {
				mdt_v[v].merge(mdt_v[v_r(u)].size_split(e->flow));
			}
		}
		(*cc)[mdt_v[v].some()].push_back(v);
	}
	return cc;
}
// O(||P||) pathlen time 
std::unique_ptr<chain_cover> minflow_reduction_cc_naive(Flowgraph<Edge::Minflow> &fg) {
	auto v_r = [](int v){return (v+1)/2;}; // fg -> original Graph
	auto v_in = [](int v){return v*2-1;};
	std::vector<std::vector<int>> chains(fg.n/2+1);
	std::vector<bool> visited(fg.n+1);
	std::vector<int> topo;
	int k = 0;
	auto cc = std::make_unique<chain_cover>();
	for(auto &[u, e]:fg.edge_out[fg.source]) {
		for(int i=0; i<e->flow; i++) {
			cc->push_back({});
			if(chains[v_r(u)].size() == 0)
				(*cc)[k].push_back(v_r(u));
			chains[v_r(u)].push_back(k);
			k++;
		}
	}
	auto dfs = [&fg, &visited, &topo, &v_r](auto dfs, int s, bool skip) {
		if(visited[s] || s == fg.sink)
			return;
		visited[s] = 1;
		for(auto &[u, e]:fg.edge_out[s]) {
			dfs(dfs, u, !skip);
		}
		if(!skip)
			topo.push_back(v_r(s));
	};
	dfs(dfs, fg.source, true);
	std::reverse(topo.begin(), topo.end());

	for(auto &v:topo) {
		for(auto &[u,e]:fg.edge_in[v_in(v)]) {
			if(u == fg.source)
				continue;
			for(int i=0; i<e->flow; i++) {
				if(chains[v].size() == 0)
					(*cc)[*chains[v_r(u)].rbegin()].push_back(v);
				chains[v].push_back(*chains[v_r(u)].rbegin());
				chains[v_r(u)].pop_back();
			}
		}
	}
	return cc;
}

bool valid_chaincover(chain_cover &cover, Graph &g) {
	std::vector<bool> visited(g.n+1);
	for(auto &cc:cover) {
		if(cc.size() == 0)
			return false;
		for(auto &u:cc) {
			assert(u >= 1);
			assert(u <= g.n);
			visited[u] = true;
		}
	}
	for(int i=1; i<=g.n; i++)
		if(!visited[i])
			return false;
	return true;
}

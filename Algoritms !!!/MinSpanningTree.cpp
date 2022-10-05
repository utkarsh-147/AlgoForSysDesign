#pragma once
#include "union_find.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>
template<typename D>
struct Edge {
int a, b;
D d;
};
template<typename D>
D minimum_spanning_tree(int N, std::vector<Edge<D>> const &edges,
std::vector<int> &ans)
/* Given a weighted undirected graph, constructs a minimum-weight spanning
* tree and returns its cost. Populates ans with the indices in edges
* of the edges in the minimum spanning tree.
* If there is no minimum spanning tree (the graph is disconnected),
* returns std::numeric_limits<D>::max().
*/
{
std::vector<int> idx(edges.size());
std::iota(idx.begin(), idx.end(), 0);
std::sort(idx.begin(), idx.end(),
[&] (int i, int j) {
return edges[i].d < edges[j].d;
});
D cost = 0;
DisjointSets uf(N);
for (int i : idx) {
Edge<D> const &e = edges[i];
if (!uf.same_set(e.a, e.b)) {
uf.unite(e.a, e.b);
cost += e.d;
ans.push_back(i);
}
}
if (uf.n_disjoint == 1)
return cost;
else {
ans.clear();
return std::numeric_limits<D>::max();
}
}
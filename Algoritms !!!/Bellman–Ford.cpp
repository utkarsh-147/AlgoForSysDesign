#pragma once
#include <limits>
#include <tuple>
#include <vector>
template<typename D>
std::vector<D> bellman_ford(
std::vector<std::vector<std::pair<int, D>>> const &edges,
int source,
std::vector<int> &prev,
bool &negative_cycle,
D *inf_ptr = nullptr)
/* Returns a vector of the shortest distances to all nodes from the source.
* Populates prev with the predecessors to each node (-1 when no predecessor).
3
* negative_cycle is set if there are negative cycles in the graph.
* Nodes numbered from 0 to N - 1 where N = edges.size().
* Distances are INF if there is no path, -INF if arbitrarily short paths exist.
* If inf_ptr is not null, puts INF into *inf_ptr.
*/
{
constexpr D INF = std::numeric_limits<D>::max() / 2;
if (inf_ptr != nullptr) *inf_ptr = INF;
int const N = edges.size();
prev.assign(N, -1);
std::vector<D> dist(N, INF);
dist[source] = 0;
// N - 1 phases for finding shortest paths,
// N phases for finding negative cycles.
negative_cycle = false;
for (int ph = 0; ph < 2 * N - 1; ++ph)
// Iterate over all edges
for (int u = 0; u < N; ++u)
if (dist[u] < INF) // prevent catching negative INF -> INF edges
for (std::pair<int, D> const &edge : edges[u]) {
int v; D w;
std::tie(v, w) = edge;
if (dist[v] > dist[u] + w) {
if (ph < N - 1) {
dist[v] = dist[u] + w;
prev[v] = u;
}
else {
negative_cycle = true;
dist[v] = -INF;
}
}
}
return dist;
}
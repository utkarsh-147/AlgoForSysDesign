#pragma once
#include <limits>
#include <tuple>
#include <vector>
template<typename D>
std::vector<std::vector<D>> floyd_warshall(
std::vector<std::vector<std::pair<int, D>>> const &edges,
bool &negative_cycle,
D *inf_ptr = nullptr)
/* Returns a matrix of the shortest distances between all nodes.
* negative_cycle is set if there are negative cycles in the graph.
* Nodes numbered from 0 to N - 1 where N = edges.size().
* Distances are INF if there is no path, -INF if arbitrarily short paths exist.
* If inf_ptr is not null, puts INF into *inf_ptr.
*/
{
constexpr D INF = std::numeric_limits<D>::max() / 2;
if (inf_ptr != nullptr) *inf_ptr = INF;
int const N = edges.size();
// Initialize distance matrix
std::vector<std::vector<D>> dist(N, std::vector<D>(N, INF));
for (int u = 0; u < N; ++u) {
dist[u][u] = 0;
for (std::pair<int, D> const &edge : edges[u]) {
int v; D w;
std::tie(v, w) = edge;
dist[u][v] = std::min(dist[u][v], w);
}
}
// Main loop
for (int k = 0; k < N; ++k)
for (int i = 0; i < N; ++i)
for (int j = 0; j < N; ++j)
if (dist[i][k] < INF && dist[k][j] < INF)
dist[i][j] = std::min(dist[i][j], dist[i][k] + dist[k][j]);
// Propagate negative cycles
negative_cycle = false;
for (int i = 0; i < N; ++i)
for (int j = 0; j < N; ++j)
for (int k = 0; k < N; ++k)
if (dist[i][k] < INF && dist[k][j] < INF && dist[k][k] < 0) {
negative_cycle = true;
dist[i][j] = -INF;
}
return dist;
}

#pragma once
#include <limits>
#include <queue>
#include <tuple>
#include <vector>
template<typename D>
std::vector<D> dijkstra(
std::vector<std::vector<std::pair<int, D>>> const &edges,
int source, int target,
std::vector<int> &prev,
D *inf_ptr = nullptr)
/* Returns a vector of the shortest distances to all nodes from the source,
* but stops whenever target is encountered (set target = -1 to compute
* all distances).
* Populates prev with the predecessors to each node (-1 when no predecessor).
* Nodes numbered from 0 to N - 1 where N = edges.size().
* Distances are INF if there is no path.
* If inf_ptr is not null, puts INF into *inf_ptr.
*/
{
constexpr D INF = std::numeric_limits<D>::max() / 2;
if (inf_ptr != nullptr) *inf_ptr = INF;
int const N = edges.size();
std::vector<bool> visited(N, false);
prev.assign(N, -1);
std::vector<D> dist(N, INF);
dist[source] = 0;
using Entry = std::pair<D, int>;
std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> Q;
Q.emplace(0, source);
while (!Q.empty()) {
D node_dist; int node;
std::tie(node_dist, node) = Q.top();
Q.pop();
if (visited[node]) continue;
visited[node] = true;
if (node == target) break;
for (std::pair<int, D> const &edge : edges[node]) {
int child; D edge_len;
std::tie(child, edge_len) = edge;
if (!visited[child] && node_dist + edge_len < dist[child]) {
dist[child] = node_dist + edge_len;
prev[child] = node;
Q.emplace(dist[child], child);
}
}
}
return dist;
}

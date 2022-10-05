#pragma once
#include <limits>
#include <queue>
#include <vector>
template<typename D>
struct Edge {
int from, to;
D cap;
D flow;
4
Edge(int from, int to, D cap)
: from(from), to(to), cap(cap), flow(0) {}
Edge() = default;
int other(int origin) const {
return origin == from ? to : from;
}
D max_increase(int origin) const {
return origin == from ? cap - flow : flow;
}
void increase(int origin, D inc) {
flow += (origin == from ? inc : -inc);
}
bool saturated(int origin) const {
return max_increase(origin) == 0;
}
};
template<typename D>
std::vector<std::vector<int>> adj_list(int N, std::vector<Edge<D>> const &edges)
{
std::vector<std::vector<int>> adj(N);
for (int i = 0; i < (int) edges.size(); ++i) {
Edge<D> const &e = edges[i];
adj[e.from].push_back(i);
adj[e.to].push_back(i);
}
return adj;
}
template<typename D>
D max_flow_body(int N, std::vector<Edge<D>> &edges,
std::vector<std::vector<int>> const &adj, int source, int sink)
{
for (Edge<D> &e : edges) e.flow = 0;
auto augment =
[&] () -> D {
std::vector<int> pred(N, -1);
std::queue<int> Q;
Q.push(source);
bool found = false;
while (!found && !Q.empty()) {
int node = Q.front();
Q.pop();
for (int i : adj[node]) {
Edge<D> const &e = edges[i];
int other = e.other(node);
if (pred[other] == -1 && !e.saturated(node)) {
Q.push(other);
pred[other] = i;
if (other == sink) {
found = true;
break;
}
}
}
}
if (found) {
D max_inc = std::numeric_limits<D>::max();
int node = sink;
while (node != source) {
Edge<D> const &e = edges[pred[node]];
max_inc = std::min(max_inc,
e.max_increase(e.other(node)));
node = e.other(node);
}
node = sink;
while (node != source) {
Edge<D> &e = edges[pred[node]];
e.increase(e.other(node), max_inc);
node = e.other(node);
}
return max_inc;
}
return 0;
};
D max_flow = 0;
D inc;
do {
inc = augment();
max_flow += inc;
} while (inc > 0);
return max_flow;
}
template<typename D>
D max_flow(int N, std::vector<Edge<D>> &edges, int source, int sink)
/* Given a directed graph with edge capacities, computes the maximum flow
* from source to sink.
* Returns the maximum flow. Mutates input: assigns flows to the edges.
*/
{
std::vector<std::vector<int>> adj = adj_list(N, edges);
return max_flow_body(N, edges, adj, source, sink);
}
5
template<typename D>
D min_cut(int N, std::vector<Edge<D>> &edges,
int source, int sink, std::vector<int> &out)
/* Given a directed weighted graph, constructs a minimum cut such that one
* side A contains source and the other side B contains sink.
* Returns the capacity of the cut (sum of weights of edges from A to B, but not
* from B to A), which is equal to the max flow from from source to sink.
* Populates out with the nodes in A.
* Mutates input: assigns flows to the edges.
*/
{
std::vector<std::vector<int>> adj = adj_list(N, edges);
D flow = max_flow_body(N, edges, adj, source, sink);
out.clear();
std::vector<bool> visited(N, false);
std::queue<int> Q;
Q.push(source);
while (!Q.empty()) {
int node = Q.front();
Q.pop();
if (visited[node]) continue;
visited[node] = true;
out.push_back(node);
for (int i : adj[node]) {
Edge<D> const &e = edges[i];
int other = e.other(node);
if (!visited[other] && !e.saturated(node))
Q.push(other);
}
}
return flow;
}
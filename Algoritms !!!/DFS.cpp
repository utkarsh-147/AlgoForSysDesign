#pragma once
#include <algorithm>
#include <functional>
#include <stack>
#include <vector>
using NodeFunc = std::function<void (int)>;
using EdgeFunc = std::function<void (int, int)>;
1
void node_nop(int) {}
void edge_nop(int, int) {}
/* Depth-first search in a directed graph.
* Conceptually, a node has one of three colours:
* · White: untouched
* · Gray: in the queue
* · Black: has been popped from the queue
* This class offers three callbacks:
* · explore(from, to): called when exploring the edge `from -> to`
* · discover(node): called when `node` becomes gray
* · finish(node): called when `node` becomes black
*/
class DFS {
public:
int const N;
std::vector<std::vector<int>> const &adj;
std::vector<bool> visited;
EdgeFunc explore;
NodeFunc discover, finish;
DFS(std::vector<std::vector<int>> const &adj)
: N(adj.size()), adj(adj), visited(N, false),
explore(edge_nop), discover(node_nop), finish(node_nop)
{}
void dfs_from(int node) {
if (visited[node]) return;
visited[node] = true;
discover(node);
for (int child : adj[node]) {
explore(node, child);
dfs_from(child);
}
finish(node);
}
void run() {
for (int node = 0; node < N; ++node)
dfs_from(node);
}
DFS &on_explore(EdgeFunc f) {
explore = f;
return *this;
}
DFS &on_discover(NodeFunc f) {
discover = f;
return *this;
}
DFS &on_finish(NodeFunc f) {
finish = f;
return *this;
}
};
std::vector<int> topological_sort(std::vector<std::vector<int>> const &adj)
/* Given a directed graph in the form of an adjacency list, return a topological
* ordering of the graph's DFS forest: an ordering of the nodes such that a
* parent always appears before all of its descendants. If the graph is acyclic,
* this is a topological ordering of the graph itself.
*/
{
std::vector<int> result;
result.reserve(adj.size());
DFS(adj).on_finish([&] (int node) { result.push_back(node); }).run();
std::reverse(result.begin(), result.end());
return result;
}
std::vector<std::vector<int>> scc_kosaraju(
std::vector<std::vector<int>> const &adj)
/* Given a directed graph in the form of an adjacency list, return a vector of
* its strongly connected components (where each component is a list of nodes).
*/
{
int const N = adj.size();
std::vector<int> order = topological_sort(adj);
std::vector<std::vector<int>> adj_T(N);
for (int i = 0; i < N; ++i)
for (int j : adj[i])
adj_T[j].push_back(i);
std::vector<std::vector<int>> comps(1);
DFS dfs(adj_T);
dfs.on_finish([&] (int node) { (comps.end() - 1)->push_back(node); });
for (int node : order) {
if (!(comps.end() - 1)->empty())
comps.emplace_back();
dfs.dfs_from(node);
}
if ((comps.end() - 1)->empty()) comps.erase(comps.end() - 1);
return comps;
}
bool is_cyclic(std::vector<std::vector<int>> const &adj)
/* Given a directed graph in the form of an adjacency list, determine whether it
* contains a cycle or not.
*/
{
// Whether a node is currently in the recursive call stack
std::vector<bool> in_stack(adj.size(), false);
// Use an exception to stop the DFS early
class CycleFound {};
try {
DFS(adj)
.on_discover([&] (int node) { in_stack[node] = true; })
.on_finish([&] (int node) { in_stack[node] = false; })
.on_explore([&] (int, int child) {
if (in_stack[child]) throw CycleFound();
})
.run();
} catch (CycleFound) {
return true;
}
return false;
}
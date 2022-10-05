#include <cmath>
#include <limits>
#include <vector>
// Usage example:
// vector<vector<double>> A = {
// {2, 1, -1},
// {-3, -1, 2},
// {-2, 1, 2}
// };
// vector<double> b = {8, -11, -3};
// vector<double> x;
// GaussianSolver<double> solver(A, b, x);
// solver.solve()
// Solution to Ax = b is now in x.
// A and b are modified into the reduced row echelon form of the system.
// Other variables: solver.consistent, solver.rank, solver.determinant 
template<typename T>
class GaussianSolver {
private:
// Add row i to row j
void add(int i, int j, T factor) {
b[j] += factor * b[i];
for (int k = 0; k < M; ++k)
A[j][k] += factor * A[i][k];
}
void swap_rows(int i, int j) {
std::swap(b[i], b[j]);
for (int k = 0; k < M; ++k)
std::swap(A[i][k], A[j][k]);
determinant = -determinant;
}
void scale(int i, T factor) {
b[i] *= factor;
for (int k = 0; k < M; ++k)
A[i][k] *= factor;
determinant /= factor;
}
public:
std::vector<std::vector<T>> &A;
std::vector<T> &b;
std::vector<T> &x;
int const N, M;
// sqrt(std::numeric_limits<T>::epsilon()) by default, epsilon is too small!
T eps;
int rank;
T determinant;
bool consistent;
std::vector<bool> uniquely_determined;
GaussianSolver(std::vector<std::vector<T>> &A,
std::vector<T> &b,
std::vector<T> &x,
T eps = sqrt(std::numeric_limits<T>::epsilon()))
: A(A), b(b), x(x), N(A.size()), M(A[0].size()), eps(eps),
rank(0), determinant(1), consistent(true)
{
uniquely_determined.assign(M, false);
}
void solve() {
// pr, pc: pivot row and column
for (int pr = 0, pc = 0; pr < N && pc < M; ++pr, ++pc) {
// Find pivot with largest absolute value
int best_r = -1;
T best = 0;
for (int r = pr; r < N; ++r)
if (std::abs(A[r][pc]) > best) {
best_r = r;
best = std::abs(A[r][pc]);
}
if (std::abs(best) <= eps) { // No pivot found
--pr; // only increase pc in the next iteration
continue;
}
// Rank = number of pivots
++rank;
// Move pivot to top and scale to 1
swap_rows(pr, best_r);
scale(pr, (T) 1 / A[pr][pc]);
A[pr][pc] = 1; // for numerical stability
// Eliminate entries below pivot
for (int r = pr + 1; r < N; ++r) {
add(pr, r, -A[r][pc]);
A[r][pc] = 0; // for numerical stability
}
}
// Eliminate entries above pivots
for (int pr = rank - 1; pr >= 0; --pr) {
// Find pivot
int pc = 0;
while (std::abs(A[pr][pc]) <= eps) ++pc;
for (int r = pr - 1; r >= 0; --r) {
add(pr, r, -A[r][pc]);
A[r][pc] = 0; // for numerical stability
}
}
// Check for inconsistency: an equation of the form 0 = 1
for (int r = N - 1; r >= rank; --r)
if (std::abs(b[r]) > eps) {
consistent = false;
return;
}
// Calculate a solution for x
// One solution is setting all non-pivot variables to 0
x.assign(M, 0);
for (int pr = 0; pr < rank; ++pr)
for (int pc = 0; pc < M; ++pc)
if (std::abs(A[pr][pc]) > eps) { // Pivot; A[pr][pc] == 1
x[pc] = b[pr];
break;
}
// Mark variables as uniquely determined or not
for (int pr = 0; pr < rank; ++pr) {
int nonzero_count = 0;
int pc = -1;
for (int c = 0; c < M; ++c)
if (std::abs(A[pr][c]) > eps) {
if (nonzero_count == 0) pc = c; // Pivot
++nonzero_count;
}
if (nonzero_count == 1)
uniquely_determined[pc] = true;
}
}
};
7.4 Simplex algorithm (linear programming)
// Simplex algorithm for linear programming.
// Written using the theory from
// Cormen, Leiserson, Rivest, Stein: Introduction to Algorithms
#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
// For debug()
#include <cstdio>
#include <iostream>
template<typename T>
class SimplexSolver {
public:
// Standard form: Maximize c x subject to A x <= b; x >= 0.
// Slack form: Basic variables xb and nonbasic variables xn. x = [xb, xn]
// Maximize z = nu + c xn subject to xb = b - A xn; x >= 0.
T const eps;
int n; // number of nonbasic variables
int m; // number of constraints
std::vector<std::vector<T>> A;
std::vector<T> b;
std::vector<T> c;
T nu;
// Holds INDICES of all variables (ranging from 0 to n + m - 1)
std::vector<int> nonbasic_vars, basic_vars;
bool feasible;
SimplexSolver(std::vector<std::vector<T>> A,
std::vector<T> b,
std::vector<T> c,
T eps = sqrt(std::numeric_limits<T>::epsilon()))
: eps(eps), n(c.size()), m(b.size()), A(A), b(b), c(c), nu(0),
nonbasic_vars(n), basic_vars(m), feasible(true)
{
assert((int) A.size() == m && (int) A[0].size() == n);
// Transform from standard form to slack form
// Initially: nonbasic_vars: 0 to n - 1, basic_vars: n to n + m - 1
std::iota(nonbasic_vars.begin(), nonbasic_vars.end(), 0);
std::iota(basic_vars.begin(), basic_vars.end(), n);
}
// xn_e: "entering" variable: nonbasic -> basic
// xb_l: "leaving" variable: basic -> nonbasic
void pivot(int e, int l) {
std::swap(nonbasic_vars[e], basic_vars[l]);
int const e_new = l, l_new = e; // Just to avoid confusion
std::vector<std::vector<T>> A_new(A);
std::vector<T> b_new(b);
std::vector<T> c_new(c);
T nu_new;
// New constraint for xn_e: replace
// xb_l = b_l - A_lj xn_j
// with
// (*) xn_e = b_l / A_le - xb_l / A_le - A_lj xn_j / A_le (j ≠ e)
b_new[e_new] = b[l] / A[l][e];
A_new[e_new][l_new] = (T) 1 / A[l][e];
for (int j = 0; j < n; ++j)
if (j != e)
A_new[e_new][j] = A[l][j] / A[l][e];
// Substitute (*) in the other constraint equations:
// In each xb_i = b_i - A_ij xn_j (i ≠ l), replace A_ie xn_e
// with A_ie (b_l / A_le - xb_l / A_le - A_lj xn_j / A_le (j ≠ e))
for (int i = 0; i < m; ++i)
if (i != l) {
b_new[i] -= A[i][e] / A[l][e] * b[l];
A_new[i][l_new] = -A[i][e] / A[l][e];
for (int j = 0; j < n; ++j)
if (j != e)
A_new[i][j] -= A[i][e] * A[l][j] / A[l][e];
}
// Substitute (*) in the objective function:
// In nu + c_j xn_j, replace c_e xn_e with
// c_e (b_l / A_le - xb_l / A_le - A_lj xn_j / A_le (j ≠ e))
nu_new = nu + c[e] * b[l] / A[l][e];
c_new[l_new] = -c[e] / A[l][e];
for (int j = 0; j < n; ++j)
if (j != e)
c_new[j] -= c[e] * A[l][j] / A[l][e];
A = A_new;
b = b_new;
c = c_new;
nu = nu_new;
}
T solve(std::vector<T> &sol) {
initial_solution();
if (feasible)
return solve_body(sol);
else
return std::numeric_limits<T>::lowest();
}
T solve_body(std::vector<T> &sol) {
T const INF = std::numeric_limits<T>::max();
while (true) {
// Find a nonbasic variable with a positive coefficient in c
int e = -1;
for (int j = 0; j < n; ++j)
if (c[j] > 0) {
e = j;
break;
}
if (e == -1) break; // c <= 0; optimal solution reached
// Find the basic variable xb_l which most severely limits how much
// we can increase xn_e
T max_increase = INF;
int l = -1;
for (int i = 0; i < m; ++i) {
T inc = A[i][e] > 0 ? b[i] / A[i][e] : INF;
if (inc < max_increase) {
max_increase = inc;
l = i;
}
}
if (l == -1)
return INF;
else
pivot(e, l);
}
// Construct solution in terms of the original variables x[0]..x[n - 1]
sol.assign(n, 0);
for (int i = 0; i < m; ++i)
if (basic_vars[i] < n)
sol[basic_vars[i]] = b[i];
// Return optimum of objective function
return nu;
}
void initial_solution() {
// Find index l of basic variable xb_l with minimum b_l
int l = -1;
T b_min = std::numeric_limits<T>::max();
for (int i = 0; i < m; ++i)
if (b[i] < b_min) {
b_min = b[i];
l = i;
}
if (b_min >= 0) return;
// Add an extra nonbasic variable x0
++n;
nonbasic_vars.push_back(n + m - 1);
// Add -x0 to the LHS of every constraint
for (std::vector<T> &row : A)
row.push_back(-1);
// Change the objective function to -x0
std::vector<T> original_c = c;
T original_nu = nu;
c.assign(n, 0);
c[n - 1] = -1;
// Perform a pivot x0 <-> xb_l.
// After this, the basic solution is feasible.
pivot(n - 1, l);
// Find an optimal solution to this auxiliary problem
std::vector<T> sol;
T ans = solve_body(sol); // ans = optimum of -x0
if (ans >= -eps) {
// Is x0 basic?
auto it = std::find(basic_vars.begin(), basic_vars.end(),
n + m - 1);
if (it != basic_vars.end())
{
// Pivot with an arbitrary non-basic variable
int e = it - basic_vars.begin();
pivot(e, 0);
}
// Find the index of x0, now a non-basic variable
int j = std::find(nonbasic_vars.begin(), nonbasic_vars.end(),
n + m - 1) - nonbasic_vars.begin();
// Erase x0 from the constraints and the list of variables
for (std::vector<T> &row : A) row.erase(row.begin() + j);
nonbasic_vars.erase(nonbasic_vars.begin() + j);
--n;
// Restore original objective function, substituting basic variables
// with their RHS
nu = original_nu;
c.assign(n, 0);
// Loop over all originally non-basic variables
for (int var = 0; var < n; ++var) {
int j = std::find(nonbasic_vars.begin(), nonbasic_vars.end(),
var) - nonbasic_vars.begin();
if (j != n)
c[j] += original_c[var];
else {
int i = std::find(basic_vars.begin(), basic_vars.end(),
var) - basic_vars.begin();
// Substitute xb_i = b_i - A_ij xn_j
nu += original_c[var] * b[i];
for (int j = 0; j < n; ++j)
c[j] -= original_c[var] * A[i][j];
}
}
}
else {
--n;
feasible = false;
}
}
void debug() const {
printf("Nonbasic vars: ");
for (int j : nonbasic_vars) printf("x[%d] ", j);
std::cout << std::endl;
printf("Basic vars: ");
for (int i : basic_vars) printf("x[%d] ", i);
std::cout << std::endl;
std::cout << "Optimize " << nu;
for (int j = 0; j < n; ++j) {
std::cout << " + (" << c[j] << ") * ";
printf("x[%d]", nonbasic_vars[j]);
}
puts("");
for (int i = 0; i < m; ++i) {
printf("x[%d] = ", basic_vars[i]);
std::cout << "(" << b[i] << ")";
for (int j = 0; j < n; ++j) {
std::cout << " - (" << A[i][j] << ") * ";
printf("x[%d]", nonbasic_vars[j]);
}
puts("");
}
puts("");
}
}

#include <limits>
template<typename P>
int int_binsearch_last(P p, int lo, int hi)
/* Takes a predicate p: int -> bool.
* Assumes that there is some I (lo <= I < hi) such that p(i) = (i <= I),
* i.e. p starts out true and then switches to false after I.
* Finds and returns I (the last number i such that p(i) is true),
* or lo - 1 if no such I exists (including if lo = hi).
*/
{
while (hi - lo > 1) {
int mid = lo + (hi - lo) / 2;
if (p(mid)) lo = mid; // mid <= I
else hi = mid; // I < mid
}
return (lo != hi && p(lo)) ? lo : lo - 1;
}
template<typename P>
int int_binsearch_first(P p, int lo, int hi)
/* Takes a predicate p: int -> bool.
* Assumes that there is some I (lo <= I < hi) such that p(i) = (i >= I).
* i.e. p starts out false and then switches to true at I.
* Finds and returns I (the first number i such that p(i) is true),
* or hi if no such I exists (including if lo = hi).
*/
{
while (hi - lo > 1) {
int mid = lo + (hi - lo - 1) / 2;
if (p(mid)) hi = mid + 1; // mid >= I, search [lo, mid + 1)
else lo = mid + 1; // I > mid, search [mid + 1, hi)
}
return (lo != hi && p(lo)) ? lo : hi;
}
template<typename P>
double double_binsearch_last(P p, double lo, double hi, int n_it)
/* Takes a predicate P: double -> bool.
* Assumes that there is some X (lo <= X <= hi) such that p(x) = (x <= X),
* i.e. p starts out true and then switches to false after X.
* Finds and returns X (the last number x such that p(x) is true),
* or -infinity if no such X exists.
* This version runs a fixed number of iterations; to get results with a given
* precision, use `while (hi - lo > eps)` instead (but this is probably a bad
* idea because of rounding errors).
*/
{
while (n_it--) {
double mid = (lo + hi) / 2;
if (p(mid)) lo = mid; // mid <= X
else hi = mid; // X < mid
}
return p(lo) ? lo : -std::numeric_limits<double>::infinity();
}
template<typename P>
double double_binsearch_first(P p, double lo, double hi, int n_it)
/* Takes a predicate P: double -> bool.
* Assumes that there is some X (lo <= X <= hi) such that p(x) = (x >= X),
* i.e. p starts out false and then switches to true at X.
* Finds and returns X (the first number x such that p(x) is true),
* or +infinity if no such X exists.
* This version runs a fixed number of iterations; to get results with a given
* precision, use `while (hi - lo > eps)` instead (but this is probably a bad
* idea because of rounding errors).
*/
{
while (n_it--) {
double mid = (lo + hi) / 2;
if (p(mid)) hi = mid; // mid >= X
else lo = mid; // X > mid
}
return p(hi) ? hi : +std::numeric_limits<double>::infinity();
}

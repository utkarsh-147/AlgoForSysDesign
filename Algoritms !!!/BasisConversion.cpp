#include <algorithm>
#include <string>
#include <vector>
std::string const DIGITS = "0123456789ABCDEF";
unsigned long long int basis_string_to_number(std::string &s, int b) {
unsigned long long int result = 0ULL;
for (char d : s) {
result = b * result
+ (std::find(DIGITS.begin(), DIGITS.end(), d) - DIGITS.begin());
}
return result;
}
std::string number_to_basis_string(unsigned long long int n, int b) {
std::vector<char> ds;
do {
ds.push_back(DIGITS[n % b]);
n = n / b;
} while (n != 0);
return std::string(ds.rbegin(), ds.rend());
}

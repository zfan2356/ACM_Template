#include <bits/stdc++.h>

using namespace std;
using ll = long long;

template<typename A, typename B>
inline std::ostream &operator<<(std::ostream &out, const std::pair <A, B> &p) {
    return out << "(" << p.first << ", " << p.second << ")";
}

template<typename T>
inline std::ostream &operator<<(std::ostream &out, const std::vector <T> &a) {
    out << "[";
    for (int i = 0; i < a.size(); i++) {
        if (i) out << ',';
        out << ' ' << a[i];
    }
    return out << " ]";
}

template<typename T>
inline std::ostream &operator<<(std::ostream &out, const std::set <T> &a) {
    return out << std::vector<T>(all(a));
}

#define DEBUG

int gcd(int a, int b) {
    return !b ? a : gcd(b, a % b);
}

void solve() {

    cout << __gcd(-3, 6) << endl;
    cout << __gcd(3, -6) << endl;
}

signed main() {
    std::ios::sync_with_stdio(0);
    std::cin.tie(0);

#ifdef DEBUG
    freopen("D:\\mypile\\acmC++\\in.txt", "r", stdin);
    freopen("D:\\mypile\\acmC++\\out.txt", "w", stdout);
#endif

    int Case = 1;

    while (Case--) {
        solve();
    }

    return 0;
}
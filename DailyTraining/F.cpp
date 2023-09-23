#include <bits/stdc++.h>

using std::cin;
using std::cout;
using std::endl;
using i64 = long long;

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
    return out << std::vector<T>(a.begin(), a.end());
}

constexpr int inf = 1e9;

void solve() {


}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int Case = 1;

    while (Case--) {
        solve();
    }

    return 0;
}
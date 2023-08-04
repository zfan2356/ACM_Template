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
    return out << std::vector<T>(all(a));
}

constexpr i64 inf = 1e18;

void solve() {
    int n;
    cin >> n;
    std::vector<int> a(n);
    i64 sum = 0;
    for (int i = 0; i < n; i++) {
        cin >> a[i];
        sum += a[i];
    }
    i64 q = sum / n;

    auto get = [&](int x) {
        i64 res = 0, ans = 0;
        for (int i = 0; i < n; i++) {
            res += a[i] - x;
            ans += abs(a[i] - x);
        }
        if (abs(res) <= n / 2 && (ans - abs(res)) % 2 == 0) {
            if( (ans - abs(res)) % 2 == 0) {
                return ans - abs(res);
            } else if ((ans + abs(res)) % 2 == 0) {
                return ans + abs(res);
            } else {
                return inf;
            }
        } else {
            return inf;
        }
    };

    i64 ans = inf;
    for (int i = -1; i <= 1; i++) {
        ans = std::min(ans, get(i + q) / 2);
    }
    cout << ans << endl;

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
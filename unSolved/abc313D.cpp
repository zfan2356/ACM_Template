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

constexpr int inf = 1e9;

void solve() {
    int n, k;
    cin >> n >> k;

    auto query = [&](std::vector<int>& p) {
        cout << "? ";
        for (auto x : p) {
            cout << " " << x;
        }
        cout << endl;
        int x;
        cin >> x;
        return x;
    };

    std::vector<int> p;
    for (int i = 1; i <= k; i++) {
        p.push_back(i);
    }
    int t = query(p);
    std::vector<int> pre(n + 1, -1);

    for (int i = k + 1; i <= n; i++) {
        p.pop_back();
        p.push_back(i);
        int t2 = query(p);
        pre[i] = t ^ t2;
    }

    for (int i = 1; i < k; i++) {

    }


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
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

void solve() {
    string l, r;
    cin >> l >> r;
    reverse(l.begin(), l.end());
    reverse(r.begin(), r.end());
    while (l.size() < r.size()) l += '0';
    auto get = [&](string a, string b) {
        int ans = 0;
        for (int i = 0; i < a.size(); i++) {
            ans += abs(a[i] - b[i]);
        }
        return ans;
    };
    int res = get(l, r);

    for (int i = r.size() - 1; i >= 0; i--) {
        if (l[i] == r[i]) {
            continue;
        }
        int ans = r[i] - l[i];
        for (int j = i - 1; j >= 0; j--) {
            ans += 9;
        }
        res = max(res, ans);
        break;
    }
    cout << res << '\n';

}

signed main() {
    std::ios::sync_with_stdio(0);
    std::cin.tie(0);

#ifdef DEBUG
    freopen("D:\\mypile\\acm\\ICPC\\in.txt", "r", stdin);
    freopen("D:\\mypile\\acm\\ICPC\\out.txt", "w", stdout);
#endif

    int Case;
    std::cin >> Case;

    while (Case--) {
        solve();
    }

    return 0;
}
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

#define DEBUG

struct E{
    std::vector<std::pair<i64, int>> s;
    int t;
};

void solve() {
    int n, q;
    cin >> n >> q;
    std::vector<E> e;

    auto get = [&](std::pair<i64, i64> p) {
        std::vector<std::pair<i64, int>> s;
        while (1) {
            auto [x, y] = p;
            if (x == y) {
                break;
            } else if (x > y) {
                i64 X = x;
                x %= y;
                if (!x) {
                    x += y;
                }
                s.push_back({(X - x) / y, -1});
            } else {
                i64 Y = y;
                y %= x;
                if (!y) {
                    y += x;
                }
                s.push_back({(Y - y) / x, -2});
            }
            p = {x, y};
        }
        s.push_back({p.first, 0});
        reverse(s.begin(), s.end());
        return s;
    };

    for (int i = 0; i < n; i++) {
        i64 a, b;
        cin >> a >> b;
        get({a, b});
    }

    std::vector<int> ans(q);
    for (int i = 0; i < q; i++) {
        i64 a, b;
        cin >> a >> b;
        auto V = get({a, b});
        e.push_back({V, -i});
        V.push_back({0, 0});
        e.push_back({V, i});
    }

    auto compare = [&](std::vector<std::pair<i64, int>> a, std::vector<std::pair<i64, int>> b) {
        if (a[0].first != b[0].first) {
            return a[0] < b[0] ? -1 : 1;
        }
        for (int i = 1; i < std::min((int)a.size(), (int)b.size()); i++) {
            if (a[i].second != b[i].second) {
                return a[i].second < b[i].second ? -1 : 1;
            }
            if (!a[i].second) {
                return 0;
            }
            if (a[i].first < b[i].first) {
                if (i + 1 > (int)a.size()) {
                    return 1;
                }
                return a[i + 1].second < b[i].second ? -1 : 1;
            }
        }
    };

    std::sort(e.begin(), e.end(), [&](E a, E b) {
        int sign = compare(a.s, b.s);
        if (sign) {
            return sign < 0;
        }
        return a.t < b.t;
    });

}

signed main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

#ifdef DEBUG
    freopen("D:\\mypile\\acm\\ICPC\\in.txt", "r", stdin);
    freopen("D:\\mypile\\acm\\ICPC\\out.txt", "w", stdout);
#endif


    int Case = 1;

    while (Case--) {
        solve();
    }

    return 0;
}
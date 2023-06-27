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
template<class Info>
struct SegmentTree {
    int n;
    std::vector<Info> info;
    SegmentTree() : n(0) {}
    SegmentTree(int n_, Info v_ = Info()) {
        init(n_, v_);
    }
    template<class T>
    SegmentTree(std::vector<T> init_) {
        init(init_);
    }
    void init(int n_, Info v_ = Info()) {
        init(std::vector(n_, v_));
    }
    template<class T>
    void init(std::vector<T> init_) {
        n = init_.size();
        info.assign(4 << std::__lg(n), Info());
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                info[p].x = init_[l];
                return;
            }
            int m = (l + r) / 2;
            build(2 * p, l, m);
            build(2 * p + 1, m, r);
            pull(p);
        };
        build(1, 0, n);
    }
    void pull(int p) {
        info[p] = info[2 * p] + info[2 * p + 1];
    }
    void modify(int p, int l, int r, int x, const Info &v) {
        if (r - l == 1) {
            info[p] = v;
            return;
        }
        int m = (l + r) / 2;
        if (x < m) {
            modify(2 * p, l, m, x, v);
        } else {
            modify(2 * p + 1, m, r, x, v);
        }
        pull(p);
    }
    void modify(int p, const Info &v) {
        modify(1, 0, n, p, v);
    }
    Info rangeQuery(int p, int l, int r, int x, int y) {
        if (l >= y || r <= x) {
            return Info();
        }
        if (l >= x && r <= y) {
            return info[p];
        }
        int m = (l + r) / 2;
        return rangeQuery(2 * p, l, m, x, y) + rangeQuery(2 * p + 1, m, r, x, y);
    }
    Info rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }
    template<class F>
    int findFirst(int p, int l, int r, int x, int y, F pred) {
        if (l >= y || r <= x || !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        int res = findFirst(2 * p, l, m, x, y, pred);
        if (res == -1) {
            res = findFirst(2 * p + 1, m, r, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findFirst(int l, int r, F pred) {
        return findFirst(1, 0, n, l, r, pred);
    }
    template<class F>
    int findLast(int p, int l, int r, int x, int y, F pred) {
        if (l >= y || r <= x || !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        int res = findLast(2 * p + 1, m, r, x, y, pred);
        if (res == -1) {
            res = findLast(2 * p, l, m, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findLast(int l, int r, F pred) {
        return findLast(1, 0, n, l, r, pred);
    }
};

constexpr int inf = 1E9;
struct Max {
    int x = -1;
};

Max operator+(Max a, Max b) {
    return {std::max(a.x, b.x)};
}

#define DEBUG

void solve() {
    vector<int> A = {1, 9, 3, 5, 8, 2, 4};
    SegmentTree<Max> s(A);
    int k;
    cin >> k;
    cout << s.findLast(0, 7, [&](Max x) {
        return x.x >= k;
    });

}

signed main() {
    std::ios::sync_with_stdio(0);
    std::cin.tie(0);

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
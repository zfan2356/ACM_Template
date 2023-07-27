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
                info[p] = init_[l];
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

struct Info {
    int cnt = 0;
    i64 sum = 0;
};

Info operator+ (Info a, Info b) {
    return {a.cnt + b.cnt, a.sum + b.sum};
}

#define DEBUG

void solve() {
    int n, q, s;
    cin >> n >> q >> s;
    std::vector<int> l(n), r(n);
    std::vector<i64> d(n), a(n), b(n);
    std::vector<i64> all;
    for (int i = 0; i < n; i++) {
        cin >> a[i] >> b[i] >> l[i] >> r[i];
        d[i] = b[i] - a[i];
        all.push_back(d[i]);
    }

    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());

    std::vector<std::vector<std::pair<int, int>>> query(n);
    std::vector<i64> ans(q);
    for (int i = 0; i < q; i++) {
        int x, y;
        cin >> x >> y;
        x--;
        query[x].push_back({y, i});
    }

    SegmentTree<Info> segmentTree(n);
    int L = s, R = s;
    bool ok = true;
    i64 pre = 0;

    auto delMin = [&]() {
        auto f = [&](Info info) {
            return info.cnt > 0;
        };
        int p = segmentTree.findFirst(0, n, f);
        segmentTree.modify(p, {0, 0});
        return all[p];
    };

    auto delMax = [&]() {
        auto f = [&](Info info) {
            return info.cnt > 0;
        };
        int p = segmentTree.findLast(0, n, f);
        segmentTree.modify(p, {0, 0});
    };

    auto get = [&] (int l, int r) {
        if (R < l || L > r) {
            ok = false;
            return ;
        }
        for (;  !(l <= L && L <= r); L += 2) {
            pre += delMin();
        }
        for (; !(l <= R && R <= r); R -= 2) {
            delMax();
        }
    };

    auto Query = [&](int x) {
        if (!ok || x < L || x > R || (x - L) & 1) {
            return -1ll;
        }
        if (x == L) {
            return pre;
        }

        int p = segmentTree.findLast(0, n, [&](Info info) {
            return info.cnt <= (x - L) / 2;
        });

//        cout << "!!!" << p << endl;
        auto info = segmentTree.rangeQuery(0, p + 1);
        auto P = segmentTree.rangeQuery(p, p + 1);

        i64 ans = info.sum;
        if (info.cnt > (x - L) / 2) {
            ans -= 1ll * (info.cnt - (x - L) / 2) * (P.sum / P.cnt);
        }
        return ans + pre;
    };

    for (int i = 0; i < n; i++) {
        if (ok) {
            pre += a[i];
            int c = std::lower_bound(all.begin(), all.end(), d[i]) - all.begin();
            segmentTree.modify(c, {1, d[i]});
            L--, R++;
            get(l[i], r[i]);
        }
        for (auto [y, id] : query[i]) {
            ans[id] = Query(y);
        }
    }

    for (int i = 0; i < q; i++) {
        cout << ans[i] << '\n';
    }

}

signed main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

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
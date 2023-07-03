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

struct HLD {
    int n;
    std::vector<int> siz, top, dep, parent, in, out, seq;
    std::vector<std::vector<int>> adj;
    int cur;

    HLD() {}
    HLD(int n) {
        init(n);
    }

    void init(int n) {
        this->n = n;
        siz.resize(n);
        top.resize(n);
        dep.resize(n);
        parent.resize(n);
        in.resize(n);
        out.resize(n);
        seq.resize(n);
        cur = 0;
        adj.assign(n, {});
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    void work(int root = 0) {
        top[root] = root;
        dep[root] = 0;
        parent[root] = -1;
        dfs1(root);
        dfs2(root);
    }

    void dfs1(int u) {
        if (parent[u] != -1) {
            adj[u].erase(std::find(adj[u].begin(), adj[u].end(), parent[u]));
        }

        siz[u] = 1;
        for (auto &v : adj[u]) {
            parent[v] = u;
            dep[v] = dep[u] + 1;
            dfs1(v);
            siz[u] += siz[v];
            if (siz[v] > siz[adj[u][0]]) {
                std::swap(v, adj[u][0]);
            }
        }
    }

    void dfs2(int u) {
        in[u] = cur++;
        seq[in[u]] = u;
        for (auto v : adj[u]) {
            top[v] = v == adj[u][0] ? top[u] : v;
            dfs2(v);
        }
        out[u] = cur;
    }

    int lca(int u, int v) {
        while (top[u] != top[v]) {
            if (dep[top[u]] > dep[top[v]]) {
                u = parent[top[u]];
            } else {
                v = parent[top[v]];
            }
        }
        return dep[u] < dep[v] ? u : v;
    }

    int dist(int u, int v) {
        return dep[u] + dep[v] - 2 * dep[lca(u, v)];
    }

    int jump(int u, int k) {
        if (dep[u] < k) {
            return -1;
        }

        int d = dep[u] - k;

        while (dep[top[u]] > d) {
            u = parent[top[u]];
        }

        return seq[in[u] - dep[u] + d];
    }

    bool isAncester(int u, int v) {
        return in[u] <= in[v] && in[v] < out[u];
    }

    int rootedParent(int u, int v) {
        std::swap(u, v);
        if (u == v) {
            return u;
        }
        if (!isAncester(u, v)) {
            return parent[u];
        }
        auto it = std::upper_bound(adj[u].begin(), adj[u].end(), v, [&](int x, int y) {
            return in[x] < in[y];
        }) - 1;
        return *it;
    }

    int rootedSize(int u, int v) {
        if (u == v) {
            return n;
        }
        if (!isAncester(v, u)) {
            return siz[v];
        }
        return n - siz[rootedParent(u, v)];
    }

    int rootedLca(int a, int b, int c) {
        return lca(a, b) ^ lca(b, c) ^ lca(c, a);
    }
};

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

constexpr int inf = 1E9;
struct Min {
    int x = inf;
};

Min operator+(Min a, Min b) {
    return {std::min(a.x, b.x)};
}

struct Max {
    int x = -inf;
};

Max operator+(Max a, Max b) {
    return {std::max(a.x, b.x)};
}

struct Sum {
    int x = 0;
};

Sum operator+(Sum a, Sum b) {
    return {a.x + b.x};
}



//#define DEBUG

struct PPP{
    int x, y, u, v, i;
};

struct PP{
    int u, v, c, w;
};

void solve() {
    int n, q;
    cin >> n >> q;
    HLD Hld(n);
    vector<Sum> W(n);
    vector<vector<int>> Color(n);
    vector<PP> edge;

    for (int i = 1; i < n; i++) {
        int u, v, w, c;
        cin >> u >> v >> c >> w;
        u--, v--;
        c--;
        if (u > v) {
            swap(u, v);
        }
        edge.push_back({u, v, c, w});
        Hld.addEdge(u, v);
        Color[c].push_back(v);
    }

    Hld.work(0);
    for (auto [u, v, c, w] : edge) {
        W[Hld.in[v]].x = w;
    }


    SegmentTree<Sum> seg1(W), seg2(n);
    vector<PPP> Q(q);
    vector<vector<int>> colorQ(n);
    for (int i = 0; i < q; i++) {
        int x, y, u, v;
        cin >> x >> y >> u >> v;
        u--, v--;
        x--;
        Q[i] = {x, y, u, v, i};
        colorQ[x].push_back(i);
    }


    auto get = [&](int u, int v, int y) {
        int ans = 0;
        while (Hld.top[u] != Hld.top[v]) {
            if (Hld.dep[Hld.top[u]] < Hld.dep[Hld.top[v]]) {
                swap(u, v);
            }

            ans += seg1.rangeQuery(Hld.in[Hld.top[u]], Hld.in[u] + 1).x;
            ans += seg2.rangeQuery(Hld.in[Hld.top[u]], Hld.in[u] + 1).x * y;
//            cout << u << ' ' << Hld.parent[Hld.top[u]] << endl;
            u = Hld.parent[Hld.top[u]];
        }

        if (Hld.dep[u] < Hld.dep[v]) {
            swap(u, v);
        }

//        cout << Hld.in[u] << ' ' << Hld.in[v] << endl;
        ans += seg1.rangeQuery(Hld.in[v], Hld.in[u] + 1).x;
        ans += seg2.rangeQuery(Hld.in[v], Hld.in[u] + 1).x * y;

        ans -= seg1.rangeQuery(Hld.in[v], Hld.in[v] + 1).x + seg2.rangeQuery(Hld.in[v], Hld.in[v] + 1).x * y;
        return ans;
    };

    vector<int> ans(q);
    for (int i = 0; i < n - 1; i++) {
        for (auto j : Color[i]) {
            seg1.modify(Hld.in[j], {0});
            seg2.modify(Hld.in[j], {1});
        }

        for (auto j : colorQ[i]) {
            auto [x, y, u, v, id] = Q[j];

            ans[id] = get(u, v, y);
        }

        for (auto j : Color[i]) {
            seg1.modify(Hld.in[j], {W[Hld.in[j]]});
            seg2.modify(Hld.in[j], {0});
        }
    }

    for (int i = 0; i < q; i++) {
        cout << ans[i] << '\n';
    }

}

signed main() {
    std::ios::sync_with_stdio(0);
    std::cin.tie(0);

    int Case = 1;

    while (Case--) {
        solve();
    }

    return 0;
}
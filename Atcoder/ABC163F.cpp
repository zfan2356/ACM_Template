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
template <typename T>
struct Fenwick {
    int n;
    std::vector<T> a;

    Fenwick(int n = 0) {
        init(n);
    }

    Fenwick(vector<T>& b) {
        init(b);
    }

    void init(int n) {
        this->n = n;
        a.assign(n, T());
    }

    void init(vector<T> &b) {
        this->n = b.size();
        a.assign(n, T());
        for (int i = 1; i <= n; i++) {
            a[i - 1] += b[i - 1];
            int j = i + (i & -i);
            if (j <= n) a[j - 1] += a[i - 1];
        }
    }

    void add(int x, T v) {
        for (int i = x + 1; i <= n; i += i & -i) {
            a[i - 1] += v;
        }
    }

    T sum(int x) {
        auto ans = T();
        for (int i = x; i > 0; i -= i & -i) {
            ans += a[i - 1];
        }
        return ans;
    }

    T rangeSum(int l, int r) {
        return sum(r) - sum(l);
    }

    int kth(T k) {
        int x = 0;
        for (int i = 1 << std::__lg(n); i; i /= 2) {
            if (x + i <= n && k >= a[x + i - 1]) {
                x += i;
                k -= a[x - 1];
            }
        }
        return x;
    }
};

/*
 * https://atcoder.jp/contests/abc163/tasks/abc163_f
 * 首先要考虑到分类, 要求的是经过k的点, 那么我们可以对路径的端点进行分类, 首先就是至少一个端点是k的路径数目, 可以直接统计出来
 * 其次就是两个端点均不为k, 但是经过k的数目, 这些路径的计数是重点
 * 我们对每一个k点考虑, 可以把路径分为两类, 一类是该点子树内取两点, 一类是一个点取自子树, 一个点取自子树外, 这两种路径均恰好经过当前点k
 * 但是如果直接计数显然会有重复, 那么该如何规定次序来避免重复呢?
 * 我们首先处理比较靠下的子树, 然后依次往上计数(直接排序)
 * 然后对于一颗子树, 我们计数完之后将其中子树内的点都打一个标记, 代表着这些点我们不会再次取作端点, 这样的话我们第一类路径就不会计数重复
 * 至于第二类路径本来就不会重复, 自然就不需要
 * 将树压缩成一维, 利用树状数组来维护, 就可以实现上述操作(学习), 因为我们要对一颗子树进行信息查询, 我们需要将一颗子树映射为一段连续的区域
 *
 *
 */
void solve() {
    int n;
    cin >> n;
    vector<vector<int>> g(n), c(n);
    for (int i = 0; i < n; i++) {
        int x;
        cin >> x;
        x--;
        c[x].push_back(i);
    }

    for (int i = 1; i < n; i++) {
        int u, v;
        cin >> u >> v;
        u--, v--;
        g[u].push_back(v);
        g[v].push_back(u);
    }

    vector<int> fa(n), dep(n), sz(n), dfn(n);
    int cnt = 0;
    function<void(int, int)> dfs = [&](int u, int p) {
        sz[u] = 1;
        dfn[u] = cnt++;
        fa[u] = p;
        if(u) {
            dep[u] = dep[p] + 1;
        }
        for (auto v : g[u]) {
            if (v == p) {
                continue;
            }
            dfs(v, u);
            sz[u] += sz[v];
        }
    };

    dfs(0, -1);
    Fenwick<int> T(n);

    auto get = [&](int col) {
        if (c[col].empty()) {
            cout << 0 << '\n';
            return ;
        }

        sort(c[col].begin(), c[col].end(), [&](int x, int y) {
            return dep[x] > dep[y];
        });

        vector<pair<int, int>> oper;
        int last = c[col].size();
        ll ans = 1ll * last * n - 1ll * (last - 1) * last / 2;
        for (auto u : c[col]) {
            int y = T.rangeSum(dfn[u] + 1, dfn[u] + sz[u]);
            int q = T.sum(n) - y;
            last--;
            y = sz[u] - y - 1;
            ans += 1ll * y * (n - sz[u] - last - q);
            ll tot = 0;
            for (auto v : g[u]) {
                if (v == fa[u]) {
                    continue;
                }
                int g = T.rangeSum(dfn[v], dfn[v] + sz[v]);
                g = sz[v] - g;
                ans += 1ll * g * tot;
                tot += g;
            }
            oper.push_back({dfn[u], y + 1});
            T.add(dfn[u], y + 1);
        }
        for (auto [x, y] : oper) {
            T.add(x, -y);
        }
        cout << ans << '\n';
    };

    for (int i = 0; i < n; i++) {
        get(i);
    }
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
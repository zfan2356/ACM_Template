```c++
vector<int> dep(n, 0x3f3f3f3f);
vector<vector<int>> fa(n, vector<int>(21));
auto bfs = [&](int root) {
    queue<int> q;
    q.push(root);
    dep[root] = 1;
    while (q.size()) {
        int u = q.front(); q.pop();
        for (auto v : g[u]) {
            if (dep[v] > dep[u] + 1) {
                dep[v] = dep[u] + 1;
                q.push(v);
                fa[v][0] = u;
                for (int k = 1; k <= 20; k++) {
                    fa[v][k] = fa[fa[v][k - 1]][k - 1];
                }
            }
        }
    }
};

auto lca = [&](int a, int b) {
    if (dep[a] < dep[b]) swap(a, b);
    for (int k = 20; k >= 0; k--) {
        if (dep[fa[a][k]] >= dep[b]) a = fa[a][k];
    }
    if (a == b) return a;
    for (int k = 20; k >= 0; k--) {
        if (fa[a][k] != fa[b][k]) {
            a = fa[a][k];
            b = fa[b][k];
        }
    }
    return fa[a][0];
};

auto dis = [&](int a, int b) {
    return dep[a] + dep[b] - 2 * dep[lca(a, b)];
};
```
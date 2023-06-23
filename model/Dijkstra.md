```c++
auto dijkstra = [&](int MAXN, int s) {
    vector<int> d(MAXN, 1e9);
    vector<bool> vis(MAXN);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> h;
    h.push({0, s});
    d[s] = 0;
    while (h.size()) {
        auto [dis, u] = h.top(); h.pop();

        if (vis[u]) continue;
        vis[u] = true;

        for (auto [v, w] : g[u]) {
            if (d[v] > d[u] + w) {
                d[v] = d[u] + w;
                h.push({d[v], v});
            }
        }
    }
    return d;
};
```
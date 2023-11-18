```c++
auto dijkstra = [&](int MAXN, int s) {
    std::vector<int> d(MAXN, 1e9);
    std::vector<bool> vis(MAXN);
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> h;
    h.push({0, s});
    d[s] = 0;
    while (h.size()) {
        auto [dis, u] = h.top(); h.pop();

        if (vis[u]) continue;
        vis[u] = true;

        for (auto [v, w] : adj[u]) {
            if (d[v] > d[u] + w) {
                d[v] = d[u] + w;
                h.push({d[v], v});
            }
        }
    }
    return d;
};
```

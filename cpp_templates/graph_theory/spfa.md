
```c++
auto spfa = [&]() {
    std::queue<int> q;
    std::vector<bool> vis(n, true);
    std::vector<int> cnt(n), dist(n);
    for (int i = 0; i < n; i++) {
        q.push(i);
    }

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        vis[u] = false;

        for (auto [v, w] : g[u]) {
            if (dist[v] > dist[u] + w) {
                dist[v] = dist[u] + w;
                cnt[v] = cnt[u] + 1;
                if (cnt[v] >= n) {
                    return true;
                }

                if (!vis[v]) {
                    vis[v] = true;
                    q.push(v);
                }
            }
        }
    }
    return false;
};
```

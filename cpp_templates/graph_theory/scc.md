```c++
struct SCC {
    int n;
    std::vector<std::vector<int>> adj;
    std::vector<int> dfn, low, stk, id, siz;
    std::vector<bool> in_stk;
    std::vector<std::pair<int, int>> edges;
    int cnt, cur;

    SCC() {}
    SCC(int n) {
        init(n);
    }

    void init(int n) {
        this->n = n;
        adj.assign(n, {});
        dfn.assign(n, -1);
        low.resize(n);
        id.resize(n);
        siz.resize(n);
        in_stk.assign(n, false);
        stk.clear();
        cnt = cur = 0;
        edges.clear();
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
    }

    void dfs(int u) {
        stk.push_back(u);
        in_stk[u] = true;
        dfn[u] = low[u] = cur++;

        for (auto v : adj[u]) {
            if (dfn[v] == -1) {
                dfs(v);
                low[u] = std::min(low[u], low[v]);
            } else if (in_stk[v]) {
                low[u] = std::min(low[u], dfn[v]);
            }
        }

        if (dfn[u] == low[u]) {
            int y;
            do {
                y = stk.back();
                stk.pop_back();
                in_stk[y] = false;
                id[y] = cnt;
                siz[cnt]++;
            } while (y != u);
            cnt++;
        }
    }

    std::pair<int, std::vector<std::pair<int, int>>> work() {
        for (int i = 0; i < n; i++) {
            if (dfn[i] == -1) {
                stk.clear();
                dfs(i);
            }
        }

        for (int u = 0; u < n; u++) {
            for (auto v : adj[u]) {
                if (id[v] != id[u]) {
                    edges.push_back({id[u], id[v]});
                }
            }
        }

        return {cnt, edges};
    }
};
```

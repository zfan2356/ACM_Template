```c++
template<int LENGTH = 21>
struct LCA {
    int root;
    std::vector<int> dep;
    std::vector<std::array<int, LENGTH>> fa;
    std::vector<std::vector<int>> g;
    LCA() {}
    LCA(int N, int Root) {
        dep.assign(N, 0x3f3f3f3f);
        g.resize(N);
        fa.resize(N);
        this->root = Root;
    }

    void addEdge(int u, int v) {
        g[u].push_back(v);
        g[v].push_back(u);
    }

    void work() {
        auto bfs = [&](int root) {
            std::queue<int> q;
            q.push(root);
            dep[root] = 1;
            while (q.size()) {
                int u = q.front(); q.pop();
                for (auto v : g[u]) {
                    if (dep[v] > dep[u] + 1) {
                        dep[v] = dep[u] + 1;
                        q.push(v);
                        fa[v][0] = u;
                        for (int k = 1; k < LENGTH; k++) {
                            fa[v][k] = fa[fa[v][k - 1]][k - 1];
                        }
                    }
                }
            }
        };
        bfs(root);
    }

    int lca(int a, int b) {
        if (dep[a] < dep[b]) {
            std::swap(a, b);
        }
        for (int k = LENGTH - 1; k >= 0; k--) {
            if (dep[fa[a][k]] >= dep[b]) {
                a = fa[a][k];
            }
        }
        if (a == b) return a;
        for (int k = LENGTH - 1; k >= 0; k--) {
            if (fa[a][k] != fa[b][k]) {
                a = fa[a][k];
                b = fa[b][k];
            }
        }
        return fa[a][0];
    }
    
    int dist(int a, int b) {
        return dep[a] + dep[b] - 2 * dep[lca(a, b)];
    }
};
```
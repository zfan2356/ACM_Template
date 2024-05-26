#### 1. tarjan求割点

```c++
int root;
auto tarjan = [&](auto self, int u) -> void {
    dfn[u] = low[u] = cur++;
    int deg = 0;
    for (auto v : adj[u]) {
        if (dfn[v] == -1) {
            ++deg;
            self(self, v);
            low[u] = std::min(low[u], low[v]);
            if (low[v] >= dfn[u] && u != root) {
                buc[u] = true;
            }
        } else low[u] = std::min(low[u], dfn[v]);
    }
    if (u == root && deg >= 2) {
        buc[u] = true;
    }
};
tarjan(tarjan, root);
```

#### 2. tarjan求桥
其实只需要改一个地方即可, 并且也不需要判断根的情况
```c++
auto tarjan = [&](auto self, int u) -> void {
    dfn[u] = low[u] = cur++;
    for (auto v : adj[u]) {
        if (dfn[v] == -1) {
            self(self, v);
            low[u] = std::min(low[u], low[v]);
            if (low[v] > dfn[u]) {
                bridge[v] = true;
            }
        } else low[u] = std::min(low[u], dfn[v]);
    }
};
```

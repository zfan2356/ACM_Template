### 树链剖分
```go
// HLD ————————————————————————————
type HLD struct {
    cur int
    siz []int
    top []int
    dep []int
    parent []int
    in []int
    out []int
    seq []int
    adj [][]int
}

func (h *HLD) Init(N int)  {
    h.cur = 0
    h.seq = make([]int, N)
    h.siz = make([]int, N)
    h.top = make([]int, N)
    h.dep = make([]int, N)
    h.parent = make([]int, N)
    h.in = make([]int, N)
    h.out = make([]int, N)
    h.adj = make([][]int, N)
}

func (h *HLD) addEdge(u int, v int)  {
    h.adj[u] = append(h.adj[u], v)
    h.adj[v] = append(h.adj[v], u)
}

func (h *HLD) dfs1(u int) {
    if h.parent[u] != -1 {
        idx := 0
        for i, v := range h.adj[u] {
            if v == h.parent[u] {
                idx = i
                break
            }
        }
        h.adj[u] = append(h.adj[u][:idx], h.adj[u][idx+1:]...)
    }

    h.siz[u] = 1
    for idx, v := range h.adj[u] {
        h.parent[v] = u
        h.dep[v] = h.dep[u] + 1
        h.dfs1(v)
        h.siz[u] += h.siz[v]
        if h.siz[v] > h.siz[h.adj[u][0]] {
            h.adj[u][idx], h.adj[u][0] = h.adj[u][0], h.adj[u][idx]
        }
    }
}

func (h *HLD) dfs2(u int)  {
    h.in[u] = h.cur
    h.cur++
    h.seq[h.in[u]] = u
    for _, v := range h.adj[u] {
        if v == h.adj[u][0] {
            h.top[v] = h.top[u]
        } else {
            h.top[v] = v
        }
        h.dfs2(v)
    }
    h.out[u] = h.cur
}

func (h *HLD) work(root int)  {
    h.top[root] = root
    h.dep[root] = 0
    h.parent[root] = -1
    h.dfs1(root)
    h.dfs2(root)
}

func (h *HLD) lca(u int, v int) int {
    for h.top[u] != h.top[v] {
        if h.dep[h.top[u]] > h.dep[h.top[v]] {
            u = h.parent[h.top[u]]
        } else {
            v = h.parent[h.top[v]]
        }
    }
    if h.dep[u] < h.dep[v] {
        return u
    }
    return v
}

func (h *HLD) dist(u int, v int) int {
    return h.dep[u] + h.dep[v] - 2 * h.dep[h.lca(u, v)]
}

```
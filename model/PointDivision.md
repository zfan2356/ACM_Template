```c++
int main() {
    int n;
    std::vector<std::vector<int>> g(n);
    std::vector<int> a(n); //节点信息
    std::vector<int> siz(n);
    std::vector<bool> vis(n, false);
    int root = -1, rootSiz = 0x3f3f3f3f;
    
    // 选取重心
    std::function<void(int, int, int)> findRoot = [&](int u, int fa, int SIZ) {
        siz[u] = 1;
        int maxv = 0;
        for (auto v : g[u]) {
            if (v == fa || vis[v]) {
                continue;
            }
            findRoot(v, u, SIZ);
            siz[u] += siz[v];
            maxv = std::max(maxv, siz[v]);
        }
        maxv = std::max(maxv, SIZ - siz[u]);
        if (rootSiz > maxv) {
            rootSiz = maxv;
            root = u;
        }
    };
    
    
    //处理从u往下走的路径信息
    std::function<void(int, int, std::array<int, 2>)> dfs = [&](int u, int fa, int pre) {
        //更新答案
        for (auto v : g[u]) {
            if (v == fa || vis[v]) {
                continue;
            }
            dfs(v, u, pre + a[v]);
        }
    };
    
    // 处理以u为重心的子树, 经过u的路径
    auto calc = [&](int u) {
        // 清空一下
        for (auto v : g[u]) {
            if (vis[v]) {
                continue;
            }
            dfs(v, u, a[u] + a[v]); // 计算从儿子v往下走的路径信息
            // 进行合并
        }
    };
    
    
    //点分治, 传入初始点以及整棵树的大小
    std::function<void(int, int)> divide = [&](int u, int SIZ) {
        root = -1, rootSiz = 0x3f3f3f3f;
        findRoot(u, -1, SIZ);
        vis[root] = true;
        calc(root);
        for (auto v : g[root]) {
            if (vis[v]) {
                continue;
            }
            divide(v, siz[v]);
        }
    };
    divide(0, n);
    
    return 0;
}


```
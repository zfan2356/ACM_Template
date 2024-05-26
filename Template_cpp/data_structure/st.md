```c++
template <typename T>
struct ST {
    int n, m;
    std::vector<std::vector<T>> f;

    ST () {}
    ST (std::vector<T> _init) {
        n = _init.size();
        m = std::__lg(n) + 2;
        f.assign(n, std::vector<T>(m, inf));
        init(_init);
    }

    void init(std::vector<T> a) {
        for (int i = 0; i < n; i++) {
            f[i][0] = a[i];
        }
        for (int j = 1; j < m; j++) {
            for (int i = 0; i + (1 << j) - 1 < n; i++) {
                f[i][j] = std::min(f[i][j - 1], f[i + (1 << j - 1)][j - 1]);
            }
        }
    }

    T query(int l, int r) {
        int t = std::__lg(r - l);
        return std::min(f[l][t], f[r - (1 << t)][t]);
    }
};
```
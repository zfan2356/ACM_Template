```c++
template <typename T>
struct WeightDSU {
    std::vector<int> f;
    std::vector<T> d;

    WeightDSU() {}
    WeightDSU(int n) {
        init(n);
    }

    void init(int n) {
        f.assign(n, -1);
        d.assign(n, T());
    }

    int find(int x) {
        if (f[x] == -1) {
            return x;
        }
        int r = find(f[x]);
        d[x] += d[f[x]];
        f[x] = r;
        return r;
    }

    bool same(int x, int y) {
        return find(x) == find(y);
    }

    void merge(int x, int y, T w) {
        w += d[x], w -= d[y];
        x = find(x), y = find(y);
        f[y] = x;
        d[y] = w;
    }

    T diff(int x, int y) {
        find(x), find(y);
        return d[y] - d[x];
    }
};
```

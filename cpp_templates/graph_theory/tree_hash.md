```c++
static const ull mask = std::chrono::steady_clock::now().time_since_epoch().count();
vector<ull> hash(n);
auto shift = [&](ull x) {
    x ^= mask;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    x ^= mask;
    return x;
};

function<void(int, int)> getHash = [&](int u, int fa) {
    hash[u] = 1;
    for (auto v : g[u]) {
        if (v == fa) {
            continue;
        }
        getHash(v, u);
        hash[u] += shift(hash[v]);
    }
};

getHash(0, -1);
```

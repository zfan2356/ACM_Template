```c++
struct MuFunction {
    std::vector<int> mu, primes;
    std::vector<bool> vis;

    MuFunction() {}
    MuFunction(int n) {
        _init(n);
    }

    void _init(int N) {
        mu.resize(N + 1);
        vis.resize(N + 1);
        mu[1] = 1;
        for (int i = 2; i <= N; i++) {
            if (!vis[i]) {
                primes.push_back(i);
                mu[i] = -1;
            }
            for (int j = 0; primes[j] * i <= N; j++) {
                vis[primes[j] * i] = true;
                if (i % primes[j] == 0) {
                    break;
                }
                mu[primes[j] * i] = -mu[i];
            }
        }
    }
    
    int getMu(int x) {
        return mu[x];
    }
};
```

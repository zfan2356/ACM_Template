```c++
struct EulerFunction {
    vector<int> phi, primes;
    vector<bool> vis;

    EulerFunction() {}
    EulerFunction(int n) {
        init(n);
    }

    void init(int n) {
        phi.resize(n + 1);
        vis.resize(n + 1);
        phi[1] = 1;
        for (int i = 2; i <= n; i++) {
            if (!vis[i]) {
                primes.push_back(i);
                phi[i] = i - 1;
            }
            for (int j = 0; primes[j] * i <= n; j++) {
                vis[primes[j] * i] = true;
                if (i % primes[j] == 0) {
                    phi[i * primes[j]] = phi[i] * primes[j];
                    break;
                }
                phi[i * primes[j]] = phi[i] * (primes[j] - 1);
            }
        }
    }

    int getEuler(int x) {
        return phi[x];
    }
};
```

```c++
auto getEuler = [&](int x) {
    int res = x;
    for (int i = 2; i <= x / i; i++) {
        if (x % i == 0) {
            res = res / i * (i - 1);
            while (x % i ==0) {
                x /= i;
            }
        }
    }
    if (x > 1) {
        res = res / x * (x - 1);
    }
    return res;
};
```
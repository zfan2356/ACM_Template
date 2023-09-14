//线性筛欧拉函数
struct EulerFunction {
    std::vector<int> phi, primes;
    std::vector<bool> vis;

    EulerFunction() {}
    EulerFunction(int n) {
        _init(n);
    }

    void _init(int N) {
        phi.resize(N + 1);
        vis.resize(N + 1);
        phi[1] = 1;
        for (int i = 2; i <= N; i++) {
            if (!vis[i]) {
                primes.push_back(i);
                phi[i] = i - 1;
            }
            for (int j = 0; primes[j] * i <= N; j++) {
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


// 求单个欧拉函数
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


//Euler降幂
std::vector<int> phi;
auto getEuler = [&](int x) {
    int res = x;
    for (int i = 2; i <= x / i; i++) {
        if (x % i == 0) {
            res = res / i * (i - 1);
            while (x % i == 0) {
                x /= i;
            }
        }
    }
    if (x > 1) {
        res = res / x * (x - 1);
    }
    return res;
};

phi.push_back(MOD);
while (phi.back() != 1) {
    phi.push_back(getEuler(phi.back()));
}

auto EulerDrop = [&](int l, int r, int MOD) {
    r = std::min(r, l - 1 + (int)phi.size());
    ll res = 1;
    bool ok = true;

    auto qpow = [&](long long a, long long b, long long p) {
        long long res = 1 % p;
        for (; b; b >>= 1) {
            if (b & 1) {
                res = 1ll * res * a;
                if (res >= p) {
                    ok = true;
                    res %= p;
                }
            }
            a = 1ll * a * a;
            if (a >= p) {
                ok = true;
                a %= p;
            }
        }
        return res;
    };
    for (int i = r; i >= l; i--) {
        ok = false;
        res = qpow(a[i], res, phi[i - l]);
        if (ok) {
            res += phi[i - l];
        }
    }
    return res % phi[0];
};

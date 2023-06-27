#include <bits/stdc++.h>

using namespace std;
using ll = long long;

template<typename A, typename B>
inline std::ostream &operator<<(std::ostream &out, const std::pair <A, B> &p) {
    return out << "(" << p.first << ", " << p.second << ")";
}

template<typename T>
inline std::ostream &operator<<(std::ostream &out, const std::vector <T> &a) {
    out << "[";
    for (int i = 0; i < a.size(); i++) {
        if (i) out << ',';
        out << ' ' << a[i];
    }
    return out << " ]";
}

template<typename T>
inline std::ostream &operator<<(std::ostream &out, const std::set <T> &a) {
    return out << std::vector<T>(all(a));
}

#define DEBUG

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

void solve() {
    int n;
    cin >> n;

    vector<int> phi(n), primes;
    vector<bool> vis(n);
    phi[1] = 1;
    for (int i = 2; i < n; i++) {
        if (!vis[i]) {
            primes.push_back(i);
            phi[i] = i - 1;
        }
        for (int j = 0; primes[j] * i < n; j++) {
            vis[primes[j] * i] = true;
            if (i % primes[j] == 0) {
                phi[i * primes[j]] = phi[i] * primes[j];
                break;
            }
            phi[i * primes[j]] = phi[i] * (primes[j] - 1);
        }
    }

}

signed main() {
    std::ios::sync_with_stdio(0);
    std::cin.tie(0);

#ifdef DEBUG
    freopen("D:\\mypile\\acm\\ICPC\\in.txt", "r", stdin);
    freopen("D:\\mypile\\acm\\ICPC\\out.txt", "w", stdout);
#endif


    int Case = 1;

    while (Case--) {
        solve();
    }

    return 0;
}
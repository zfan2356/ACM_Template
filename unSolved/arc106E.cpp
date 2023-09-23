#include <bits/stdc++.h>

using std::cin;
using std::cout;
using std::endl;
using i64 = long long;

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
    return out << std::vector<T>(a.begin(), a.end());
}

constexpr int inf = 1e9;


void solve() {
    int n, k;
    cin >> n >> k;
    std::vector<int> a(n);
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }
    const int M = 2 * n * k;
    std::vector<int> vis(M + 1);
    for (int i = 1; i <= M; i++) {
        for (int j = 0; j < n; j++) {
            if ((i - 1) / a[j] % 2 == 0) {
                vis[i] |= (1 << j);
            }
        }
    }
    const int S = 1 << n;
    // 前x天作为右部, 这时二分图是否存在完美匹配, 左边有n * k个点,
    // 利用hall定理, 枚举左边的集合s, 然后和右边所有有连边的点的数目必须>= popcount(s)
    // 令dp[i] 代表点集的状态为i时, 和右边没有连边的点的个数
    auto check = [&](int x) {
        std::vector<int> dp(S);
        for (int i = 1; i <= x; i++) {
            dp[vis[i] ^ (S - 1)]++;
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < S; j++) {
                if (j >> i & 1)
            }
        }

        for (int i = 0; i < s; i++) {
            if (__builtin_popcount(i) * k > x - dp[S - 1 - i]) {
                return false;
            }
        }
        return true;
    };

    int l = 1, r = M;
    while (l < r) {
        int mid = (l + r) / 2;
        if (check(mid)) {
            r = mid;
        } else {
            l = mid + 1;
        }
    }
    cout << r << '\n';
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int Case = 1;

    while (Case--) {
        solve();
    }

    return 0;
}
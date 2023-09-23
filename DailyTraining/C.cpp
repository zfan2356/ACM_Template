#include <bits/stdc++.h>
using namespace std;
const int N = 60;
int g[N], c[N << 1], tt, a[N], b[N], inv[N], n, f[N][N << 1];
const int P = 998244353;
int power(int a, int b) {
    int c = 1;
    for (; b; b >>= 1) {
        if (b & 1) c = 1ll * c * a % P;
        a = 1ll * a * a % P;
    }
    return c;
}
int main() {
    scanf("%d", &n);
    for (int i = 1; i <= n; ++i) {
        scanf("%d%d", a + i, b + i);
        c[++tt] = a[i];
        c[++tt] = ++b[i];
    }
    sort(c + 1, c + tt + 1);
    tt = unique(c + 1, c + tt + 1) - c - 1;
    for (int i = 1; i <= n; ++i) {
        a[i] = lower_bound(c + 1, c + tt + 1, a[i]) - c;
        b[i] = lower_bound(c + 1, c + tt + 1, b[i]) - c;
    }
    for (int i = 1; i <= n; ++i) inv[i] = power(i, P - 2);
    a[0] = 1, b[0] = tt + 1;
    for (int i = 1; i <= tt; ++i) f[0][i] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = a[i]; j < b[i]; ++j) {
            int len = c[j + 1] - c[j];
            g[1] = len;
            for (int k = 2; k <= i; ++k) g[k] = 1ll * g[k - 1] * (len + k - 1) % P * inv[k] % P;
            for (int k = i - 1; k >= 0; --k) {
                f[i][j] += 1ll * f[k][j + 1] * g[i - k] % P;//C(len + i - k - 1, i - k)
                f[i][j] %= P;
                if (j < a[k] || j >= b[k]) break;
            }
        }
        for (int j = tt - 1; j; --j) {
            f[i][j] += f[i][j + 1];
            f[i][j] %= P;
        }
    }
    int ans = f[n][1];
    cout << ans << '\n';
    for (int i = 1; i <= n; ++i) ans = 1ll * ans * power(c[b[i]] - c[a[i]], P - 2) % P;
    printf("%d", ans);
}
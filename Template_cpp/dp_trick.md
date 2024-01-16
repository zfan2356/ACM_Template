### 数位dp常见套路

> 求有多少个小于 $N$ 的数, 数位和为 $k$, 且数位和可以整除当前的数

```c++
std::vector dp(15, std::vector(150, std::vector(150, -1LL)));
auto dfs = [&](auto self, int u, int s, int x, int lim) -> i64 {
    if (u == -1) {
        return s == k && x == 0;
    }
    if (!lim && dp[u][s][x] != -1) {
        return dp[u][s][x];
    }
    int up = lim ? digit[u] : 9;
    i64 res = 0;
    for (int i = 0; i <= up; i++) {
        res += self(self, u - 1, s + i, (x * 10 + i) % k, lim && (i == up));
    }
    if (!lim) dp[u][s][x] = res;
    return res;
};
```

$dp[u][s][x]$ 代表的是, 当前进行到第 $u$ 位, 和为 $s$, 并且和 $mod\ k = x$ 的数量

dfs 中 $limit$ 代表是否到达限制, $dijit[]$ 中存储的是 $N$ 的数位

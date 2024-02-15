## 扫描线

应用: 求 $n$ 个矩形的面积并

```c++
#include <bits/stdc++.h>

using i64 = long long;

constexpr int N = 2E5 + 10;

constexpr int inf = 1E9 + 10;

struct L {
    int l, r, h, v;
};
std::vector<int> all;
std::vector<L> line;

i64 sum[N * 4], mn[N * 4], tag[N * 4];

void pull(int p) {
    mn[p] = std::min(mn[p * 2], mn[p * 2 + 1]);
    if (mn[p] == mn[p * 2]) {
        sum[p] = sum[2 * p];
    } else {
        sum[p] = 0;
    }
    if (mn[p] == mn[2 * p + 1]) sum[p] += sum[2 * p + 1];
}

void build(int p, int l, int r) {
    if (r - l == 1) {
        sum[p] = all[r] - all[l];
        return ;
    }
    int m = (l + r) / 2;
    build(2 * p, l, m);
    build(2 * p + 1, m, r);
    pull(p);
}

void apply(int p) {
    if (tag[p] != 0) {
        tag[2 * p] += tag[p];
        tag[2 * p + 1] += tag[p];
        mn[2 * p] += tag[p];
        mn[2 * p + 1] += tag[p];
        tag[p] = 0;
    }
}

void rangeApply(int p, int l, int r, int x, int y, int v) {
    if (l >= y || r <= x) {
        return ;
    }
    if (x <= l && r <= y) {
        tag[p] += v;
        mn[p] += v;
        return ;
    }
    int m = (l + r) / 2;
    apply(p);
    rangeApply(2 * p, l, m, x, y, v);
    rangeApply(2 * p + 1, m, r, x, y, v);
    pull(p);
}

i64 rangeQuery(int p, int l, int r, int x, int y) {
    if (l >= y || r <= x) {
        return 0LL;
    }
    if (x <= l && r <= y) {
        i64 len = all[r] - all[l];
        return mn[p] ? len : len - sum[p];
    }
    int m = (l + r) / 2;
    apply(p);
    return rangeQuery(2 * p, l, m, x, y) + rangeQuery(2 * p + 1, m, r, x, y);
}

void solve() {
    int n;
    std::cin >> n;

    for (int i = 0; i < n; i++) {
        int x1, y1, x2, y2;
        std::cin >> x1 >> y1 >> x2 >> y2;
        all.push_back(x1);
        all.push_back(x2);
        line.push_back({x1, x2, y1, 1});
        line.push_back({x1, x2, y2, -1});
    }
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());

    std::sort(line.begin(), line.end(), [&](L a, L b) {
        return a.h < b.h;
    });

    int M = all.size();
    build(1, 0, M - 1);

    i64 ans = 0;
    for (int i = 1; i < n * 2; i++) {
        int x = std::lower_bound(all.begin(), all.end(), line[i - 1].l) - all.begin();
        int y = std::lower_bound(all.begin(), all.end(), line[i - 1].r) - all.begin();
        rangeApply(1, 0, M - 1, x, y, line[i - 1].v);
        ans += rangeQuery(1, 0, M - 1, 0, M - 1) * (line[i].h - line[i - 1].h);
    }

    std::cout << ans << '\n';


}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int T = 1;

    while (T--) {
        solve();
    }

    return 0;
}
```

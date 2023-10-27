#include <bits/stdc++.h>

using i64 = long long;

struct Node {
    int cnt = 0;
    Node *l = nullptr, *r = nullptr;
    Node() {}
};

void pull(Node *&t) {
    t->cnt = 0;
    if (t->l) t->cnt += t->l->cnt;
    if (t->r) t->cnt += t->r->cnt;
}

void modify(Node *&t, int l, int r, int x, int val) {
    if (t == nullptr) {
        t = new Node;
    }
    if (r - l == 1) {
        t->cnt += val;
        return ;
    }

    int m = (l + r) / 2;
    if (x < m) {
        modify(t->l, l, m, x, val);
    } else {
        modify(t->r, m, r, x, val);
    }
    pull(t);
}

int query(Node *t, int l, int r, int x, int y) {
    if (t == nullptr || l >= y || r <= x) {
        return 0;
    }
    if (x <= l && r <= y) {
        return t->cnt;
    }
    int m = (l + r) / 2;
    return query(t->l, l, m, x, y) + query(t->r, m, r, x, y);
}

constexpr int inf = 1e9;

struct T {
    int x, r, f;
};

void solve() {
    int n, k;
    std::cin >> n >> k;

    std::vector<T> a(n);
    std::vector<i64> all;
    for (int i = 0; i < n; i++) {
        int x, r, f;
        std::cin >> x >> r >> f;
        a[i] = {x, r, f};
        all.push_back(x), all.push_back(x + r), all.push_back(x - r);
    }

    std::sort(a.begin(), a.end(), [&](T a, T b) {
        return a.r > b.r;
    });
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());

    const int M = all.size(), N = 1e4 + 10;

    std::vector<Node *> seg(N);

    i64 ans = 0;
    for (int i = 0; i < n; i++) {
        auto [x, r, f] = a[i];
        int L = std::lower_bound(all.begin(), all.end(), x - r) - all.begin();
        int R = std::lower_bound(all.begin(), all.end(), x + r) - all.begin();
        x = std::lower_bound(all.begin(), all.end(), x) - all.begin();
//        std::cout << x << "\n";

        for (int j = std::max(1, a[i].f - k); j <= std::min(N - 1, a[i].f + k); j++) {
            if (seg[j]) {
                ans += query(seg[j], 0, M, L, R + 1);
//                std::cout << L << " " << R << '\n';
            }
        }
        modify(seg[a[i].f], 0, M, x, 1);
    }

//    std::cout << query(seg[1], 0, M, 7, 9) << '\n';

    std::cout << ans << '\n';
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
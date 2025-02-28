```c++
/*
 * 主席树模板, 支持历史版本单点修改, 历史版本单点查询
 */
struct Node {
    Node *l = nullptr;
    Node *r = nullptr;
    int v = 0;
    Node(Node *t) {
        if (t) {
            *this = *t;
        }
    }
};

Node *modify(Node *t, int l, int r, int x, int val) {
    t = new Node(t);
    if (r - l == 1) {
        t->v = val;
        return t;
    }
    int m = (l + r) / 2;
    if (x < m) {
        t->l = modify(t->l, l, m, x, val);
    } else {
        t->r = modify(t->r, m, r, x, val);
    }
    return t;
}

int query(Node *t, int l, int r, int x) {
    if (r - l == 1) {
        return t->v;
    }
    int m = (l + r) / 2;
    if (x < m) {
        return query(t->l, l, m, x);
    } else {
        return query(t->r, m, r, x);
    }
}

void solve() {
    int n, m;
    std::cin >> n >> m;

    std::vector<int> a(n);
    for (int i = 0; i < n; i++) {
        std::cin >> a[i];
    }

    std::vector<Node *> tree(m + 1);
    std::function<Node *(Node *, int, int)> build = [&](Node *t, int l, int r) {
        t = new Node(t);
        if (r - l == 1) {
            t->v = a[l];
            return t;
        }
        int m = (l + r) / 2;
        t->l = build(t->l, l, m), t->r = build(t->r, m, r);
        return t;
    };

    tree[0] = build(tree[0], 0, n);

    for (int i = 1; i <= m; i++) {
        int v, op, x, val;
        std::cin >> v >> op;

        if (op == 1) {
            std::cin >> x >> val;
            x--;
            tree[i] = modify(tree[v], 0, n, x, val);
        } else {
            std::cin >> x;
            x--;
            std::cout << query(tree[v], 0, n, x) << '\n';
            tree[i] = tree[v];
        }
    }
}
```

```c++
/*
 权值主席树: 静态区间第k小
    其本质是可以实现静态查询一个区间内小于某个数的值有多少个

其实权值线段树可以实现的操作, 如果套一个主席树, 就可以拓展到在线区间查询中使用
 */

struct Node {
    Node *l = nullptr;
    Node *r = nullptr;
    int cnt = 0;

    Node(Node *t) {
        if (t) {
            *this = *t;
        }
    }
};

Node *add(Node *t, int l, int r, int x) {
    t = new Node(t);
    t->cnt += 1;
    if (r - l == 1) {
        return t;
    }
    int m = (l + r) / 2;
    if (x < m) {
        t->l = add(t->l, l, m, x);
    } else {
        t->r = add(t->r, m, r, x);
    }
    return t;
}

int query(Node *u, Node *v, int l, int r, int k) {
    if (r - l == 1) {
        return l;
    }
    int m = (l + r) / 2;
    Node *lv = v->l;
    Node *lu = u->l;
    int sum = (lv->cnt) - (lu->cnt);
    if (sum >= k) {
        return query(u->l, v->l, l, m, k);
    } else {
        return query(u->r, v->r, m, r, k - sum);
    }
}

void solve() {
    int n, m;
    std::cin >> n >> m;

    std::vector<int> a(n);
    for (int i = 0; i < n; i++) {
        std::cin >> a[i];
    }

    std::vector<int> b(a);
    std::sort(b.begin(), b.end());
    b.erase(std::unique(b.begin(), b.end()), b.end());
    const int M = b.size();

    std::vector<Node *> tree(n + 1);
    std::function<Node *(Node *, int, int)> build = [&](Node *t, int l, int r) {
        t = new Node(t);
        if (r - l == 1) {
            return t;
        }
        int m = (l + r) / 2;
        t->l = build(t->l, l, m), t->r = build(t->r, m, r);
        return t;
    };
    tree[0] = build(tree[0], 0, M);

    for (int i = 1; i <= n; i++) {
        int A = std::lower_bound(b.begin(), b.end(), a[i - 1]) - b.begin();
        tree[i] = add(tree[i - 1], 0, M, A);
    }

    for (int i = 0; i < m; i++) {
        int l, r, k;
        std::cin >> l >> r >> k;

        int p = query(tree[l - 1], tree[r], 0, M, k);
        std::cout << b[p] << '\n';
    }
}
```

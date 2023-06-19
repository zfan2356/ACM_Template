```c++
template <typename T>
struct Fenwick {
    int n;
    std::vector<T> a;   //下标从0开始

    Fenwick(int n = 0) {
        init(n);
    }

    Fenwick(vector<T>& b) {
        init(b);
    }

    void init(int n) {
        this->n = n;
        a.assign(n, T());
    }

    void init(vector<T> &b) {   //根据b O(n) 建树
        this->n = b.size();
        a.assign(n, T());
        for (int i = 1; i <= n; i++) {
            a[i - 1] += b[i - 1];
            int j = i + (i & -i);
            if (j <= n) a[j - 1] += a[i - 1];
        }
    }

    void add(int x, T v) {
        for (int i = x + 1; i <= n; i += i & -i) {
            a[i - 1] += v;
        }
    }

    T sum(int x) {      //求[0, x)的前缀和
        auto ans = T();
        for (int i = x; i > 0; i -= i & -i) {
            ans += a[i - 1];
        }
        return ans;
    }

    T rangeSum(int l, int r) {  //求[l, r)的区间和
        return sum(r) - sum(l);
    }

    int kth(T k) {
        int x = 0;
        for (int i = 1 << std::__lg(n); i; i /= 2) {
            if (x + i <= n && k >= a[x + i - 1]) {
                x += i;
                k -= a[x - 1];
            }
        }
        return x;
    }
};
```

### 单点修改, 区间查询
```c++
void solve() {
    int n, m;
    cin >> n >> m;
    
    vector<int> a(n);
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }
    
    Fenwick<int> f(a);

    while (m--) {
        int t, x, k;
        cin >> t >> x >> k;
        if (t == 1) {
            x--;
            f.add(x, k);
        }
        else {
            x--, k--;
            cout << f.rangeSum(x, k + 1) << '\n';
        }
    }
}
```

### 单点查询, 区间修改
```c++
void solve() {
    int n, m;
    cin >> n >> m;

    vector<int> a(n);
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    for (int i = n - 1; i > 0; i--) {
        a[i] -= a[i - 1];
    }

    Fenwick<int> f(a);

    while (m--) {
        int t;
        cin >> t;
        if (t == 1) {
            int x, y, k;
            cin >> x >> y >> k;
            x--, y--;
            f.add(x, k), f.add(y + 1, -k);
        } else {
            int x;
            cin >> x;
            x--;
            cout << f.sum(x + 1) << '\n';
        }
    }
}
```

### 区间修改, 区间查询
```c++
void solve() {
    int n, m;
    cin >> n >> m;
    vector<int> a(n + 1), b(n + 1);
    for (int i = 1; i <= n; i++) {
        cin >> a[i];
    }
    for (int i = n; i; i--) {
        a[i] -= a[i - 1];
        b[i] = i * a[i];
    }

    Fenwick<int> tr1(a), tr2(b);

    while (m--) {
        char op[2];
        cin >> op;
        if (op[0] == 'C') {
            int x, y, v;
            cin >> x >> y >> v;
            tr1.add(x, v), tr1.add(y + 1, -v);
            tr2.add(x, v * x), tr2.add(y + 1, -(v * (y + 1)));
        }
        else {
            int x, y;
            cin >> x >> y;
            ll ans = (y + 1) * tr1.sum(y + 1) - x * tr1.sum(x);
            ans -= tr2.rangeSum(x, y + 1);
            cout << ans << endl;
        }
    }
}
```

### 查询全局第k小
```c++
void solve() {
    int n;
    cin >> n;

    vector<int> a(n);
    vector<int> all;
    for (int i = 0; i < n; i++) {
        cin >> a[i];
        all.push_back(a[i]);
    }
    sort(all.begin(), all.end());
    all.erase(unique(all.begin(), all.end()), all.end());

    for (int i = 0; i < n; i++) {
        a[i] = lower_bound(all.begin(), all.end(), a[i]) - all.begin();
    }
    int N = all.size();
    Fenwick<int> f(N);

    for (int i = 0; i < n; i++) {
        f.add(a[i], 1);
    }

    int m; cin >> m;
    for (int i = 0; i < m; i++) {
        int x;
        cin >> x;
        f.kth(x);   // 第x小的数, 对应值为all[f.kth(x) - 1];
    }

}
```

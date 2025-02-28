```c++
template<typename T>
struct Splay {
    int root, tot;
    std::vector<std::array<int, 2>> tr;
    std::vector<int> cnt, siz, fa;
    std::vector<T> val;

    Splay(int n) {
        tr.resize(n + 1);
        fa.resize(n + 1);
        cnt.resize(n + 1);
        siz.resize(n + 1);
        val.resize(n + 1);
        init();
    }

    void init() {
        std::fill(fa.begin(), fa.end(), 0);
        std::fill(tr.begin(), tr.end(), (std::array<int, 2>){0, 0});
        root = 0;
        tot = 0;
        cnt.clear();
        siz.clear();
        val.clear();
    }

    void clear(int x) {
        tr[x][0] = tr[x][1] = fa[x] = siz[x] = val[x] = cnt[x] = 0;
    }
    void push(int x) {
        siz[x] = siz[tr[x][0]] + siz[tr[x][1]] + cnt[x];
    }
    bool get(int x) {
        return x == tr[fa[x]][1];
    }

    void rotate(int x) {
        int y = fa[x], z = fa[y], type = get(x);
        tr[y][type] = tr[x][type ^ 1];
        if (tr[x][type ^ 1]) {
            fa[tr[x][type ^ 1]] = y;
        }
        tr[x][type ^ 1] = y;
        fa[y] = x, fa[x] = z;
        if (z) {
            tr[z][y == tr[z][1]] = x;
        }
        push(y);
    }

    void splay(int x) {
        for (int f = fa[x]; f = fa[x], f; rotate(x)) {
            if (fa[f]) {
                rotate(get(x) == get(f) ? f : x);
            }
        }
        root = x;
    }

    void insert(T k) {
        if (!root) {
            val[++tot] = k;
            cnt[tot]++;
            root = tot;
            push(root);
            return ;
        }

        int cur = root, f = 0;
        while (1) {
            if (val[cur] == k) {
                cnt[cur]++;
                push(cur), push(f);
                splay(cur);
                break;
            }
            f = cur;
            cur = tr[cur][val[cur] < k];
            if (!cur) {
                val[++tot] = k;
                cnt[tot]++;
                fa[tot] = f;
                tr[f][val[f] < k] = tot;
                push(tot), push(f);
                splay(tot);
                break;
            }
        }
    }

    int rank(int k) {
        int res = 0, cur = root;
        while (1) {
            if (k < val[cur]) {
                cur = tr[cur][0];
            } else {
                res += siz[tr[cur][0]];
                if (k == val[cur]) {
                    splay(cur);
                    return res + 1;
                }
                res += cnt[cur];
                cur = tr[cur][1];
            }
        }
    }

    int kth(int k) {
        int cur = root;
        while (1) {
            if (tr[cur][0] && k <= siz[tr[cur][0]]) {
                cur = tr[cur][0];
            } else {
                k -= cnt[cur] + siz[tr[cur][0]];
                if (k <= 0) {
                    splay(cur);
                    return val[cur];
                }
                cur = tr[cur][1];
            }
        }
    }

    int precursor() {   //前驱
        int cur = tr[root][0];
        if (!cur) {
            return cur;
        }

        while (tr[cur][1]) {
            cur = tr[cur][1];
        }
        splay(cur);
        return cur;
    }
    T precursor(T v) {
        insert(v);
        T res = val[precursor()];
        Delete(v);
        return res;
    }

    int successor() {   //后继
        int cur = tr[root][1];
        if (!cur) {
            return cur;
        }

        while (tr[cur][0]) {
            cur = tr[cur][0];
        }
        splay(cur);
        return cur;
    }
    T successor(T v) {
        insert(v);
        T res = val[successor()];
        Delete(v);
        return res;
    }

    void Delete(int k) {
        rank(k);
        if (cnt[root] > 1) {
            cnt[root]--;
            push(root);
            return ;
        }

        if (!tr[root][0] && !tr[root][1]) {
            clear(root);
            root = 0;
        } else if (!tr[root][0]) {
            int cur = root;
            root = tr[root][1];
            fa[root] = 0;
            clear(cur);
        } else if (!tr[root][1]) {
            int cur = root;
            root = tr[root][0];
            fa[root] = 0;
            clear(cur);
        } else {
            int cur = root;
            int x = precursor();
            fa[tr[cur][1]] = x;
            tr[x][1] = tr[cur][1];
            clear(cur);
            push(root);
        }
    }
};
```


#### Splay区间操作
```c++
template<typename T>
struct Splay {
    int root, tot;
    std::vector<std::array<int, 2>> tr;
    std::vector<int> siz, fa, tag, cnt;
    std::vector<T> val;

    Splay() {}
    Splay(int n) {
        init(n);
    }
    Splay(std::vector<T> _init) {
        init(_init);
    }

    void init(int n) {
        tr.assign(n + 1, (std::array<int, 2>){0, 0});
        fa.assign(n + 1, 0);
        siz.assign(n + 1, 0);
        tag.assign(n + 1, 0);
        cnt.assign(n + 1, 0);
        val.assign(n + 1, T());
        root = tot = 0;
    }

    void init(std::vector<T> _init) {
        int n = _init.size();
        init(n);
        std::function<int(int, int, int)> build = [&](int l, int r, int f) {
            if (l > r) {
                return 0;
            }
            int m = (l + r) / 2, p = ++tot;
            fa[p] = f;
            cnt[p]++;
            val[p] = _init[m - 1];
            siz[p]++;
            tr[p] = {build(l, m - 1, p), build(m + 1, r, p)};
            push(p);
            return p;
        };
        root = build(1, n, 0);
    }

    void push(int p) {
        siz[p] = cnt[p];
        siz[p] += siz[tr[p][0]] + siz[tr[p][1]];
    }
    void pull(int p) {
        if (p && tag[p]) {
            tag[tr[p][1]] ^= 1;
            tag[tr[p][0]] ^= 1;
            std::swap(tr[p][0], tr[p][1]);
            tag[p] = 0;
        }
    }
    bool get(int p) {
        return p == tr[fa[p]][1];
    }

    void rotate(int x) {
        int y = fa[x], z = fa[y], type = get(x);
        tr[y][type] = tr[x][type ^ 1];
        if (tr[x][type ^ 1]) {
            fa[tr[x][type ^ 1]] = y;
        }
        tr[x][type ^ 1] = y;
        fa[y] = x, fa[x] = z;
        if (z) {
            tr[z][y == tr[z][1]] = x;
        }
        push(y);
    }

    void splay(int x, int g) {
        for (int f; (f = fa[x]) != g; rotate(x)) {
            if (fa[f] != g) {
                rotate(get(x) == get(f) ? f : x);
            }
        }
        if(g == 0) {
            root = x;
        }
    }

    int find(int x) {
        int cur = root;
        while (1) {
            pull(cur);
            if (tr[cur][0] && x <= siz[tr[cur][0]]) {
                cur = tr[cur][0];
            } else {
                x -= siz[tr[cur][0]] + 1;
                if (x <= 0) {
                    return cur;
                }
                cur = tr[cur][1];
            }
        }
    }

    void reverse(int l, int r) {    //区间翻转
        r += 2;
        l = find(l), r = find(r);

        splay(l, 0);
        splay(r, l);

        int cur = tr[root][1];
        cur = tr[cur][0];
        tag[cur] ^= 1;
    }

    void print(int p) { //输出中序遍历结果
        pull(p);
        if (tr[p][0]) {
            print(tr[p][0]);
        }
        if (val[p] != -inf && val[p] != inf) {
            cout << val[p] << " ";
        }
        if (tr[p][1]) {
            print(tr[p][1]);
        }
    }
};
```

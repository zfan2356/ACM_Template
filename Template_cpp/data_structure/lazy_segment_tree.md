## 带懒标记线段树

```c++
template<class Info, class Tag>
struct LazySegmentTree {
    int n;
    std::vector<Info> info;
    std::vector<Tag> tag;
    LazySegmentTree() : n(0) {}
    LazySegmentTree(int n_, Info v_ = Info()) {
        init(n_, v_);
    }
    template<class T>
    LazySegmentTree(std::vector<T> init_) {
        init(init_);
    }
    void init(int n_, Info v_ = Info()) {
        init(std::vector(n_, v_));
    }
  
    template<class T>
    void init(std::vector<T> init_) {
        n = init_.size();
        info.assign(4 << std::__lg(n), Info());
        tag.assign(4 << std::__lg(n), Tag());
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                info[p] = init_[l];
                return;
            }
            int m = (l + r) / 2;
            build(2 * p, l, m);
            build(2 * p + 1, m, r);
            pull(p);
        };
        build(1, 0, n);
    }
    void pull(int p) {
        info[p] = info[2 * p] + info[2 * p + 1];
    }
    void apply(int p, const Tag &v) {
        info[p].apply(v);
        tag[p].apply(v);
    }
    void push(int p) {
        apply(2 * p, tag[p]);
        apply(2 * p + 1, tag[p]);
        tag[p] = Tag();
    }
  
    void modify(int p, int l, int r, int x, const Info &v) {
        if (r - l == 1) {
            info[p] = v;
            return;
        }
        int m = (l + r) / 2;
        push(p);
        if (x < m) {
            modify(2 * p, l, m, x, v);
        } else {
            modify(2 * p + 1, m, r, x, v);
        }
        pull(p);
    }
    void modify(int p, const Info &v) {
        modify(1, 0, n, p, v);
    }
  
    Info rangeQuery(int p, int l, int r, int x, int y) {
        if (l >= y || r <= x) {
            return Info();
        }
        if (l >= x && r <= y) {
            return info[p];
        }
        int m = (l + r) / 2;
        push(p);
        return rangeQuery(2 * p, l, m, x, y) + rangeQuery(2 * p + 1, m, r, x, y);
    }
    Info rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }
  
    void rangeApply(int p, int l, int r, int x, int y, const Tag &v) {
        if (l >= y || r <= x) {
            return;
        }
        if (l >= x && r <= y) {
            apply(p, v);
            return;
        }
        int m = (l + r) / 2;
        push(p);
        rangeApply(2 * p, l, m, x, y, v);
        rangeApply(2 * p + 1, m, r, x, y, v);
        pull(p);
    }
    void rangeApply(int l, int r, const Tag &v) {
        return rangeApply(1, 0, n, l, r, v);
    }
  
    template<class F>
    int findFirst(int p, int l, int r, int x, int y, F pred) {
        if (l >= y || r <= x || !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        push(p);
        int res = findFirst(2 * p, l, m, x, y, pred);
        if (res == -1) {
            res = findFirst(2 * p + 1, m, r, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findFirst(int l, int r, F pred) {
        return findFirst(1, 0, n, l, r, pred);
    }
  
    template<class F>
    int findLast(int p, int l, int r, int x, int y, F pred) {
        if (l >= y || r <= x || !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        push(p);
        int res = findLast(2 * p + 1, m, r, x, y, pred);
        if (res == -1) {
            res = findLast(2 * p, l, m, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findLast(int l, int r, F pred) {
        return findLast(1, 0, n, l, r, pred);
    }
};
 
struct Tag {
    long long add = 0;
    void apply(const Tag &t) & {
        add += t.add;
    }
};
 
struct Info {
    long long mn = 1E18;
    void apply(const Tag &t) & {
        mn += t.add;
    }
};
 
Info operator+(const Info &a, const Info &b) {
    return {std::min(a.mn, b.mn)};
}

```

#### 技巧: 区间历史最值

```c++
struct Tag {
    i64 t1 = 0, t2 = 0;
    i64 s1 = INT_MIN, s2 = INT_MIN;

    void apply_t(const Tag &t) & {
        if (s1 == INT_MIN) {
            t2 = max(t2, t1 + t.t2);
            t1 += t.t1;
        } else {
            s2 = max(s2, s1 + t.t2);
            s1 += t.t1;
        }
    }

    void apply_s(const Tag &t) & {
        s2 = max(s2, t.s2);
        s1 = t.s1;
    }

    void apply(const Tag &t) & {
        apply_t(t);
        if (t.s1 != INT_MIN || t.s2 != INT_MIN) {
            apply_s(t);
        }
    }
};


struct Info {
    i64 a1 = 0;  //区间最大值
    i64 a2 = 0;  //区间历史最大值
    void apply(const Tag &t) & {
        a2 = max(a2, a1 + t.t2);
        a1 += t.t1;
        if (t.s1 != INT_MIN || t.s2 != INT_MIN) {
            a2 = max(a2, t.s2);
            a1 = t.s1;
        }
    }
};

Info operator+(const Info &a, const Info &b) {
    return {std::max(a.a1, b.a1), std::max(a.a2, b.a2)};
}

for (int i = 0; i < m; i++) {
    int op, l, r;
    i64 x;
    cin >> op >> l >> r;
    l--;
    if (op == 1) {  //区间最大值
        cout << lazySegmentTree.rangeQuery(l, r).a1 << endl;
    } else if (op == 2){    //区间历史最大值
        cout << lazySegmentTree.rangeQuery(l, r).a2 << endl;
    } else if (op == 3) {   //区间加
        cin >> x;
        lazySegmentTree.rangeApply(l, r, {x, x, INT_MIN, INT_MIN});
    } else {        //区间赋值
        cin >> x;
        lazySegmentTree.rangeApply(l, r, {0, 0, x, x});
    }
}
```

## 带懒标记动态开点权值线段树

```c++
struct Node {
    i64 v = inf, add = 0;
    Node *l = nullptr, *r = nullptr;

    void apply(int x) { // 区间赋值
        v = x;
        add = x;
    }
};

void pull(Node *&t) {  //区间求min
    t->v = std::min(t->l->v, t->r->v);
}

void make(Node *&p, int l, int r) { // 创建新节点时为其赋值, query是[l, r] 的值
    if (p == nullptr) {
        p = new Node();
        p->v = Query(l, r);
    }
}

void Apply(Node *&t, int l, int r) {
    if (t->add) {
        int m = (l + r) / 2;
        make(t->l, l, m);
        t->l->apply(t->add);
      
        make(t->r, m, r);
        t->r->apply(t->add);
        t->add = 0;
    }
}

void rangeApply(Node *&p, int l, int r, int x, int y, int v) {
    make(p, l, r);
    if (l >= y || r <= x) {
        return ;
    }
    if (x <= l && r <= y) {
        p->apply(v);
        return ;
    }
    Apply(p, l, r);
    int m = (l + r) / 2;
    rangeApply(p->l, l, m, x, y, v);
    rangeApply(p->r, m, r, x, y, v);
    pull(p);
}

int rangeQuery(Node *&p, int l, int r, int x, int y) {
    if (l >= y || r <= x) {
        return inf;
    }
    if (p == nullptr) {
        make(p, l, r);
    }
    if (x <= l && r <= y) {
        return p->v;
    }
    int m = (l + r) / 2;
    Apply(p, l, r);
    return std::min(rangeQuery(p->l, l, m, x, y), rangeQuery(p->r, m, r, x, y));
}

```

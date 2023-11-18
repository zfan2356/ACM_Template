```c++
std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
constexpr int inf = 1e9;

template<typename T>
struct Tree {
    struct Node {
        Node *l;
        Node *r;
        int siz, cnt, tag;
        T val;
        u64 w;

        Node(T v = 0) : l{nullptr}, r{nullptr}, siz{1}, cnt{1}, val{v}, w{rng()} {}
    };
    Node* root = nullptr;

    void pull(Node* t) {
        if (t != nullptr && t->tag != 0) {
            std::swap(t->l, t->r);
            if (t->l != nullptr) {
                t->l->tag ^= 1;
            }
            if (t->r != nullptr) {
                t->r->tag ^= 1;
            }
            t->tag = 0;
        }
    }

    void push(Node* t) {
        t->siz = t->cnt;
        if (t->l != nullptr) {
            t->siz += t->l->siz;
        }
        if (t->r != nullptr) {
            t->siz += t->r->siz;
        }
    }

    std::pair<Node*, Node*> split(Node* t, T v) {   //按值分裂
        if (t == nullptr) return {nullptr, nullptr};
        pull(t);
        if (t->val <= v) {
            auto [l, r] = split(t->r, v);
            t->r = l;
            push(t);
            return {t, r};
        } else {
            auto [l, r] = split(t->l, v);
            t->l = r;
            push(t);
            return {l, t};
        }
    }

    std::array<Node*, 3> split_by_rank(Node* t, int rk) {   //按排名分裂
        if (t == nullptr) return {nullptr, nullptr, nullptr};
        int siz_l = t->l == nullptr ? 0 : t->l->siz;
        pull(t);
        if (rk <= siz_l) {
            auto [l, m, r] = split_by_rank(t->l, rk);
            t->l = r;
            push(t);
            return {l, m, t};
        } else if (rk <= siz_l + t->cnt) {
            Node* l = t->l;
            Node* r = t->r;
            t->l = t->r = nullptr;
            return {l, t, r};
        } else {
            auto [l, m, r] = split_by_rank(t->r, rk - siz_l - t->cnt);
            t->r = l;
            push(t);
            return {t, m, r};
        }
    }

    std::pair<Node*, Node*> split_by_siz(Node* t, int siz) {
        if (t == nullptr) {
            return {nullptr, nullptr};
        }
        pull(t);
        int siz_l = t->l == nullptr ? 0 : t->l->siz;
        if (siz <= siz_l) {
            auto [l, r] = split_by_siz(t->l, siz);
            t->l = r;
            push(t);
            return {l, t};
        } else {
            auto [l, r] = split_by_siz(t->r, siz - siz_l - 1);
            t->r = l;
            push(t);
            return {t, r};
        }
    }

    Node* merge(Node* u, Node* v) { // u中所有val < v中所有val
       if (u == nullptr && v == nullptr) {
           return nullptr;
       } else if (u == nullptr) {
           return v;
       } else if (v == nullptr) {
           return u;
       } else {
           pull(u), pull(v);
           if (u->w < v->w) {
               u->r = merge(u->r, v);
               push(u);
               return u;
           } else {
               v->l = merge(u, v->l);
               push(v);
               return v;
           }
       }
    }

    void reverse(int l, int r) {
        auto L = split_by_siz(root, l - 1);
        auto M = split_by_siz(L.second, r - l + 1);
        M.first->tag ^= 1;
        root = merge(L.first, merge(M.first, M.second));
    }

    void insert(T v) {
        if (root == nullptr) {
            root = new Node(v);
            return;
        }
        auto [l, r] = split(root, v);
        auto ltr = split(l, v - 1);
        Node* node;
        if (ltr.second == nullptr) {
            node = new Node(v);
        } else {
            ltr.second->cnt++;
            push(ltr.second);
        }

        Node* tree = merge(ltr.first, ltr.second == nullptr ? node : ltr.second);
        root = merge(tree, r);
    }

    void delete_(T v){
        auto [l, r] = split(root, v);
        auto ltr = split(l, v - 1);

        if (ltr.second->cnt > 1) {
            ltr.second->cnt--;
            push(ltr.second);
            ltr.first = merge(ltr.first, ltr.second);
        } else {
            if (l == ltr.second) {
                l = nullptr;
            }
            delete ltr.second;
            ltr.second = nullptr;
        }
        root = merge(ltr.first, r);
    }

    int rank(Node* t, T v) {
        auto [l, r] = split(t, v - 1);
        int res = (l == nullptr ? 0 : l->siz) + 1;
        root = merge(l, r);
        return res;
    }

    T kth(Node* t, int rk) {
        auto [l, m, r] = split_by_rank(t, rk);
        T res = m->val;
        root = merge(merge(l, m), r);
        return res;
    }

    T precursor(T v) {
        auto [l, r] = split(root, v - 1);
        T res = kth(l, l->siz);
        root = merge(l, r);
        return res;
    }

    T successor(T v) {
        auto [l, r] = split(root, v);
        T res = kth(r, 1);
        root = merge(l, r);
        return res;
    }

    void print(Node *t) {
        if (t == nullptr) {
            return;
        }
        pull(t);
        print(t->l);
        std::cout << t->val << " ";
        print(t->r);
    }
};
```
template <typename T>
struct Fenwick {
    int n;
    std::vector<T> a;

    Fenwick(int n = 0) {
        init(n);
    }

    Fenwick(std::vector<T>& b) {
        init(b);
    }

    void init(int n) {
        this->n = n;
        a.assign(n, T());
    }

    void init(std::vector<T> &b) {
        this->n = b.size();
        a.assign(n, T());
        for (int i = 1; i <= n; i++) {
            a[i - 1] += b[i - 1];
            int j = i + (i & -i);
            if (j <= n) {
                a[j - 1] += a[i - 1];
            }
        }
    }

    void add(int x, T v) {
        for (int i = x + 1; i <= n; i += i & -i) {
            a[i - 1] += v;
        }
    }

    T sum(int x) {
        auto ans = T();
        for (int i = x; i > 0; i -= i & -i) {
            ans += a[i - 1];
        }
        return ans;
    }

    T rangeSum(int l, int r) {
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


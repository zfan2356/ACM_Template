#### MTT模板

```c++
namespace P1 {
    typedef int DAT;
    constexpr DAT P = 469762049, 998244353, 1004535809;//三个模数, 复制粘贴即可
    template<typename T>
    constexpr T Pow(T a, DAT b = P - 2) {
        T res = 1;
        for (; b; b /= 2, a *= a)if (b % 2)res *= a;
        return res;
    }
    struct Z {
        DAT z;
        constexpr DAT norm(DAT x) {
            if (x < 0)x += P;
            if (x >= P)x -= P;
            return x;
        }
        constexpr Z(DAT x = 0) : z(norm(x)) {}
        constexpr DAT val() const { return z; }
        constexpr Z operator-() const { return Z(P - z); }
        constexpr Z inv() const {/*assert(z!=0);*/return Pow(*this); }
        constexpr Z &operator*=(const Z &r) {
            z = (i64) z * r.z % P;
            return *this;
        }// int
        //Z&operator*=(const Z&r){u64 res=(u64)z*r.z-(u64)((long double)z/P*r.z+0.5L)*P;z=(res<P?res:res+P);return *this;}// long long
        constexpr Z &operator+=(const Z &r) {
            z = norm(z + r.z);
            return *this;
        }
        constexpr Z &operator-=(const Z &r) {
            z = norm(z - r.z);
            return *this;
        }
        constexpr Z &operator/=(const Z &r) { return (*this) *= r.inv(); }
        constexpr friend Z operator*(const Z &l, const Z &r) {
            Z res = l;
            return res *= r;
        }
        constexpr friend Z operator+(const Z &l, const Z &r) {
            Z res = l;
            return res += r;
        }
        constexpr friend Z operator-(const Z &l, const Z &r) {
            Z res = l;
            return res -= r;
        }
        constexpr friend Z operator/(const Z &l, const Z &r) {
            Z res = l;
            return res /= r;
        }
        friend constexpr bool operator==(Z lhs, Z rhs) {
            return lhs.val() == rhs.val();
        }
        friend constexpr bool operator!=(Z lhs, Z rhs) {
            return lhs.val() != rhs.val();
        }
    };

    constexpr Z G(3), INVG = Pow(G);// generator
    typedef std::vector <Z> Poly;
    int EX2(int n) { return 1 << (32 - __builtin_clz(n - 1)); }
    std::vector<int> Rev{0};
    Z Invlen(1);

    void fft(Poly &a, int type = 1) {// a.size == 2^k   type==1 -> dft   type==-1 -> idft
        int n = a.size();
        if (n != Rev.size()) {
            Rev.resize(n);
            int k = 1 << (__builtin_ctz(n) - 1);
            for (int i = 1; i < n; ++i)Rev[i] = (Rev[i / 2] / 2) | ((i & 1) * k);
            Invlen = Z(n).inv();
        }
        for (int i = 1; i < n; ++i) {
            if (i < Rev[i]) std::swap(a[i], a[Rev[i]]);
        }
        for (int i = 1; i < n; i *= 2) {
            Z wn = Pow((type == 1 ? G : INVG), ((P - 1) / (2 * i))), w(1), t;
            for (int j = 0; j < n; j += 2 * i, w = 1)
                for (int k = j; k < i + j; ++k, w *= wn)
                    t = a[k + i] * w, a[k + i] = a[k] - t, a[k] += t;
        }
        if (type < 0)for (int i = 0; i < n; ++i)a[i] *= Invlen;
    }

    void dot(Poly &a, const Poly &b) {// a.size == b.size
        for (int i = 0; i < int(a.size()); ++i)a[i] *= b[i];
    }

    Poly operator*(Poly a, Poly b) {//for all size
        if (a.size() == 0 || b.size() == 0)return Poly();
        int n = a.size() + b.size() - 1;
        if (std::min(a.size(), b.size()) <= 8) {
            Poly c(n);
            for (int i = 0; i < a.size(); ++i)
                for (int j = 0; j < b.size(); ++j)
                    c[i + j] += a[i] * b[j];
            return c;
        }
        int m = EX2(n);
        a.resize(m); fft(a);
        b.resize(m); fft(b);
        dot(a, b);
        fft(a, -1); a.resize(n);
        return a;
    }
}
namespace MTT {
    typedef int DAT;
    constexpr DAT P = 1e9 + 7;// 注意赋值
    template<typename T>
    constexpr T Pow(T a, DAT b = P - 2) {
        T res = 1;
        for (; b; b /= 2, a *= a)if (b % 2)res *= a;
        return res;
    }//必须是封装后的类型
    struct Z {
        DAT z;
        constexpr DAT norm(DAT x) {
            if (x < 0)x += P;
            if (x >= P)x -= P;
            return x;
        } // -P <= x < 2P
        constexpr Z(DAT x = 0) : z(norm(x)) {}
        constexpr DAT val() const { return z; }
        constexpr Z operator-() const { return Z(P - z); }
        constexpr Z inv() const {/*assert(z!=0);*/return Pow(*this); }
        constexpr Z &operator*=(const Z &r) {
            z = (i64) z * r.z % P;
            return *this;
        }// int
        //Z&operator*=(const Z&r){u64 res=(u64)z*r.z-(u64)((long double)z/P*r.z+0.5L)*P;z=(res<P?res:res+P);return *this;}// long long
        constexpr Z &operator+=(const Z &r) {
            z = norm(z + r.z);
            return *this;
        }
        constexpr Z &operator-=(const Z &r) {
            z = norm(z - r.z);
            return *this;
        }
        constexpr Z &operator/=(const Z &r) { return (*this) *= r.inv(); }
        constexpr friend Z operator*(const Z &l, const Z &r) {
            Z res = l;
            return res *= r;
        }
        constexpr friend Z operator+(const Z &l, const Z &r) {
            Z res = l;
            return res += r;
        }
        constexpr friend Z operator-(const Z &l, const Z &r) {
            Z res = l;
            return res -= r;
        }
        constexpr friend Z operator/(const Z &l, const Z &r) {
            Z res = l;
            return res /= r;
        }
        friend constexpr bool operator==(Z lhs, Z rhs) {
            return lhs.val() == rhs.val();
        }
        friend constexpr bool operator!=(Z lhs, Z rhs) {
            return lhs.val() != rhs.val();
        }
    };

    typedef std::vector <Z> Poly;

    int EX2(int n) { return 1 << (32 - __builtin_clz(n - 1)); }

    void dot(Poly &a, const Poly &b) {// a.size == b.size
        for (int i = 0; i < int(a.size()); ++i)a[i] *= b[i];
    }

    template<typename T>
    constexpr void exgcd(T a, T b, T &x, T &y) {
        if (b == 0) { x = 1, y = 0; }
        else {
            exgcd(b, a % b, y, x);
            y -= (a / b) * x;
        }
    }

    Poly operator*(Poly a, Poly b) {//for all size
        if (a.size() == 0 || b.size() == 0)return Poly();
        int n = a.size() + b.size() - 1;
        if (std::min(a.size(), b.size()) <= 12) {
            Poly c(n);
            for (int i = 0; i < a.size(); ++i)
                for (int j = 0; j < b.size(); ++j)
                    c[i + j] += a[i] * b[j];
            return c;
        }

        P1::Poly a1(n);
        for (int i = 0; i < a.size(); ++i)a1[i] = a[i].val();
        P1::Poly b1(n);
        for (int i = 0; i < b.size(); ++i)b1[i] = b[i].val();
        a1 = a1 * b1;

        P2::Poly a2(n);
        for (int i = 0; i < a.size(); ++i)a2[i] = a[i].val();
        P2::Poly b2(n);
        for (int i = 0; i < b.size(); ++i)b2[i] = b[i].val();
        a2 = a2 * b2;

        P3::Poly a3(n);
        for (int i = 0; i < a.size(); ++i)a3[i] = a[i].val();
        P3::Poly b3(n);
        for (int i = 0; i < b.size(); ++i)b3[i] = b[i].val();
        a3 = a3 * b3;

        Poly res(n);

        typedef __int128 LLL;
        static constexpr LLL p = (LLL) P1::P * (LLL) P2::P * (LLL) P3::P;
        static constexpr LLL t1 =
                LLL(P1::Pow(P1::Z(P2::P % P1::P)).val()) * LLL(P1::Pow(P1::Z(P3::P % P1::P)).val()) * P2::P * P3::P % p;
        static constexpr LLL t2 =
                LLL(P2::Pow(P2::Z(P1::P % P2::P)).val()) * LLL(P2::Pow(P2::Z(P3::P % P2::P)).val()) * P1::P * P3::P % p;
        static constexpr LLL t3 =
                LLL(P3::Pow(P3::Z(P1::P % P3::P)).val()) * LLL(P3::Pow(P3::Z(P2::P % P3::P)).val()) * P1::P * P2::P % p;
        for (int i = 0; i < n; ++i)res[i] = (t1 * a1[i].val() + t2 * a2[i].val() + t3 * a3[i].val()) % p % P;
        return res;
    }

    Poly inv2k(const Poly &f) {// f.size == 2^k
        int n = f.size();
        //assert(n==EX2(n));///
        //assert(f[0].val()!=0);///
        Poly g{f[0].inv()};
        g.reserve(n);
        for (int k = 2; k <= n; k *= 2) {
            Poly ff(f.begin(), f.begin() + k);
            ff = (g * g) * ff;
            g.resize(k);
            for (int i = 0; i < k / 2; ++i)g[i] = 2 * g[i] - ff[i];
            for (int i = k / 2; i < k; ++i)g[i] = -ff[i];
        }
        return g;
    }

    Poly inv(Poly f) {// for all size
        int n = f.size(), m = EX2(n);
        f.resize(m);
        f = inv2k(f);
        f.resize(n);
        return f;
    }

    Poly dx(const Poly &f) {// for all size 求导
        if (f.size() == 0)return Poly();
        Poly res(f.size() - 1);
        for (int i = 0; i < res.size(); ++i)res[i] = f[i + 1] * (i + 1);
        return res;
    }

    Poly Inv({Z(1)});

    Poly fdx(const Poly &f) {// for all size    积分
        Poly res(f.size() + 1);
        if (f.size() > Inv.size())
            for (int i = Inv.size() + 1; i <= f.size(); ++i)
                Inv.push_back((P - P / i) * Inv[P % i - 1]);
        for (int i = 0; i < f.size(); ++i)res[i + 1] = f[i] * Inv[i];
        return res;
    }

    Poly ln2k(const Poly &f) {// f.size == 2^k
        //assert(f[0].val()==1);///
        Poly res = dx(f) * inv2k(f);
        res.resize(f.size() - 1);
        return fdx(res);
    }

    Poly exp2k(const Poly &f) {// f.size == 2^k
        int n = f.size();
        Poly g({Z(1)});
        g.reserve(n);
        for (int k = 2; k <= n; k *= 2) {
            //assert(g[0].val()==1);///
            g.resize(k);
            Poly h = ln2k(g);
            h[0] = 1;
            for (int i = 1; i < k; ++i)h[i] = f[i] - h[i];
            g = g * h;
            //assert(g[0].val()==1);///
            g.resize(k);
            //assert(g[0].val()==1);///
        }
        return g;
    }

    Poly sqrt2k(const Poly &f, Z q = Z(1)) {// f.size == 2^k
        int n = f.size();
        Poly g{q};
        for (int k = 2; k <= n; k *= 2) {
            g.resize(k);
            Poly h = Poly(f.begin(), f.begin() + k) * inv2k(g);
            for (int i = 0; i < k; ++i)g[i] = (g[i] + h[i]) * ((P + 1) / 2);
        }
        return g;
    }

    std::pair <Poly, Poly> mod(const Poly &a, const Poly &b) {// for all size
        const int n = a.size() - 1, m = b.size() - 1, t = n - m + 1;
        assert(t > 0);
        Poly f(t);
        for (int i = 0; i < t; ++i)f[i] = a[n - i];
        Poly g(t);
        for (int i = 0; i < std::min(t, m + 1); ++i)g[i] = b[m - i];
        g = inv(g) * f;
        Poly q(t);
        for (int i = 0; i < t; ++i)q[i] = g[t - i - 1];
        f = b * q;
        Poly r(m);
        for (int i = 0; i < m; ++i)r[i] = a[i] - f[i];
        return {q, r};
    }
}
```
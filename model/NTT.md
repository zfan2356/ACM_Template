```c++
#define fp(i, a, b) for (int i = (a), i##_ = (b) + 1; i < i##_; ++i)
#define fd(i, a, b) for (int i = (a), i##_ = (b) - 1; i > i##_; --i)
#define MUL(a, b) (i64(a) * (b) % P)
#define ADD(a, b) (((a) += (b)) >= P ? (a) -= P : 0) // (a += b) %= P
#define SUB(a, b) (((a) -= (b)) < 0 ? (a) += P: 0)  // ((a -= b) += P) %= P

/*
template<typename A, typename B>
inline std::ostream &operator<<(std::ostream &out, const std::pair <A, B> &p) {
    return out << "(" << p.first << ", " << p.second << ")";
}

template<typename T>
inline std::ostream &operator<<(std::ostream &out, const std::vector <T> &a) {
    out << "[";
    for (int i = 0; i < a.size(); i++) {
        if (i) out << ',';
        out << ' ' << a[i];
    }
    return out << " ]";
}

template<typename T>
inline std::ostream &operator<<(std::ostream &out, const std::set <T> &a) {
    return out << std::vector<T>(all(a));
}
 */

const int P = 998244353;
const int N = 3e5 + 10;
using Poly = std::vector<int>;
using MultiPoly = std::vector<Poly>;

//快速幂
int qpow(i64 a, int b = P - 2, i64 x = 1) {
    for (; b; b >>= 1, a = a * a % P) {
        if (b & 1) {
            x = x * a % P;
        }
    }
    return x;
}


class Cipolla { //二次剩余
    int P, I2{};
    using pll = std::pair<i64, i64>;
#define X first
#define Y second
    i64 mul(i64 a, i64 b) const { return a * b % P; }

    pll mul(pll a, pll b) const {
        return { (a.X * b.X + I2 * a.Y % P * b.Y) % P, (a.X * b.Y + a.Y * b.X) % P };
    }

    template<class T> T POW(T a, int b, T x) {
        for (; b; b >>= 1, a = mul(a, a)) {
            if (b & 1) x = mul(x, a);
        }
        return x;
    }
public:
    Cipolla(int p = 0) : P(p) {}
    std::pair<int, int> sqrt(int n) {
        int a = rand(), x;
        if (!(n %= P)) return { 0, 0};
        if (POW(n, (P - 1) >> 1, 1) == P - 1) return { -1, -1};

        while (POW(I2 = ((i64)a * a - n + P) % P, (P - 1) >> 1, 1) == 1) {
            a = rand();
        }
        x = (int)POW(pll{ a, 1 }, (P + 1) >> 1, { 1, 0 }).X;

        if (2 * x > P) {
            x = P - x;
        }
        return { x, P - x };
    }
#undef X
#undef Y
};

//预处理L以内的逆元(0 ~ L-1)
Poly getInv(int L) {
    Poly inv(L);
    inv[1] = 1;
    fp(i, 2, L - 1) {
        inv[i] = MUL((P - P / i), inv[P % i]);
    }
    return inv;
}

auto inv = getInv(N);

namespace NTT {
    const int g = 3;
    Poly Omega(int L) {
        int wn = qpow(g, P / L);
        Poly w(L); w[L >> 1] = 1;
        fp(i, L / 2 + 1, L - 1) {
            w[i] = MUL(w[i - 1], wn);
        }
        fd(i, L / 2 - 1, 1) {
            w[i] = w[i << 1];
        }
        return w;
    }
    auto W = Omega(1 << 23); // 注意这边的size，如果大于3e5，改成23；
    void DIF(int* a, int n) {
        for (int k = n >> 1; k; k >>= 1)
            for (int i = 0, y; i < n; i += k << 1)
                for (int j = 0; j < k; ++j) {
                    y = a[i + j + k];
                    a[i + j + k] = MUL(a[i + j] - y + P, W[k + j]);
                    ADD(a[i + j], y);
                }
    }
    void IDIT(int* a, int n) {
        for (int k = 1; k < n; k <<= 1)
            for (int i = 0, x, y; i < n; i += k << 1)
                for (int j = 0; j < k; ++j) {
                    x = a[i + j];
                    y = MUL(a[i + j + k], W[k + j]);
                    a[i + j + k] = x - y < 0 ? x - y + P : x - y, ADD(a[i + j], y);
                }
        int Inv = P - (P - 1) / n;
        fp(i, 0, n - 1) {
            a[i] = MUL(a[i], Inv);
        }
        std::reverse(a + 1, a + n);
    }
}
namespace FWT {
    void FWTor(Poly& a, bool rev) {
        int n = a.size();
        for (int l = 2, m = 1; l <= n; l <<= 1, m <<= 1)
            for (int j = 0; j < n; j += l) fp(i, 0, m - 1) {
                    if (!rev) {
                        a[i + j + m] = ADD(a[i + j + m], a[i + j]);
                    } else {
                        a[i + j + m] = SUB(a[i + j + m], a[i + j]);
                    }
                }
    }
    void FWTand(Poly& a, bool rev) {
        int n = a.size();
        for (int l = 2, m = 1; l <= n; l <<= 1, m <<= 1)
            for (int j = 0; j < n; j += l) fp(i, 0, m - 1) {
                    if (!rev) {
                        a[i + j] = ADD(a[i + j], a[i + j + m]);
                    } else {
                        a[i + j] = SUB(a[i + j], a[i + j + m]);
                    }
                }
    }
    void FWTxor(Poly& a, bool rev) {
        int n = a.size(), inv2 = (P + 1) >> 1;
        for (int l = 2, m = 1; l <= n; l <<= 1, m <<= 1)
            for (int j = 0; j < n; j += l) fp(i, 0, m - 1) {
                    int x = a[i + j], y = a[i + j + m];
                    if (!rev) {
                        a[i + j] = ADD(x, y);
                        a[i + j + m] = SUB(x, y);
                    }
                    else {
                        a[i + j] = MUL(ADD(x, y), inv2);
                        a[i + j + m] = MUL(SUB(x, y), inv2);
                    }
                }
    }
}

namespace Polynomial {
    // size确定以及NTT乘法
    int norm(int n) {
        return 1 << ((int)log2(n - 1) + 1);
    }
    void norm(Poly& a) {
        if (!a.empty()) {
            a.resize(norm(a.size()), 0);
        } else {
            a = {0};
        }
    }
    void DFT(Poly& a) {
        NTT::DIF(a.data(), a.size());
    }

    void IDFT(Poly& a) {
        NTT::IDIT(a.data(), a.size());
    }

    Poly& dot(Poly& a, Poly& b) {
        fp(i, 0, a.size() - 1) {
            a[i] = MUL(a[i], b[i]);
        }
        return a;
    }

    // 和整数的乘除运算
    Poly& operator*=(Poly& a, int b) {
        for (auto& x : a) {
            x = MUL(x, b);
        }
        return a;
    }
    Poly operator*(Poly a, int b) {
        return a *= b;
    }
    Poly operator*(int a, Poly b) {
        return b * a;
    }
    Poly& operator/=(Poly& a, int b) {
        return a *= qpow(b);
    }
    Poly operator/(Poly a, int b) {
        return a /= b;
    }

    // 多项式之间的加减运算
    Poly& operator+=(Poly& a, Poly b) {
        a.resize(std::max(a.size(), b.size()));
        fp(i, 0, b.size() - 1) {
            ADD(a[i], b[i]);
        }
        return a;
    }
    Poly operator+(Poly a, Poly b) {
        return a += b;
    }
    Poly& operator-=(Poly& a, Poly b) {
        a.resize(std::max(a.size(), b.size()));
        fp(i, 0, b.size() - 1) {
            SUB(a[i], b[i]);
        }
        return a;
    }
    Poly operator-(Poly a, Poly b) {
        return a -= b;
    }

    // 多项式乘法
    Poly operator*(Poly a, Poly b) {
        int n = a.size() + b.size() - 1, L = norm(n);
        if (a.size() <= 30 || b.size() <= 30) {
            Poly c(n);
            fp(i, 0, a.size() - 1) {
                fp(j, 0, b.size() - 1) {
                    c[i + j] = (c[i + j] + (i64)a[i] * b[j]) % P;
                }
            }
            return c;
        }
        a.resize(L), b.resize(L);
        DFT(a), DFT(b), dot(a, b), IDFT(a);
        a.resize(n);
        return a;
    }

    // 多项式逆元
    Poly Inv2k(Poly a) { // |a| = 2 ^ k
        int n = a.size(), m = n >> 1;
        if (n == 1) {
            return { qpow(a[0]) };
        }
        Poly b = Inv2k(Poly(a.begin(), a.begin() + m)), c = b;
        b.resize(n);
        DFT(a), DFT(b);
        dot(a, b);
        IDFT(a);
        fp(i, 0, n - 1) {
            a[i] = i < m ? 0 : P - a[i];
        }
        DFT(a);
        dot(a, b);
        IDFT(a);
        move(c.begin(), c.end(), a.begin());
        return a;
    }

    Poly Inv(Poly a) {
        int n = a.size();
        norm(a);
        a = Inv2k(a);
        a.resize(n);
        return a;
    }

    // 多项式除法/取模
    Poly operator/(Poly a, Poly b) {
        int k = a.size() - b.size() + 1;
        if (k < 0) {
            return { 0 };
        }
        reverse(a.begin(), a.end());
        reverse(b.begin(), b.end());
        b.resize(k);
        a = a * Inv(b);
        a.resize(k);
        reverse(a.begin(), a.end());
        return a;
    }
    std::pair<Poly, Poly> operator%(Poly a, const Poly& b) {
        Poly c = a / b;
        a -= b * c;
        a.resize(b.size() - 1);
        return { c, a };
    }

    // 多项式求导积分
    Poly deriv(Poly a) {
        fp(i, 1, a.size() - 1) {
            a[i - 1] = MUL(i, a[i]);
        }
        return a.pop_back(), a;
    }
    Poly integ(Poly a) {
        a.push_back(0);
        fd(i, a.size() - 1, 1) {
            a[i] = MUL(inv[i], a[i - 1]);
        }
        return a[0] = 0, a;
    }

    // 取ln
    Poly Ln(Poly a) {
        int n = a.size();
        a = deriv(a) * Inv(a);
        a.resize(n - 1);
        return integ(a);
    }

    // 取exp
    Poly Exp(Poly a) {
        int n = a.size(), k = norm(n);
        Poly b = { 1 }, c, d;
        a.resize(k);
        for (int L = 2; L <= k; L <<= 1) {
            d = b;
            b.resize(L);
            c = Ln(b);
            c.resize(L);
            fp(i, 0, L - 1) {
                c[i] = a[i] - c[i] + (a[i] < c[i] ? P : 0);
            }
            ADD(c[0], 1);
            DFT(b), DFT(c);
            dot(b, c);
            IDFT(b);
            move(d.begin(), d.end(), b.begin());
        }
        b.resize(n);
        return b;
    }

    // 开根
    Poly Sqrt(Poly a) {
        int n = a.size(), k = norm(n);
        a.resize(k);
        Poly b = { (new Cipolla(P))->sqrt(a[0]).first, 0 }, c;
        for (int L = 2; L <= k; L <<= 1) {
            b.resize(L);
            c = Poly(a.begin(), a.begin() + L) * Inv2k(b);
            fp(i, L / 2, L - 1) {
                b[i] = MUL(c[i], (P + 1) / 2);
            }
        }
        b.resize(n);
        return b;
    }

    // 多项式快速幂
    Poly Pow1(Poly& a, int b) {  // a[0] = 1, 循环卷积
        return Exp(Ln(a) * b);
    }
    Poly Pow2(Poly& a, int b) {
        int n = (a.size() - 1) * b + 1, L = norm(n);
        a.resize(L);
        DFT(a);
        fp(i, 0, L - 1) {
            a[i] = qpow(a[i], b);
        }
        IDFT(a);
        return a;
    }
    Poly Pow(Poly a, int b1, int b2) { // b1 = b % P, b2 = b % phi(P) and b >= n if a[0] > 0
        int n = a.size(), d = 0, k;
        while (d < n && !a[d]) {
            ++d;
        }
        if ((i64)d * b1 >= n) {
            return Poly(n);
        }
        a.erase(a.begin(), a.begin() + d);
        k = qpow(a[0]), norm(a *= k);
        a = Pow1(a, b1) * qpow(k, P - 1 - b2);
        a.resize(n), d *= b1;
        fd(i, n - 1, 0) {
            a[i] = i >= d ? a[i - d] : 0;
        }
        return a;
    }

    Poly Sin(Poly &a) {
        int i = qpow(3, (P - 1) / 4);
        Poly x(a * i);
        return (Exp(x) - Exp((P - 1) * x)) * qpow(2 * i % P);
    }

    Poly Cos(Poly &a) {
        int i = qpow(3, (P - 1) / 4);
        Poly x(a * i);
        return (Exp(x) + Exp((P - 1) * x)) * qpow(2);
    }

    Poly ASin(Poly &a) {
        int i = qpow(3, (P - 1) / 4);
        return (P - 1) * i % P * Ln(i * a + Sqrt(Poly{1} - a * a));
    }

    Poly ATan(Poly &a) {
        int i = qpow(3, (P - 1) / 4);
        return i * qpow(2) % P * (Ln(Poly{1} - i * a) - Ln(Poly{1} + i * a));
    }

    // Get [x ^ k](f / g)
    int divAt(Poly f, Poly g, i64 k) {
        int n = std::max(f.size(), g.size()), m = norm(n);
        for (; k; k >>= 1) {
            f.resize(m * 2, 0), DFT(f);
            g.resize(m * 2, 0), DFT(g);
            fp(i, 0, 2 * m - 1) {
                f[i] = MUL(f[i], g[i ^ 1]);
            }
            fp(i, 0, m - 1) {
                g[i] = MUL(g[2 * i], g[2 * i + 1]);
            }
            g.resize(m);
            IDFT(f), IDFT(g);
            for (int i = 0, j = k & 1; i < n; i++, j += 2) {
                f[i] = f[j];
            }
            f.resize(n), g.resize(n);
        }
        return f[0];
    }

    // Get a[k] by a[n] = sum c[i] * a[n - i]
    int LinearRecur(Poly a, Poly c, i64 k) {
        c[0] = P - 1, a = a * c;
        a.resize(c.size() - 1);
        return divAt(a, c, k);
    }

    //Binary convolution for &^|
    Poly operator|(Poly a, Poly b) {
        int n = std::max(a.size(), b.size()), N = norm(n);
        a.resize(N);
        FWT::FWTor(a, false);

        b.resize(N);
        FWT::FWTor(b, false);

        Poly A(N);
        fp(i, 0, N - 1) {
            A[i] = MUL(a[i], b[i]);
        }
        FWT::FWTor(A, true);
        return A;
    }

    Poly operator&(Poly a, Poly b) {
        int n = std::max(a.size(), b.size()), N = norm(n);
        a.resize(N);
        FWT::FWTand(a, false);

        b.resize(N);
        FWT::FWTand(b, false);

        Poly A(N);
        fp(i, 0, N - 1) {
            A[i] = MUL(a[i], b[i]);
        }

        FWT::FWTand(A, true);
        return A;
    }

    Poly operator^(Poly a, Poly b) {
        int n = std::max(a.size(), b.size()), N = norm(n);
        a.resize(N);
        FWT::FWTxor(a, false);

        b.resize(N);
        FWT::FWTxor(b, false);

        Poly A(N);
        fp(i, 0, N - 1) {
            A[i] = MUL(a[i], b[i]);
        }

        FWT::FWTxor(A, true);
        return A;
    }
    
}

using namespace Polynomial;

struct Comb {
    int n;
    std::vector<int> fac, invfac, inv;

    Comb() : n{0}, fac{1}, invfac{1}, inv{0} {}
    Comb(int n) : n{n} {
        init();
    }

    void init() {
        fac.resize(n + 1);
        invfac.resize(n + 1);
        inv.resize(n + 1);

        inv[1] = 1;
        for (int i = 2; i <= n; i++) {
            inv[i] = P - 1ll * inv[P % i] * (P / i) % P;
        }

        fac[0] = invfac[0] = 1;
        for (int i = 1; i <= n; i++) {
            fac[i] = 1ll * fac[i - 1] * i % P;
        }
        invfac[n - 1] = qpow(fac[n - 1], P - 2);
        for (int i = n - 1; i; i--) {
            invfac[i - 1] = 1ll * invfac[i] * i % P;
        }
    }
    
    int binom(int n, int m) {
        if (n < m || m < 0 || n < 0) {
            return 0;
        }
        return 1ll * fac[n] * invfac[m] % P * invfac[n - m] % P;
    }
} comb;

```


```c++
template<class T>
#define constexpr
constexpr T power(T a, i64 b) {
    T res = 1;
    for (; b; b /= 2, a *= a) {
        if (b % 2) {
            res *= a;
        }
    }
    return res;
}
 
template<int P>
struct MInt {
    int x;
    constexpr MInt() : x{} {}
    constexpr MInt(i64 x) : x{norm(x % getMod())} {}
     
    static int Mod;
    constexpr static int getMod() {
        if (P > 0) {
            return P;
        } else {
            return Mod;
        }
    }
    constexpr static void setMod(int Mod_) {
        Mod = Mod_;
    }
    constexpr int norm(int x) const {
        if (x < 0) {
            x += getMod();
        }
        if (x >= getMod()) {
            x -= getMod();
        }
        return x;
    }
    constexpr int val() const {
        return x;
    }
    explicit constexpr operator int() const {
        return x;
    }
    constexpr MInt operator-() const {
        MInt res;
        res.x = norm(getMod() - x);
        return res;
    }
    constexpr MInt inv() const {
        assert(x != 0);
        return power(*this, getMod() - 2);
    }
    constexpr MInt &operator*=(MInt rhs) & {
        x = 1LL * x * rhs.x % getMod();
        return *this;
    }
    constexpr MInt &operator+=(MInt rhs) & {
        x = norm(x + rhs.x);
        return *this;
    }
    constexpr MInt &operator-=(MInt rhs) & {
        x = norm(x - rhs.x);
        return *this;
    }
    constexpr MInt &operator/=(MInt rhs) & {
        return *this *= rhs.inv();
    }
    friend constexpr MInt operator*(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res *= rhs;
        return res;
    }
    friend constexpr MInt operator+(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res += rhs;
        return res;
    }
    friend constexpr MInt operator-(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res -= rhs;
        return res;
    }
    friend constexpr MInt operator/(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res /= rhs;
        return res;
    }
    friend constexpr std::istream &operator>>(std::istream &is, MInt &a) {
        i64 v;
        is >> v;
        a = MInt(v);
        return is;
    }
    friend constexpr std::ostream &operator<<(std::ostream &os, const MInt &a) {
        return os << a.val();
    }
    friend constexpr bool operator==(MInt lhs, MInt rhs) {
        return lhs.val() == rhs.val();
    }
    friend constexpr bool operator!=(MInt lhs, MInt rhs) {
        return lhs.val() != rhs.val();
    }
};
 
template<>
int MInt<0>::Mod = 1;
 
template<int V, int P>
constexpr MInt<P> CInv = MInt<P>(V).inv();
 
const int P = 998244353;
using Z = MInt<P>;
 
std::vector<int> rev;
template<int P>
std::vector<MInt<P>> roots{0, 1};
 
template<int P>
constexpr MInt<P> findPrimitiveRoot() {
    MInt<P> i = 2;
    int k = __builtin_ctz(P - 1);
    while (true) {
        if (power(i, (P - 1) / 2) != 1) {
            break;
        }
        i += 1;
    }
    return power(i, (P - 1) >> k);
}
 
template<int P>
constexpr MInt<P> primitiveRoot = findPrimitiveRoot<P>();
 
template<>
constexpr MInt<998244353> primitiveRoot<998244353> {31};
 
template<int P>
constexpr void dft(std::vector<MInt<P>> &a) {
    int n = a.size();
     
    if (int(rev.size()) != n) {
        int k = __builtin_ctz(n) - 1;
        rev.resize(n);
        for (int i = 0; i < n; i++) {
            rev[i] = rev[i >> 1] >> 1 | (i & 1) << k;
        }
    }
     
    for (int i = 0; i < n; i++) {
        if (rev[i] < i) {
            std::swap(a[i], a[rev[i]]);
        }
    }
    if (roots<P>.size() < n) {
        int k = __builtin_ctz(roots<P>.size());
        roots<P>.resize(n);
        while ((1 << k) < n) {
            auto e = power(primitiveRoot<P>, 1 << (__builtin_ctz(P - 1) - k - 1));
            for (int i = 1 << (k - 1); i < (1 << k); i++) {
                roots<P>[2 * i] = roots<P>[i];
                roots<P>[2 * i + 1] = roots<P>[i] * e;
            }
            k++;
        }
    }
    for (int k = 1; k < n; k *= 2) {
        for (int i = 0; i < n; i += 2 * k) {
            for (int j = 0; j < k; j++) {
                MInt<P> u = a[i + j];
                MInt<P> v = a[i + j + k] * roots<P>[k + j];
                a[i + j] = u + v;
                a[i + j + k] = u - v;
            }
        }
    }
}
 
template<int P>
constexpr void idft(std::vector<MInt<P>> &a) {
    int n = a.size();
    std::reverse(a.begin() + 1, a.end());
    dft(a);
    MInt<P> inv = (1 - P) / n;
    for (int i = 0; i < n; i++) {
        a[i] *= inv;
    }
}
 
template<int P = 998244353>
struct Poly : public std::vector<MInt<P>> {
    using Value = MInt<P>;
     
    Poly() : std::vector<Value>() {}
    explicit constexpr Poly(int n) : std::vector<Value>(n) {}
     
    explicit constexpr Poly(const std::vector<Value> &a) : std::vector<Value>(a) {}
    constexpr Poly(const std::initializer_list<Value> &a) : std::vector<Value>(a) {}
     
    template<class InputIt, class = std::_RequireInputIter<InputIt>>
    explicit constexpr Poly(InputIt first, InputIt last) : std::vector<Value>(first, last) {}
     
    template<class F>
    explicit constexpr Poly(int n, F f) : std::vector<Value>(n) {
        for (int i = 0; i < n; i++) {
            (*this)[i] = f(i);
        }
    }
     
    constexpr Poly shift(int k) const {
        if (k >= 0) {
            auto b = *this;
            b.insert(b.begin(), k, 0);
            return b;
        } else if (this->size() <= -k) {
            return Poly();
        } else {
            return Poly(this->begin() + (-k), this->end());
        }
    }
    constexpr Poly trunc(int k) const {
        Poly f = *this;
        f.resize(k);
        return f;
    }
    constexpr friend Poly operator+(const Poly &a, const Poly &b) {
        Poly res(std::max(a.size(), b.size()));
        for (int i = 0; i < a.size(); i++) {
            res[i] += a[i];
        }
        for (int i = 0; i < b.size(); i++) {
            res[i] += b[i];
        }
        return res;
    }
    constexpr friend Poly operator-(const Poly &a, const Poly &b) {
        Poly res(std::max(a.size(), b.size()));
        for (int i = 0; i < a.size(); i++) {
            res[i] += a[i];
        }
        for (int i = 0; i < b.size(); i++) {
            res[i] -= b[i];
        }
        return res;
    }
    constexpr friend Poly operator-(const Poly &a) {
        std::vector<Value> res(a.size());
        for (int i = 0; i < int(res.size()); i++) {
            res[i] = -a[i];
        }
        return Poly(res);
    }
    constexpr friend Poly operator*(Poly a, Poly b) {
        if (a.size() == 0 || b.size() == 0) {
            return Poly();
        }
        if (a.size() < b.size()) {
            std::swap(a, b);
        }
        int n = 1, tot = a.size() + b.size() - 1;
        while (n < tot) {
            n *= 2;
        }
        if (((P - 1) & (n - 1)) != 0 || b.size() < 128) {
            Poly c(a.size() + b.size() - 1);
            for (int i = 0; i < a.size(); i++) {
                for (int j = 0; j < b.size(); j++) {
                    c[i + j] += a[i] * b[j];
                }
            }
            return c;
        }
        a.resize(n);
        b.resize(n);
        dft(a);
        dft(b);
        for (int i = 0; i < n; ++i) {
            a[i] *= b[i];
        }
        idft(a);
        a.resize(tot);
        return a;
    }
    constexpr friend Poly operator*(Value a, Poly b) {
        for (int i = 0; i < int(b.size()); i++) {
            b[i] *= a;
        }
        return b;
    }
    constexpr friend Poly operator*(Poly a, Value b) {
        for (int i = 0; i < int(a.size()); i++) {
            a[i] *= b;
        }
        return a;
    }
    constexpr friend Poly operator/(Poly a, Value b) {
        for (int i = 0; i < int(a.size()); i++) {
            a[i] /= b;
        }
        return a;
    }
    constexpr Poly &operator+=(Poly b) {
        return (*this) = (*this) + b;
    }
    constexpr Poly &operator-=(Poly b) {
        return (*this) = (*this) - b;
    }
    constexpr Poly &operator*=(Poly b) {
        return (*this) = (*this) * b;
    }
    constexpr Poly &operator*=(Value b) {
        return (*this) = (*this) * b;
    }
    constexpr Poly &operator/=(Value b) {
        return (*this) = (*this) / b;
    }
    constexpr Poly deriv() const {
        if (this->empty()) {
            return Poly();
        }
        Poly res(this->size() - 1);
        for (int i = 0; i < this->size() - 1; ++i) {
            res[i] = (i + 1) * (*this)[i + 1];
        }
        return res;
    }
    constexpr Poly integr() const {
        Poly res(this->size() + 1);
        for (int i = 0; i < this->size(); ++i) {
            res[i + 1] = (*this)[i] / (i + 1);
        }
        return res;
    }
    constexpr Poly inv(int m) const {
        Poly x{(*this)[0].inv()};
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (Poly{2} - trunc(k) * x)).trunc(k);
        }
        return x.trunc(m);
    }
    constexpr Poly log(int m) const {
        return (deriv() * inv(m)).integr().trunc(m);
    }
    constexpr Poly exp(int m) const {
        Poly x{1};
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (Poly{1} - x.log(k) + trunc(k))).trunc(k);
        }
        return x.trunc(m);
    }
    constexpr Poly pow(int k, int m) const {
        int i = 0;
        while (i < this->size() && (*this)[i] == 0) {
            i++;
        }
        if (i == this->size() || 1LL * i * k >= m) {
            return Poly(m);
        }
        Value v = (*this)[i];
        auto f = shift(-i) * v.inv();
        return (f.log(m - i * k) * k).exp(m - i * k).shift(i * k) * power(v, k);
    }
    constexpr Poly sqrt(int m) const {
        Poly x{1};
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x + (trunc(k) * x.inv(k)).trunc(k)) * CInv<2, P>;
        }
        return x.trunc(m);
    }
    constexpr Poly mulT(Poly b) const {
        if (b.size() == 0) {
            return Poly();
        }
        int n = b.size();
        std::reverse(b.begin(), b.end());
        return ((*this) * b).shift(-(n - 1));
    }
    constexpr std::vector<Value> eval(std::vector<Value> x) const {
        if (this->size() == 0) {
            return std::vector<Value>(x.size(), 0);
        }
        const int n = std::max(x.size(), this->size());
        std::vector<Poly> q(4 * n);
        std::vector<Value> ans(x.size());
        x.resize(n);
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                q[p] = Poly{1, -x[l]};
            } else {
                int m = (l + r) / 2;
                build(2 * p, l, m);
                build(2 * p + 1, m, r);
                q[p] = q[2 * p] * q[2 * p + 1];
            }
        };
        build(1, 0, n);
        std::function<void(int, int, int, const Poly &)> work = [&](int p, int l, int r, const Poly &num) {
            if (r - l == 1) {
                if (l < int(ans.size())) {
                    ans[l] = num[0];
                }
            } else {
                int m = (l + r) / 2;
                work(2 * p, l, m, num.mulT(q[2 * p + 1]).resize(m - l));
                work(2 * p + 1, m, r, num.mulT(q[2 * p]).resize(r - m));
            }
        };
        work(1, 0, n, mulT(q[1].inv(n)));
        return ans;
    }
};
 
template<int P = 998244353>
Poly<P> berlekampMassey(const Poly<P> &s) {
    Poly<P> c;
    Poly<P> oldC;
    int f = -1;
    for (int i = 0; i < s.size(); i++) {
        auto delta = s[i];
        for (int j = 1; j <= c.size(); j++) {
            delta -= c[j - 1] * s[i - j];
        }
        if (delta == 0) {
            continue;
        }
        if (f == -1) {
            c.resize(i + 1);
            f = i;
        } else {
            auto d = oldC;
            d *= -1;
            d.insert(d.begin(), 1);
            MInt<P> df1 = 0;
            for (int j = 1; j <= d.size(); j++) {
                df1 += d[j - 1] * s[f + 1 - j];
            }
            assert(df1 != 0);
            auto coef = delta / df1;
            d *= coef;
            Poly<P> zeros(i - f - 1);
            zeros.insert(zeros.end(), d.begin(), d.end());
            d = zeros;
            auto temp = c;
            c += d;
            if (i - temp.size() > f - oldC.size()) {
                oldC = temp;
                f = i;
            }
        }
    }
    c *= -1;
    c.insert(c.begin(), 1);
    return c;
}
 
 
template<int P = 998244353>
MInt<P> linearRecurrence(Poly<P> p, Poly<P> q, i64 n) {
    int m = q.size() - 1;
    while (n > 0) {
        auto newq = q;
        for (int i = 1; i <= m; i += 2) {
            newq[i] *= -1;
        }
        auto newp = p * newq;
        newq = q * newq;
        for (int i = 0; i < m; i++) {
            p[i] = newp[i * 2 + n % 2];
        }
        for (int i = 0; i <= m; i++) {
            q[i] = newq[i * 2];
        }
        n /= 2;
    }
    return p[0] / q[0];
}

struct Comb {
    int n;
    std::vector<Z> _fac;
    std::vector<Z> _invfac;
    std::vector<Z> _inv;
    
    Comb() : n{0}, _fac{1}, _invfac{1}, _inv{0} {}
    Comb(int n) : Comb() {
        init(n);
    }
    
    void init(int m) {
        if (m <= n) return;
        _fac.resize(m + 1);
        _invfac.resize(m + 1);
        _inv.resize(m + 1);
        
        for (int i = n + 1; i <= m; i++) {
            _fac[i] = _fac[i - 1] * i;
        }
        _invfac[m] = _fac[m].inv();
        for (int i = m; i > n; i--) {
            _invfac[i - 1] = _invfac[i] * i;
            _inv[i] = _invfac[i] * _fac[i - 1];
        }
        n = m;
    }
    
    Z fac(int m) {
        if (m > n) init(2 * m);
        return _fac[m];
    }
    Z invfac(int m) {
        if (m > n) init(2 * m);
        return _invfac[m];
    }
    Z inv(int m) {
        if (m > n) init(2 * m);
        return _inv[m];
    }
    Z binom(int n, int m) {
        if (n < m || m < 0) return 0;
        return fac(n) * invfac(m) * invfac(n - m);
    }
} comb;
```
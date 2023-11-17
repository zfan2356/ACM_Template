```c++
template<class T = double>
struct __Complex {
    T x, y;
    __Complex() = default;
    __Complex(const T& x, const T& y) : x(x), y(y) {}
    __Complex& operator+=(const __Complex& b) {
        x += b.x;
        y += b.y;
        return *this;
    }
    __Complex& operator-=(const __Complex& b) {
        x -= b.x;
        y -= b.y;
        return *this;
    }
    __Complex& operator*=(const __Complex& b) {
        __Complex temp = *this;
        temp.x = x * b.x - y * b.y;
        temp.y = x * b.y + y * b.x;
        *this = temp;
        return *this;
    }
    __Complex& operator*=(const double& b) {
        this -> x *= b;
        this -> y *=b;
        return *this;
    }
    __Complex operator+(const __Complex& b) {
        __Complex a = *this;
        a += b;
        return a;
    }
    __Complex operator-(const __Complex& b) {
        __Complex a = *this;
        a -= b;
        return a;
    }
    __Complex operator*(const __Complex& b) {
        __Complex a = *this;
        a *= b;
        return a;
    }
    __Complex conj() {
        return __Complex(x, -y);
    }
    friend std::ostream& operator<<(std::ostream& os, const __Complex& a) {
        os << a.x << " " << a.y;
        return os;
    }
};
using Complex = __Complex<>;
const long double PI = acos(-1.0);
const long double PI2 = PI / 2;
const int CONQUER_BIT = 16;
const int CONQUER_MASK = (1 << CONQUER_BIT) - 1;
const Complex I(0, 1);
std::vector<Complex> r;
int preLg;
void pre(const int lg) {
    r.resize(1 << lg);
    for (int i = preLg ; i < lg ; i++) {
        int L = 1 << i;
        r[L] = Complex(cos(PI2 / L), sin(PI2 / L));
        for (int j = L + 1 ; j < (L << 1) ; j++) {
            r[j] = r[j - L] * r[L];
        }
    }
}
template<class T> long long Tint2(const T val, const int mod) {
    long long v = val;
    return ((v < 0) ? (mod + (v - 1) / 2 % mod) : (v + 1) / 2) % mod;
}
template<class T> struct Poly {
    std::vector<T> a;
    Poly(const int size) {
        a.resize(size);
    }
    T &operator[](const int x) {
        return a[x];
    }
    void resize(const int n) {
        a.resize(n);
    }
    int size() {
        return a.size();
    }
    void FFT() {
        int n = a.size();
        for (int i = n ; i >= 2 ; i >>= 1) {
            int L = i >> 1;
            for (int j = 0 ; j != L ; j++) {
                Complex x = a[j], y = a[j + L];
                a[j] = x + y;
                a[j + L] = x - y;
            }
            for (int j = i, m = 1 ; j != n ; j += i, m++) {
                Complex rt = r[m];
                for (int k = 0 ; k != L ; k++) {
                    Complex x = a[j + k], y = a[j + k + L] * rt;
                    a[j + k] = x + y;
                    a[j + k + L] = x - y;
                }
            }
        }
    }
    void IFFT() {
        int n = a.size();
        for (int i = 2 ; i <= n ; i <<= 1) {
            int L = i >> 1;
            for (int j = 0 ; j != L ; j++) {
                Complex x = a[j], y = a[j + L];
                a[j] = x + y;
                a[j + L] = x - y;
            }
            for (int j = i, m = 1 ; j != n ; j += i, m++) {
                Complex rt = r[m];
                for (int k = 0 ; k != L ; k++) {
                    Complex x = a[j + k], y = a[j + k + L];
                    a[j + k] = x + y;
                    a[j + k + L] = (x - y) * rt;
                }
            }
        }
        double inv = 1.0 / n;
        for (int i = 0 ; i < n ; i++) {
            a[i] *= inv;
        }
        reverse(begin(a) + 1, end(a));
    }
    void mul(Poly &x, const int mod) {
        int n = 1, u = a.size(), m = x.size(), lg = 0, len = u + m - 1;
        while (n < len) {
            n <<= 1;
            lg++;
        }
        if (lg > preLg) {
            pre(lg);
            preLg = lg;
        }
        Poly<Complex> P(n), Q(n);
        for (int i = 0 ; i < u ; i++) {
            P[i] = Complex(a[i] & CONQUER_MASK, a[i] >> CONQUER_BIT);
        }
        P.FFT();
        Poly<Complex> _Q(P);
        for (int i = 0 ; i < m ; i++) {
            Q[i] = Complex(x[i] & CONQUER_MASK, x[i] >> CONQUER_BIT);
        }
        Q.FFT();
        P[0] *= Q[0].x * 2;
        _Q[0] *= Q[0].y * 2;
        for (int d = 0 ; d < lg ; d++) {
            int L = 1 << d, msk = L - 1;
            for (int i = L ; i < (L << 1) ; i++) {
                Complex &p = Q[i], q = Q[i ^ msk].conj();
                Complex a = (p + q), b = (q - p) * I;
                P[i] *= a;
                _Q[i] *= b;
            }
        }
        P.IFFT();
        _Q.IFFT();
        resize(len);
        for (int i = 0 ; i < len ; i++) {
            long long cur = (Tint2(_Q[i].y, mod) << (CONQUER_BIT << 1))
                            + (Tint2(_Q[i].x + P[i].y, mod) << CONQUER_BIT)
                            + (Tint2(P[i].x, mod));
            a[i] = cur % mod;
        }
    }
}; // Poly
```
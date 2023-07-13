```c++
template<typename T>
struct ExGcd {
    T a, b, x, y, r;

    ExGcd() {}
    ExGcd(T _a, T _b) {
        init(_a, _b);
    }

    void init(T _a, T _b) {
        this->a = _a;
        this->b = _b;
        this->r = exgcd(_a, _b, this->x, this->y);
    }

    T exgcd(T a, T b, T &x, T &y){	// 扩展欧几里得算法, 求x, y，使得ax + by = gcd(a, b)
        if(!b) {
            return x = 1, y = 0, a;
        }
        T r = exgcd(b, a % b, x, y);
        tie(x, y) = make_tuple(y, x - (a / b) * y);
        return r;
    }

    T getInv(T _a, T _b, T _m) {   //ax = b (mod m)中x的值
        init(_a, _m);
        if (_b % r != 0) {
            return -1;
        } else {
            return (_b / r * x % _m + _m) % _m;
        }
    }
};

template<typename T>
T MUL(T a, T b, T mod) {
    T ans = 0;
    for (; b; a = (a + a) % mod, b >>= 1) {
        if (b & 1) {
            ans = (ans + a) % mod;
        }
    }
    return ans;
}

template<typename T>
T ExCRT(int N, std::vector<T>& A, std::vector<T>& M) {   // x = A_i (mod M_i)
    auto MOD = [&](T a, T b) {
        return ((a % b) + b) % b;
    };

    for (int i = 0; i < N; i++) {
        ExGcd<T> e(M[0], M[i]);
        T c = MOD(A[i] - A[0], M[i]);
        if (c % e.r != 0) {
            return -1;
        }

//        e.x = MOD(e.x * c / e.r, abs(M[i] / e.r));
        e.x = MUL(e.x, c / e.r, abs(M[i] / e.r));
        A[0] = e.x * M[0] + A[0];
        M[0] = abs(M[0] / e.r * M[i]);
        A[0] = MOD(A[0], M[0]);
    }
    return A[0];
}

```
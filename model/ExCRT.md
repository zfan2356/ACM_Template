```c++
template<typename T>
T exgcd(T a, T b, T &x, T &y) {
    if(!b) {
        x = 1, y = 0;
        return a;
    }
    T r = exgcd(b, a % b, x, y);
    tie(x, y) = make_tuple(y, x - (a / b) * y);
    return r;
}

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
        T x, y;
        T r = exgcd(M[0], M[i], x, y);
        T c = MOD(A[i] - A[0], M[i]);
        if (c % r != 0) {
            return -1;
        }

        x = MOD(x * c / r, abs(M[i] / r));
//        x = MUL(x, c / r, abs(M[i] / r));
        A[0] = x * M[0] + A[0];
        M[0] = abs(M[0] / r * M[i]);
        A[0] = MOD(A[0], M[0]);
    }
    return A[0];
}

```
```c++
const int mod = 1e9 + 7;
template<typename T>
struct Matrix {
    int m, n;
    vector<vector<T>> data;

    Matrix() : m(0), n(0) {
        data.resize(0);
    }

    Matrix(int row, int col, T num = 0) : m(row), n(col) {
        data.resize(m, vector<T>(n, num));
    }

    Matrix(int n) : m(n), n(n) {
        data.resize(n, vector<T>(n, T()));
        for (int i = 0; i < n; i++) {
            data[i][i] = 1;
        }
    }

    Matrix operator+(const Matrix &b) const {
        if (m != b.m || n != b.n) {
            return Matrix();
        }
        Matrix ans(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans.data[i][j] = (data[i][j] + b.data[i][j]) % mod;
            }
        }
        return ans;
    }

    Matrix operator-(const Matrix &b) const {
        if (m != b.m || n != b.n) {
            return Matrix();
        }
        Matrix ans(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans.data[i][j] = (data[i][j] - b.data[i][j]) % mod;
            }
        }
        return ans;
    }

    Matrix operator*(const T &b) const {
        Matrix ans(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans.data[i][j] = (data[i][j] * b) % mod;
            }
        }
        return ans;
    }

    Matrix operator*(const Matrix &b) const {
        if (n != b.m) {
            return Matrix();
        }
        Matrix ans(m, b.n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < b.n; j++) {
                for (int k = 0; k < n; k++) {
                    ans.data[i][j] = (ans.data[i][j] + data[i][k] * b.data[k][j] % mod) % mod;
                }
            }
        }
        return ans;
    }

    Matrix operator%(const T &b) const {
        Matrix ans(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans.data[i][j] = data[i][j] % b;
            }
        }
        return ans;
    }

    vector<int> &operator[](int rank) {
        return data[rank];
    }

    T Sum() const {
        T ans = 0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans = (ans + data[i][j]) % mod;
            }
        }
        return ans;
    }

    T Mul() const {
        T ans = 1;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans = ans * data[i][j] % mod;
            }
        }
        return ans;
    }


    Matrix Inv() {
        Matrix res(n), a = *this;
        for (int i = 0, r; i < n; ++i) {
            r = i;
            for (int j = i + 1; j < n; ++j) {
                if (a.data[j][i] > a.data[r][i]) {
                    r = j;
                }
            }
            swap(a.data[i], a.data[r]);
            swap(res.data[i], res.data[r]);

            if (!a.data[i][i]) {
                return res.data[0][0] = -1, res;
            }

            T Invaii = qpow(a.data[i][i], mod - 2);
            for (int k = 0; k < n; ++k) {
                a.data[i][k] = a.data[i][k] * Invaii % mod;
            }
            for (int k = 0; k < n; ++k) {
                res.data[i][k] = res.data[i][k] * Invaii % mod;
            }
            for (int j = 0; j < n; ++j)
                if (j != i) {
                    T tmp = a.data[j][i];
                    for (int k = i; k < n; ++k) {
                        a.data[j][k] = (a.data[j][k] - tmp * a.data[i][k] % mod + mod) % mod;
                    }
                    for (int k = 0; k < n; ++k) {
                        res.data[j][k] = (res.data[j][k] - tmp * res.data[i][k] % mod + mod) % mod;
                    }
                }
        }
        return res;
    }

    friend Matrix qpow(Matrix a, int b) {
        if (a.m != a.n) {
            return Matrix();
        }

        Matrix ans(a.n);
        for (; b; a = a * a, b >>= 1) {
            if (b & 1) {
                ans = ans * a;
            }
        }
        return ans;
    }

    friend istream &operator>>(istream &in, Matrix &x) {
        in >> x.m >> x.n;
        x.data.resize(x.m, vector<int>(x.n));
        for (int i = 0; i < x.m; i++) {
            for (int j = 0; j < x.n; j++) {
                in >> x.data[i][j];
            }
        }
        return in;
    }

    friend ostream &operator<<(ostream &out, const Matrix &x) {
        for (int i = 0; i < x.m; i++) {
            for (int j = 0; j < x.n; j++) {
                out << x.data[i][j] << (j == x.n - 1 ? '\n' : ' ');
            }
        }
        return out;
    }

    bool operator==(const Matrix &b) {
        if (m != b.m || n != b.n) {
            return false;
        }

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (data[i][j] != b.data[i][j]) return false;
            }
        }
        return true;
    }

    bool operator!=(const Matrix &b) {
        if (m != b.m || n != b.n) {
            return true;
        }

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (data[i][j] != b.data[i][j]) return true;
            }
        }
        return false;
    }
};

```
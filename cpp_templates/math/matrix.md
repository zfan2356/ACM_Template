```c++
template<typename T>
struct Matrix {
    int m, n;
    std::vector<std::vector<T>> data;

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

    Matrix(const Matrix& other) : m(other.m), n(other.n), data(other.data) {}

    Matrix operator+(const Matrix &b) const {
        if (m != b.m || n != b.n) {
            return Matrix();
        }
        Matrix ans(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans.data[i][j] = data[i][j] + b.data[i][j];
            }
        }
        return ans;
    }

    Matrix operator=(const Matrix& other) {
        if (this == &other) {
            return *this;
        }

        m = other.m;
        n = other.n;
        data = other.data;

        return *this;
    }


    Matrix operator-(const Matrix &b) const {
        if (m != b.m || n != b.n) {
            return Matrix();
        }
        Matrix ans(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans.data[i][j] = data[i][j] - b.data[i][j];
            }
        }
        return ans;
    }

    Matrix operator*(const T &b) const {
        Matrix ans(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans.data[i][j] = data[i][j] * b;
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
                    ans.data[i][j] += data[i][k] * b.data[k][j];
                }
            }
        }
        return ans;
    }

    std::vector<T> &operator[](int rank) {
        return data[rank];
    }

    T Sum() const {
        T ans = 0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans += data[i][j];
            }
        }
        return ans;
    }

    T Mul() const {
        T ans = 1;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ans *= data[i][j];
            }
        }
        return ans;
    }

    friend Matrix qpow(Matrix &a, T b) {
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

    friend std::istream &operator>>(std::istream &in, Matrix &x) {
        in >> x.m >> x.n;
        x.data.resize(x.m, vector<int>(x.n));
        for (int i = 0; i < x.m; i++) {
            for (int j = 0; j < x.n; j++) {
                in >> x.data[i][j];
            }
        }
        return in;
    }

    friend std::ostream &operator<<(std::ostream &out, const Matrix &x) {
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

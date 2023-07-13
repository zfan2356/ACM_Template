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

```
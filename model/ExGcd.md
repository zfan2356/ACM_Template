```c++
template<typename T>
T exgcd(T a, T b, T &x, T &y) {
    if(!b) {
        x = 1, y = 0;
        return a;
    }
    T r = exgcd(b, a % b, x, y);
    std::tie(x, y) = std::make_tuple(y, x - (a / b) * y);
    return r;
}

auto getInv = [&](int a, int b, int m) { //ax = b (mod m)中x的值
    int x, y;
    int r = exgcd(a, m, x, y);
    if (b % r != 0) {
        return -1;
    }
    return (b / r * x % m + m) % m;
};

```
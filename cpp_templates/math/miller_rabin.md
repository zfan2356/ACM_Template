```c++
namespace Miller_Rabin{
//    typedef __int128 i64;		//如果范围较大可以使用__int128

    __int128 read() {
        __int128 x = 0, f = 1;
        char ch = getchar();
        while (!isdigit(ch) && ch != '-')ch = getchar();
        if (ch == '-')f = -1, ch = getchar();
        while (isdigit(ch))x = x * 10 + ch - '0', ch = getchar();
        return f * x;
    }

    void print(__int128 x) {
        if (x < 0) {
            putchar('-');
            x = -x;
        }
        if (x > 9)
            print(x / 10);
        putchar(x % 10 + '0');
    }
    i64 gcd(i64 a, i64 b) {
        return b ? gcd(b, a % b) : a;
    }
    i64 mul(i64 a, i64 b, i64 mod) {		//龟速乘
        i64 ret = 0;
        while (b) {
            if (b & 1) ret = (ret + a) % mod;
            a = (a + a) % mod;
            b >>= 1;
        }
        return ret;
    }
    i64 pow(i64 a, i64 b, i64 mod) {		//快速幂
        i64 ret = 1;
        while (b) {
            if (b & 1) ret = mul(ret, a, mod);
            a = mul(a, a, mod);
            b >>= 1;
        }
        return ret;
    }
    bool check(i64 a, i64 n) {
        i64 x = n - 1, y;
        int t = 0;
        while ((x & 1) == 0) {
            x >>= 1;
            t++;
        }
        x = pow(a, x, n);
        for (int i = 1; i <= t; i++) {
            y = mul(x, x, n);
            if (y == 1 && x != 1 && x != n - 1) return true;
            x = y;
        }
        if (y != 1) return true;
        return false;
    }
    bool isprime(i64 n) {
        if (n == 2) return true;
        if (n == 1 || !(n & 1)) return false;
        const int arr[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
        for (int i = 0; i < 12; i++) {
            if (arr[i] >= n) break;
            if (check(arr[i], n)) return false;
        }
        return true;
    }
}
```
